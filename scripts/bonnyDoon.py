#!/usr/bin/env python
"""Run signal-to-reference alignments
"""
from __future__ import print_function
import pandas as pd

from signalAlignLib import *
from alignmentAnalysisLib import CallMethylation, get_first_sequence
from variantCallingLib import scan_for_proposals
from multiprocessing import Process, Queue, current_process, Manager
from serviceCourse.file_handlers import FolderHandler
from argparse import ArgumentParser
from random import shuffle
from shutil import copyfile
from operator import itemgetter

STEP = 5


def parse_args():
    parser = ArgumentParser(description=__doc__)

    parser.add_argument('--file_directory', '-d', action='store',
                        dest='files_dir', required=True, type=str, default=None,
                        help="directory with MinION fast5 reads to align")
    parser.add_argument('--ref', '-r', action='store',
                        dest='ref', required=True, type=str,
                        help="reference sequence to align to, in FASTA")
    parser.add_argument('--in_template_hmm', '-T', action='store', dest='in_T_Hmm',
                        required=False, type=str, default=None,
                        help="input HMM for template events, if you don't want the default")
    parser.add_argument('--in_complement_hmm', '-C', action='store', dest='in_C_Hmm',
                        required=False, type=str, default=None,
                        help="input HMM for complement events, if you don't want the default")
    parser.add_argument('--templateHDP', '-tH', action='store', dest='templateHDP', default=None,
                        help="template serialized HDP file")
    parser.add_argument('--complementHDP', '-cH', action='store', dest='complementHDP', default=None,
                        help="complement serialized HDP file")
    parser.add_argument('--degenerate', '-x', action='store', dest='degenerate', default="variant",
                        help="Specify degenerate nucleotide options: "
                             "variant -> {ACGT}, twoWay -> {CE} threeWay -> {CEO}")
    parser.add_argument('--stateMachineType', '-smt', action='store', dest='stateMachineType', type=str,
                        default="threeState", help="decide which model to use, threeState by default")
    parser.add_argument('--threshold', '-t', action='store', dest='threshold', type=float, required=False,
                        default=None, help="posterior match probability threshold, Default: 0.01")
    parser.add_argument('--diagonalExpansion', '-e', action='store', dest='diag_expansion', type=int,
                        required=False, default=None, help="number of diagonals to expand around each anchor")
    parser.add_argument('--constraintTrim', '-m', action='store', dest='constraint_trim', type=int,
                        required=False, default=None, help='amount to remove from an anchor constraint')
    parser.add_argument('--target_regions', '-q', action='store', dest='target_regions', type=str,
                        required=False, default=None, help="tab separated table with regions to align to")
    parser.add_argument('---un-banded', '-ub', action='store_false', dest='banded',
                        default=True, help='flag, turn off banding')
    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False,
                        default=4, type=int, help="number of jobs to run concurrently")
    parser.add_argument('--nb_files', '-n', action='store', dest='nb_files', required=False,
                        default=500, type=int, help="maximum number of reads to align")
    # todo help string
    parser.add_argument('--cycles', dest='cycles', default=1, required=False, type=int)

    parser.add_argument('--output_location', '-o', action='store', dest='out',
                        required=True, type=str, default=None,
                        help="directory to put the alignments")
    # todo help string
    parser.add_argument('--corrected', dest='corrected', required=False, default='corrected.fa')

    # todo added by tpesout (from master, with new signalAlign API)
    parser.add_argument("--bwt", action='store', dest="bwt", default=None, required=False,
                        help="path to BWT files. example: ../ref.fasta")

    # todo these are used by process_reference_fasta, but maybe aren't needed for this script
    # parser.add_argument('--ambig_char', '-X', action='store', required=False, default="X", type=str, dest='ambig_char',
    #                     help="Character to substitute at positions, default is 'X'.")
    # parser.add_argument("--motif", action="store", dest="motif_key", default=None)

    args = parser.parse_args()
    return args

def resolvePath(p):
    if p is None:
        return None
    elif p.startswith("/"):
        return p
    else:
        return os.path.abspath(p)


def group_sites_in_window2(sites, window=6):
    def collect_group(start):
        i = start
        g = [sites[start]]
        while sites[i + 1] - sites[i] < window:
            g.append(sites[i + 1])
            i += 1
            if len(sites) <= i + 1:
                break
        return g, i + 1

    sites.sort()
    groups = []
    i = 0
    while i + 1 < len(sites):
        g, i = collect_group(i)
        groups.append(g)
    return groups


def make_degenerate_reference(input_fasta, start, forward_sequence_path, backward_sequence_path,
                              block_size=1, step=6):
    """
    input_sequence: string, input nucleotide sequence
    out_path: string, path to directory to put new sequences with substituted degenerate characters
    block_size: not implemented
    step: number of bases between degenerate characters
    :return (subbed sequence, complement subbed sequence)
    """

    input_sequence = get_first_sequence(input_fasta)

    complement_sequence = reverse_complement(dna=input_sequence, reverse=False, complement=True)

    t_seq = list(input_sequence)
    c_seq = list(complement_sequence)

    positions = xrange(start, len(input_sequence), step)
    for position in positions:
        t_seq[position] = "X"
        c_seq[position] = "X"

    t_seq = ''.join(t_seq)
    c_seq = ''.join(c_seq)

    sequence_length = len(t_seq)

    with open(forward_sequence_path, 'w') as f:
        f.write("{seq}".format(seq=t_seq))
    with open(backward_sequence_path, 'w') as f:
        f.write("{seq}".format(seq=c_seq))

    return True, sequence_length


def aligner(work_queue, done_queue):
    try:
        for f in iter(work_queue.get, 'STOP'):
            alignment = SignalAlignment(**f)
            alignment.run()
    except Exception, e:
        done_queue.put("%s failed with %s" % (current_process().name, e.message))


def run_methyl_caller(work_queue, done_queue):
    try:
        for f in iter(work_queue.get, 'STOP'):
            c = CallMethylation(**f)
            c.write()
    except Exception, e:
        done_queue.put("%s failed with %s" % (current_process().name, e.message))


def load_data(file_path):
    data = pd.read_table(file_path,
                         usecols=(0, 1, 2, 3, 4, 5, 6),
                         names=['site', 'strand', 'pA', 'pC', 'pG', 'pT', 'read'],
                         dtype={'site': np.int64,
                                'strand': np.str,
                                'pC': np.float64,
                                'pmC': np.float64,
                                'phmC': np.float64,
                                'read': np.str,
                                })
    return data


def symbol_to_base(symbol):
    return ["A", "C", "G", "T"][symbol]


def rc_probs(probs):
    return [probs[3], probs[2], probs[1], probs[0]]


def update_reference(data, reference_sequence, min_depth=0, get_sites=False):
    d = load_data(data)

    ref = get_first_sequence(reference_sequence)
    ref = list(ref)

    candidate_sites = []
    add_to_candidates = candidate_sites.append

    for g, x in d.groupby("site"):
        marginal_forward_p = pd.Series(0, ['pA', 'pC', 'pG', 'pT'])
        marginal_backward_p = pd.Series(0, ['pA', 'pC', 'pG', 'pT'])
        assert(len(x['site'].unique()) == 1)
        site = x['site'].unique()[0]

        if len(x['read']) < min_depth:
            continue

        for i, read in x.iterrows():
            if ((read['read'].endswith(".forward.tsv") and read['strand'] == 't') or
                    (read['read'].endswith(".backward.tsv") and read['strand'] == 'c')):
                direction = True
            else:
                direction = False

            if direction:
                marginal_forward_p += read[['pA', 'pC', 'pG', 'pT']]
            else:
                marginal_backward_p += read[['pA', 'pC', 'pG', 'pT']]

        marginal_prob = marginal_forward_p + rc_probs(marginal_backward_p)

        normed_marginal_probs = marginal_prob.map(lambda x: x / sum(marginal_prob))
        called_base = normed_marginal_probs.argmax()[1]
        #called_base = marginal_prob.map(lambda x: x / sum(marginal_prob)).argmax()[1]

        if called_base != ref[site]:
            if get_sites is False:
                print("Changing {orig} to {new} at {site}".format(orig=ref[site], new=called_base, site=site))
                ref[site] = called_base
            else:
                print("Proposing edit at {site} from {orig} to {new}, \n{probs}"
                      "".format(orig=ref[site], new=called_base, site=site, probs=normed_marginal_probs))
                difference = normed_marginal_probs.max() - normed_marginal_probs["p" + ref[site]]
                print(difference)
                add_to_candidates((site, difference))

    if get_sites is True:
        return candidate_sites
    else:
        return ''.join(ref)


def main(args):
    # parse args
    args = parse_args()

    command_line = " ".join(sys.argv[:])
    print("bonnyDoon - Command Line: {cmdLine}\n".format(cmdLine=command_line), file=sys.stderr)

    # get absolute paths to inputs
    args.files_dir           = resolvePath(args.files_dir)
    args.ref                 = resolvePath(args.ref)
    args.out                 = resolvePath(args.out)
    args.bwt                 = resolvePath(args.bwt)
    args.in_T_Hmm            = resolvePath(args.in_T_Hmm)
    args.in_C_Hmm            = resolvePath(args.in_C_Hmm)
    args.templateHDP         = resolvePath(args.templateHDP)
    args.complementHDP       = resolvePath(args.complementHDP)
    args.target_regions      = resolvePath(args.target_regions)

    start_message = """
#   Starting BonnyDoon Error-Correction
#   Aligning files from: {fileDir}
#   Aligning to reference: {reference}
#   Using BWT: {bwt}
#   Aligning maximum of {nbFiles} files
#   Using model: {model}
#   Using banding: {banding}
#   Aligning to regions in: {regions}
#   Non-default template HMM: {inThmm}
#   Non-default complement HMM: {inChmm}
#   Template HDP: {tHdp}
#   Complement HDP: {cHdp}
    """.format(fileDir=args.files_dir, reference=args.ref, bwt=args.bwt, nbFiles=args.nb_files, banding=args.banded,
               inThmm=args.in_T_Hmm, inChmm=args.in_C_Hmm, model=args.stateMachineType, regions=args.target_regions,
               tHdp=args.templateHDP, cHdp=args.complementHDP)

    print(start_message, file=sys.stdout)
    # cull the MinION files
    fast5s = cull_fast5_files(args.files_dir, args.nb_files)

    # get the (input) reference sequence
    if not os.path.isfile(args.ref):
        print("bonnyDoon - Did not find valid reference file", file=sys.stderr)
        sys.exit(1)

    # make a working folder in the specified directory
    temp_folder = FolderHandler()
    temp_dir_path = temp_folder.open_folder(os.path.join(args.out, "tempFiles_errorCorrection"))

    # TODO this is what's in master
    # this is used by signal align, it's a map of reference contigs.  currently this script only supports one contig
    reference_map = process_reference_fasta(fasta=args.ref, motif_key=None, work_folder=temp_folder, sub_char=None)

    # TODO this was old reference mgmt
    # this is used by scan_for_proposals (and modifies
    reference_sequence_path = args.ref
    # unpack the reference sequence
    reference_sequence_string = get_first_sequence(reference_sequence_path)

    # get bwa index
    if args.bwt is not None:
        print("bonnyDoon - using provided BWT %s" % args.bwt)
        bwa_ref_index = args.bwt
    else:
        print("bonnyDoon - indexing reference at %s" % args.ref, file=sys.stderr)
        bwa_ref_index = get_bwa_index(args.ref, temp_dir_path)
        print("bonnyDoon - indexing reference, done", file=sys.stderr)

    # alignment args are the parameters to the HMM/HDP model, and don't change
    alignment_args = {
        "path_to_EC_refs": None,
        "destination": temp_dir_path,
        "stateMachineType": args.stateMachineType,
        "bwa_index": bwa_ref_index,
        "in_templateHmm": args.in_T_Hmm,
        "in_complementHmm": args.in_C_Hmm,
        "in_templateHdp": args.templateHDP,
        "in_complementHdp": args.complementHDP,
        # todo these break signalAlignLib:SignalAlignment ctor
        # "banded": args.banded,
        # "sparse_output": True,
        "threshold": args.threshold,
        "diagonal_expansion": args.diag_expansion,
        "constraint_trim": args.constraint_trim,
        "target_regions": None,
        "degenerate": degenerate_enum(args.degenerate),
    }
    #TODO you could save alignments by altering the above "destination" parameter

    # get the sites that have proposed edits
    print("\n\nbonnyDoon - scanning for proposals with %d fast5s and step %d" % (len(fast5s), STEP))
    output_files = scan_for_proposals(temp_folder, STEP, reference_map, reference_sequence_string, fast5s, alignment_args, args.nb_jobs)
    print("\nbonnyDoon - got {} output files:".format(len(output_files)))
    for output_file in output_files:
        print("\t{}".format(output_file))
    print("\n\nbonnyDoon - fin\n")

    return

if __name__ == "__main__":
    sys.exit(main(sys.argv))
