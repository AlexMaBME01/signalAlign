#!/usr/bin/env python
"""Run signal-to-reference alignments
"""
from __future__ import print_function

from signalAlignLib import *
from alignmentAnalysisLib import CallMethylation, get_first_sequence
from variantCallingLib import make_reference_files_and_alignment_args, run_service, aligner, variant_caller
from serviceCourse.file_handlers import FolderHandler
from argparse import ArgumentParser
import glob

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
    parser.add_argument("--kmer_size", action='store', dest="kmer_size", default=5, required=False,
                        help="size of kmers in fast5 file")

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


def build_fast5_to_read_id_dict(fast5_locations):
    fast5_to_read_id = dict()
    for fast5 in fast5_locations:
        npr = NanoporeRead(fast5, False)
        read_id = npr.read_label
        fast5_id = os.path.basename(fast5)[:-6]
        fast5_to_read_id[fast5_id] = read_id
    return fast5_to_read_id



def discover_single_nucleotide_probabilities(working_folder, step, reference_map, reference_sequence_string,
                                             list_of_fast5s, alignment_args, workers,
                                             output_directory=None, use_saved_alignments=True, save_alignments=True):
    # I'm hacking together the new (improved?) signal align API and the previous version of the api (from when the
    # bonnyDoon script was last working).  The reference map groups by contigs (and is needed by the current SignalAlign
    # API). this script uses the reference sequence string
    # TODO either make this function properly handle the reference_map param, or change SignalAlign to handle single ref string
    if len(reference_map) != 1:
        print("[error] scan_for_proposals must have only one entry in reference map (got %d)" % len(reference_map))
        sys.exit(1)
    reference_map_contig_name = list(reference_map.keys())[0]
    alignment_args['reference_map'] = reference_map
    single_contig_reference_map = reference_map[reference_map_contig_name]

    reference_sequence_length = len(reference_sequence_string)
    assert reference_sequence_length > 0, "Got empty string for reference sequence."

    # proposals will contain the sites that we're going to change to N
    #todo remove this
    proposals = []

    fast5_to_read = build_fast5_to_read_id_dict(list_of_fast5s)
    print("[info] built map of fast5 identifiers to read ids with {} elements".format(len(fast5_to_read)))

    #todo we want to "join" all the steps per fast5

    for s in xrange(step):
        print("\n[info] starting step %d" % s)
        saved_step_dir = os.path.join(working_folder.path, "step_{}".format(s))
        scan_positions = range(s, reference_sequence_length, step)
        #tpesout: changed this function to update the values in single_contig_reference_map to fit new signalAlign API
        check = make_reference_files_and_alignment_args(working_folder, reference_sequence_string,
                                                        single_contig_reference_map, n_positions=scan_positions)
        assert check, "Problem making degenerate reference for step {step}".format(step=s)

        # do or get alignments
        if use_saved_alignments and os.path.isdir(saved_step_dir):
            alignments = [x for x in glob.glob(os.path.join(saved_step_dir, "*.tsv")) if os.stat(x).st_size != 0]
            print("[info] using {} saved alignments in {}".format(len(alignments), saved_step_dir))
        else:
            print("[info] running aligner on %d fast5 files with %d workers" % (len(list_of_fast5s), workers))
            run_service(aligner, list_of_fast5s, alignment_args, workers, "in_fast5")
            alignments = [x for x in glob.glob(os.path.join(working_folder.path, "*.tsv")) if os.stat(x).st_size != 0]
        alignment_count = len(alignments)

        if alignment_count == 0:
            print("[error] Didn't find any alignment files here {}".format(working_folder.path))
            sys.exit(1)
        else:
            print("[info] Found %d alignment files (%d input fast5s) here %s" %
                  (alignment_count, len(list_of_fast5s), working_folder.path))

        marginal_probability_prefix = working_folder.add_file_path("marginals.{step}".format(step=s))

        proposal_args = {
            "sequence": None,
            "out_file_prefix": marginal_probability_prefix,
            # removed to force use of offset and kmer length
            # "positions": {"forward": scan_positions, "backward": scan_positions},
            "step_offset": s,
            "degenerate_type": alignment_args["degenerate"],
            "kmer_length": step
        }

        print("[info] running variant_caller on %d alignments files with %d workers" % (alignment_count, workers))
        run_service(variant_caller, alignments, proposal_args, workers, "alignment_file")

        # remove or save old alignments
        files = glob.glob(working_folder.path + "*.tsv")
        if save_alignments:
            if not os.path.isdir(saved_step_dir): os.mkdir(saved_step_dir)
            print("[info] saving {} alignment files into {}".format(len(files), saved_step_dir))
            for f in files:
                os.rename(f, os.path.join(saved_step_dir, os.path.basename(f)))
        else:
            print("[info] deleting {} alignment files".format(len(files)))
            for f in files:
                os.remove(f)
        print("[info] step %d completed\n" % s)

    # per fast5, we want to coalesce all step calling into one file (named by read)
    if output_directory is None:
        output_directory = os.path.join(working_folder.path, "reads")
    if not os.path.isdir(output_directory): os.mkdir(output_directory)
    print("[info] writing output to {}".format(output_directory))
    output_files = list()
    # iterate over input fast5s
    for fast5_id in fast5_to_read.iterkeys():
        # get files
        files = glob.glob(os.path.join(working_folder.path,
                                       "marginals*{}*{}".format(fast5_id, CallMethylation.FILE_EXTENSION)))
        if len(files) != step:
            print("[error] input fast5 '{}' yielded {} output files, expected {}".format(fast5_id, len(files), step))
            if len(files) == 0:
                continue

        # read all lines in all files
        output_lines = list()
        for file in files:
            with open(file, 'r') as input:
                for line in input:
                    line = line.split("\t")
                    line[0] = int(line[0])
                    output_lines.append(line)
        # sort based on position
        output_lines.sort(key=lambda x: x[0])

        # write output
        output_filename = "{}.tsv".format(fast5_to_read[fast5_id])
        output_file = os.path.join(output_directory, output_filename)
        with open(output_file, 'w') as output:
            output.write("## fast5_input: {}.fast5\n".format(fast5_id))
            output.write("## read_id: {}\n".format(fast5_to_read[fast5_id]))
            output.write("## contig: {}\n".format(reference_map_contig_name))
            output.write("#CHROM\tPOS\tpA\tpC\tpG\tpT\n".format(reference_map_contig_name))
            for line in output_lines:
                line = [reference_map_contig_name, str(line[0]), line[2], line[3], line[4], line[5]]
                output.write("\t".join(line) + "\n")
        #save
        output_files.append(output_file)

    # document and return
    print("[info] wrote {} output files ({} input fast5s) in {}"
          .format(len(output_files), len(fast5_to_read), output_directory))
    return output_files

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
#   Kmer size: {kmerSize}
    """.format(fileDir=args.files_dir, reference=args.ref, bwt=args.bwt, nbFiles=args.nb_files, banding=args.banded,
               inThmm=args.in_T_Hmm, inChmm=args.in_C_Hmm, model=args.stateMachineType, regions=args.target_regions,
               tHdp=args.templateHDP, cHdp=args.complementHDP, kmerSize=args.kmer_size)

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
    print("\n\nsingleNucleotideProbabilities - scanning for proposals with %d fast5s" % len(fast5s))
    output_files = discover_single_nucleotide_probabilities(temp_folder, args.kmer_size, reference_map,
                                                            reference_sequence_string, fast5s, alignment_args,
                                                            args.nb_jobs, output_directory=args.out)
    print("\nsingleNucleotideProbabilities - got {} output files:".format(len(output_files)))
    for output_file in output_files:
        print("\t{}".format(output_file))
    print("\n\nsingleNucleotideProbabilities - fin\n")

    return

if __name__ == "__main__":
    sys.exit(main(sys.argv))
