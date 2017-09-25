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
import shutil


def parse_args():
    parser = ArgumentParser(description=__doc__)

    parser.add_argument('--file_directory', '-d', action='store',
                        dest='files_dir', required=True, type=str, default=None,
                        help="directory with MinION fast5 reads to align")
    parser.add_argument('--ref', '-r', action='store',
                        dest='ref', required=True, type=str,
                        help="reference sequence to align to, in FASTA")
    parser.add_argument('--output_location', '-o', action='store', dest='out',
                        required=True, type=str, default=None,
                        help="directory to put the alignments")
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
    parser.add_argument("--bwt", action='store', dest="bwt", default=None, required=False,
                        help="path to BWT files. example: ../ref.fasta")
    parser.add_argument("--kmer_size", action='store', dest="kmer_size", default=5, required=False,
                        help="size of kmers in fast5 file")
    parser.add_argument("--step_size", action='store', dest="step_size", default=10, required=False,
                        help="distance between positions of uncertainty")

    parser.add_argument("--validate", action='store', dest='validation_file', default=None, required=False,
                        help="validate an output file as compared to its fast5 file (only performs this action)")

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


def validate_snp_directory(snp_directory, reference_sequence_path, print_summary=False, rewrite_files=True):
    # prep
    all_identities = list()
    all_lengths = list()
    all_identity_ratios = list()
    full_reference_sequence = get_first_sequence(reference_sequence_path).upper()
    files = glob.glob(os.path.join(snp_directory, "*.tsv"))
    print("\n[singleNucleotideProbabilities] Validating {} files in {}\n".format(len(files), snp_directory))

    # for rewriting
    if rewrite_files:
        rewrite_files = snp_directory[:-1] if snp_directory.endswith("/") else snp_directory
        rewrite_files += ".improved"
        if not os.path.isdir(rewrite_files): os.mkdir(rewrite_files)
        print("[singleNucleotideProbabilities] Rewriting files into {}\n".format(rewrite_files))
        rewritten_files = 0
    else:
        rewrite_files = None

    # look at all files
    orig_cnt = 0
    prob_cnt = 0
    rewr_cnt = 0
    norw_cnt = 0

    for file in files:
        identity, length, problem = validate_snp_file(file, full_reference_sequence,
                                                      print_sequences=False, print_summary=print_summary)
        ratio = 1.0 * identity / length
        all_identities.append(identity)
        all_lengths.append(length)
        all_identity_ratios.append(ratio)

        if rewrite_files is not None:
            if problem:
                print("{}: not rewriting problem file".format(os.path.basename(file)))
                prob_cnt += 1
                continue
            # make a new file
            new_file = os.path.join(rewrite_files, os.path.basename(file))
            if ratio < .5:
                rewrite_snp_file(file, new_file)
                rewritten_files += 1
                identity, length, _ = validate_snp_file(new_file, full_reference_sequence,
                                                     print_sequences=False, print_summary=False)
                ratio = 1.0 * identity / length
                is_better = ratio > .5
                print("{}:\trewritten file has {} length, {} identity, and ratio {}:" .format(
                    os.path.basename(file), length, identity, ratio))
                if is_better:
                    rewr_cnt += 1
                    print("{}:\trewritten file is acceptable".format( os.path.basename(file)))
                else:
                    os.remove(new_file)
                    norw_cnt += 1
                    print("{}:\trewritten file is still not acceptable.  removed from destination.".format( os.path.basename(file)))
            else:
                shutil.copyfile(file, new_file)
                orig_cnt += 1

    # printing results
    print("\n[singleNucleotideProbabilities] Summary of {} files:".format(len(files)))
    print("\tAVG Identity:       {}".format(np.mean(all_identities)))
    print("\tAVG Length:         {}".format(np.mean(all_lengths)))
    print("\tAVG Identity Ratio: {}".format(np.mean(all_identity_ratios)))
    print("\tIdentity Ratio:     {}".format(1.0 * sum(all_identities) / sum(all_lengths)))
    if rewrite_files is not None:
        total_cnt = orig_cnt + prob_cnt + rewr_cnt + norw_cnt
        print("[singleNucleotideProbabilities] Summary of rewriting:")
        print("\tIncl - Original:    {} ({}%)".format(orig_cnt, int(100 * orig_cnt / total_cnt)))
        print("\tIncl - Rewritten:   {} ({}%)".format(rewr_cnt, int(100 * rewr_cnt / total_cnt)))
        print("\tRem  - SNP Prob:    {} ({}%)".format(prob_cnt, int(100 * prob_cnt / total_cnt)))
        print("\tRem  - Bad Reverse: {} ({}%)".format(norw_cnt, int(100 * norw_cnt / total_cnt)))
        print("\n[singleNucleotideProbabilities] rewrote {} ({}%) files.  Rerunning validation."
              .format(rewritten_files, int(100 * rewritten_files / len(files))))
        validate_snp_directory(rewrite_files, reference_sequence_path, print_summary=False, rewrite_files=False)



def rewrite_snp_file(snp_file, output_file_location):
    with open(snp_file, 'r') as input, open(output_file_location, 'w') as output:
        for line in input:
            if line.startswith("#"):
                output.write(line)
            else:
                line = line.strip().split("\t")
                new_line = list()
                new_line.append(line[0]) #ch -> cg
                new_line.append(line[1]) #id -> id
                new_line.append(line[5]) #pT -> pA
                new_line.append(line[4]) #pG -> pC
                new_line.append(line[3]) #pC -> pG
                new_line.append(line[2]) #pA -> pT
                output.write("\t".join(new_line) + "\n")


def validate_snp_file(snp_file, full_reference_sequence, print_sequences=False, print_summary=False):
    identifier = os.path.basename(snp_file)
    consensus_sequence = list()
    header_positions = list()
    header_characters = list()
    first_pos = None
    last_pos = None
    problem = False
    duplicated_positions = 0
    unspecified_positions = 0
    with open(snp_file, 'r') as snp:
        for line in snp:
            if line.startswith("##"):
                continue
            elif line.startswith("#"):
                line = line.split("\t")
                i = 0
                for l in line:
                    if l.startswith("p"):
                        header_positions.append(i)
                        header_characters.append(l[1].upper())
                    i += 1
            else:
                line = line.split("\t")
                # positions
                pos = int(line[1])
                # set first_position (for reference matching)
                if first_pos is None: first_pos = pos
                # for cases where positions are duplicated or missing
                if last_pos is not None:
                    if last_pos >= pos:
                        duplicated_positions += 1
                        continue
                    while last_pos + 1 < pos:
                        consensus_sequence.append("-")
                        unspecified_positions += 1
                        last_pos += 1
                last_pos = pos
                #consensus
                max_prob = -1.0
                max_prob_idx = None
                idx = 0
                for pos in header_positions:
                    prob = float(line[pos].strip())
                    if prob > max_prob:
                        max_prob = prob
                        max_prob_idx = idx
                    idx += 1
                consensus_sequence.append(header_characters[max_prob_idx])

    # get sequences
    consensus_sequence = "".join(consensus_sequence)
    reference_sequence = full_reference_sequence[min(first_pos, last_pos):max(first_pos, last_pos)+1]

    # this is our quality metric
    length = 0
    identity = 0
    for c, r in zip(consensus_sequence, reference_sequence):
        length += 1
        if c == r: identity += 1

    # for separating results
    if print_sequences or print_summary: print("")

    # sanity check
    if duplicated_positions > 0:
        print("{}: Found {} duplicated positions!"
              .format(identifier, duplicated_positions))
    if unspecified_positions * 100 > length:
        print("{}: Found {} unspecified positions ({}% of total length)"
              .format(identifier, unspecified_positions, int(100.0 * unspecified_positions / length)))
        problem = True

    # printing full sequences
    if print_sequences:
        print("{}: Whole Sequences:".format(identifier))
        print("\treference:  {}".format(reference_sequence))
        print("\tconsensus:  {}".format(consensus_sequence))
        print("\tcomplement: {}".format(reverse_complement(consensus_sequence, complement=True, reverse=False)))

    if print_summary:
        print("%s:\tlength:%8d\t\tidentity:%8d\t\tratio:%f" % (identifier, length, identity, 1.0*identity/length))

    return identity, length, problem



def discover_single_nucleotide_probabilities(working_folder, kmer_length, reference_map, reference_sequence_string,
                                             list_of_fast5s, alignment_args, workers, step_size,
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

    # read fast5s and extract read ids
    fast5_to_read = build_fast5_to_read_id_dict(list_of_fast5s)
    print("[info] built map of fast5 identifiers to read ids with {} elements".format(len(fast5_to_read)))

    for s in xrange(step_size):
        print("\n[info] starting step %d / %d" % (s + 1, step_size))
        saved_step_dir = os.path.join(working_folder.path, "step_{}".format(s))
        scan_positions = range(s, reference_sequence_length, step_size)
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
            "step_size": step_size,
            "step_offset": s,
            "degenerate_type": alignment_args["degenerate"],
            "kmer_length": kmer_length
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
        if len(files) != step_size:
            print("[error] input fast5 '{}' yielded {} output files, expected {}".format(fast5_id, len(files), step_size))
            if len(files) == 0:
                continue

        # read all lines in all files
        output_lines = list()
        template_count = 0
        complement_count = 0
        for file in files:
            if "forward" in file.split("."):
                template_count += 1
            if "backward" in file.split("."):
                complement_count += 1
            with open(file, 'r') as input:
                for line in input:
                    line = line.split("\t")
                    line[0] = int(line[0])
                    output_lines.append(line)
        # sort based on position
        output_lines.sort(key=lambda x: x[0])

        # template/complement and sanity checks
        reverse = complement_count > template_count
        strand_identifier = "complement" if reverse else "template"
        if template_count == 0 and complement_count == 0:
            print("[warn] {}: could not infer template or complement".format(fast5_id))
            strand_identifier = "unknown"
        if template_count != 0 and complement_count != 0:
            print("[warn] {}: got {} template and {} complement calls after variant calling"
                  .format(fast5_id, template_count, complement_count))
            strand_identifier += " (by majority)"

        # write output
        output_filename = "{}.tsv".format(fast5_to_read[fast5_id])
        output_file = os.path.join(output_directory, output_filename)
        with open(output_file, 'w') as output:
            output.write("## fast5_input: {}.fast5\n".format(fast5_id))
            output.write("## read_id: {}\n".format(fast5_to_read[fast5_id]))
            output.write("## contig: {}\n".format(reference_map_contig_name))
            output.write("## strand: {}\n".format(strand_identifier))
            output.write("#CHROM\tPOS\tpA\tpC\tpG\tpT\n".format(reference_map_contig_name))
            for line in output_lines:
                if reverse:
                    line = [reference_map_contig_name, str(line[0]), line[2], line[3], line[4], line[5]]
                else:
                    line = [reference_map_contig_name, str(line[0]), line[5], line[4], line[3], line[2]]
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
    print("[singleNucleotideProbabilities] Command Line: {cmdLine}\n".format(cmdLine=command_line), file=sys.stderr)


    # first: see if we want to validate and return
    if args.validation_file is not None:
        if os.path.isfile(args.validation_file):
            full_reference_sequence = get_first_sequence(args.ref).upper()
            validate_snp_file(args.validation_file, full_reference_sequence, print_sequences=True, print_summary=True)
        elif os.path.isdir(args.validation_file):
            validate_snp_directory(args.validation_file, args.ref, print_summary=True)
        else:
            print("[error] got invalid validation location: {}".format(args.validation_file))
        return 0

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

    # assert integers
    args.step_size = int(args.step_size)
    args.kmer_size = int(args.kmer_size)

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
#   Step size: {stepSize}
    """.format(fileDir=args.files_dir, reference=args.ref, bwt=args.bwt, nbFiles=args.nb_files, banding=args.banded,
               inThmm=args.in_T_Hmm, inChmm=args.in_C_Hmm, model=args.stateMachineType, regions=args.target_regions,
               tHdp=args.templateHDP, cHdp=args.complementHDP, kmerSize=args.kmer_size, stepSize=args.step_size)

    print(start_message, file=sys.stdout)
    # cull the MinION files
    fast5s = cull_fast5_files(args.files_dir, args.nb_files)

    # get the (input) reference sequence
    if not os.path.isfile(args.ref):
        print("[singleNucleotideProbabilities] Did not find valid reference file", file=sys.stderr)
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
        print("[singleNucleotideProbabilities] using provided BWT %s" % args.bwt)
        bwa_ref_index = args.bwt
    else:
        print("[singleNucleotideProbabilities] indexing reference at %s" % args.ref, file=sys.stderr)
        bwa_ref_index = get_bwa_index(args.ref, temp_dir_path)
        print("[singleNucleotideProbabilities] indexing reference, done", file=sys.stderr)

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
        "threshold": args.threshold,
        "diagonal_expansion": args.diag_expansion,
        "constraint_trim": args.constraint_trim,
        "target_regions": None,
        "degenerate": degenerate_enum("variant"),
        "remove_temp_folder": False
    }
    #TODO you could save alignments by altering the above "destination" parameter

    # get the sites that have proposed edits
    print("\n\n[singleNucleotideProbabilities] scanning for proposals with %d fast5s" % len(fast5s))
    output_files = discover_single_nucleotide_probabilities(temp_folder, args.kmer_size, reference_map,
                                                            reference_sequence_string, fast5s, alignment_args,
                                                            args.nb_jobs, args.step_size, output_directory=args.out)
    print("\n[singleNucleotideProbabilities] got {} output files:".format(len(output_files)))
    i = 0
    for output_file in output_files:
        print("\t{}".format(output_file))
        i += 1
        if i > 10 and len(output_files) > 10:
            print("\t...")
            break

    #validation
    if len(output_files) != 0:
        validate_snp_directory(os.path.dirname(output_files[0]), args.ref, print_summary=True)

    print("\n\n[singleNucleotideProbabilities] fin\n")

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
