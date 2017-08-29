#!/usr/bin/env python
"""Library for calling variants
"""
from __future__ import print_function
import sys
import os
import glob
import pandas as pd
import numpy as np
from random import shuffle
from signalAlignLib import SignalAlignment
from alignmentAnalysisLib import CallMethylation
from multiprocessing import Process, Queue, current_process, Manager
from serviceCourse.parsers import read_fasta
from serviceCourse.sequenceTools import reverse_complement


def randomly_select_alignments(path_to_alignments, max_alignments_to_use):
    alignments = [x for x in glob.glob(path_to_alignments) if os.stat(x).st_size != 0]
    if len(alignments) == 0:
        print("[error] Didn't find any alignment files here {}".format(path_to_alignments))
        sys.exit(1)
    shuffle(alignments)

    if len(alignments) < max_alignments_to_use:
        return alignments
    else:
        return alignments[:max_alignments_to_use]


def get_forward_mask(list_of_alignments, suffix):
    mask = []
    for alignment in list_of_alignments:
        if alignment.endswith(".backward.tsv{}".format(suffix)):
            mask.append(False)
        else:
            mask.append(True)
    return mask


def get_alignments_labels_and_mask(path_to_alignments, max, suffix=""):
    alignments = randomly_select_alignments(path_to_alignments, max)
    mask = get_forward_mask(alignments, suffix)
    return alignments, mask


def get_reference_sequence(path_to_fasta):
    seqs = []

    for header, comment, sequence in read_fasta(path_to_fasta):
        seqs.append(sequence)

    assert len(seqs) > 0, "Didn't find any sequences in the reference file"

    if len(seqs) > 1:
        print("[NOTICE] Found more than one sequence in the reference file, using the first one")

    return seqs[0]


def make_degenerate_reference(input_sequence, positions, forward_sequence_path, backward_sequence_path,
                              block_size=1):
    """
    input_sequence: string, input nucleotide sequence
    out_path: string, path to directory to put new sequences with substituted degenerate characters
    block_size: not implemented, will be the size of the Ns to add (eg. NN = block_size 2)
    :return (subbed sequence, complement subbed sequence)
    """

    complement_sequence = reverse_complement(dna=input_sequence, reverse=False, complement=True)

    if positions is not None:
        t_seq = list(input_sequence)
        c_seq = list(complement_sequence)
        for position in positions:
            t_seq[position] = "X"
            c_seq[position] = "X"
        t_seq = ''.join(t_seq)
        c_seq = ''.join(c_seq)
    else:
        t_seq = input_sequence
        c_seq = complement_sequence

    with open(forward_sequence_path, 'w') as f:
        f.write("{seq}".format(seq=t_seq))
    with open(backward_sequence_path, 'w') as f:
        f.write("{seq}".format(seq=c_seq))

    return True


def load_variant_call_data(file_path):
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


def call_sites_with_marginal_probs(data, reference_sequence_string, min_depth=0, get_sites=False):
    d = load_variant_call_data(data)

    reference_sequence_list = list(reference_sequence_string)

    candidate_sites = []
    add_to_candidates = candidate_sites.append

    # print("#POS\tA\tC\tG\tT\tDIFF")
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

        normed_marginal_probs = marginal_prob.map(lambda y: y / sum(marginal_prob))
        called_base = normed_marginal_probs.argmax()[1]
        difference = normed_marginal_probs.max() - normed_marginal_probs["p" + reference_sequence_list[site]]
        add_to_candidates((site, difference))

        # print("{}\t{}\t{}\t{}\t{}\t{}".format(site, normed_marginal_probs['pA'], normed_marginal_probs['pC'],
        #                                   normed_marginal_probs['pG'], normed_marginal_probs['pT'], difference))
        # if called_base != reference_sequence_list[site]:
        #     if get_sites is False:
        #         print("Changing {orig} to {new} at {site} depth {depth}"
        #               "".format(orig=reference_sequence_list[site], new=called_base, site=site, depth=len(x['read'])))
        #         reference_sequence_list[site] = called_base
        #     else:
        #         print("Proposing edit at {site} from {orig} to {new}, \n{probs}"
        #               "".format(orig=reference_sequence_list[site], new=called_base, site=site,
        #                         probs=normed_marginal_probs))
        #         difference = normed_marginal_probs.max() - normed_marginal_probs["p" + reference_sequence_list[site]]
        #         print(difference)
        #         add_to_candidates((site, difference))

    if get_sites is True:
        return candidate_sites
    else:
        return ''.join(reference_sequence_list)


def aligner(work_queue, done_queue):
    try:
        total_alignments = 0
        failed_alignments = 0
        for f in iter(work_queue.get, 'STOP'):
            try:
                alignment = SignalAlignment(**f)
                alignment.run()
            except Exception, e:
                name = current_process().name
                message = e.message
                if message is None or len(message) == 0:
                    message = "exception:" + str(e)
                error = "aligner '{}' failed with: {}".format(name, message)
                print("[aligner] " + error)
                done_queue.put(error)
                failed_alignments += 1
            total_alignments += 1
        print("[aligner] '%s' completed %d alignments with %d failures"
              % (current_process().name, total_alignments, failed_alignments))
    except Exception, e:
        name = current_process().name
        message = e.message
        if message is None or len(message) == 0:
            message = "exception:" + str(e)
        error = "CRITICAL: aligner '{}' failed with: {}".format(name, message)
        print("[aligner] " + error)
        done_queue.put(error)


def variant_caller(work_queue, done_queue):
    try:
        i = 0
        for f in iter(work_queue.get, 'STOP'):
            print("[variant_caller] '{}' invoking CallMethylation with {}".format(current_process().name, f))
            c = CallMethylation(**f)
            c.write()
            i += 1
        print("[variant_caller] '%s' completed %d variant calls" % (current_process().name, i))
    except Exception, e:
        name = current_process().name
        message = e.message
        if message is None or len(message) == 0:
            message = "exception:" + str(e)
        error = "variant_caller '{}' failed with: {}".format(name, message)
        print("[variant_caller] " + error)
        done_queue.put(error)


def run_service(service, service_iterable, service_arguments, workers, iterable_argument):
    args = service_arguments.keys()
    args.append(iterable_argument)
    print("[debug] running service {} with {} workers, {} iterables, and arguments {}"
          .format(service, workers, len(service_iterable), args))
    # setup workers for multiprocessing
    work_queue = Manager().Queue()
    done_queue = Manager().Queue()

    jobs = []
    for x in service_iterable:
        args = dict({iterable_argument: x},
                    **service_arguments)
        work_queue.put(args)

    for w in xrange(workers):
        p = Process(target=service, args=(work_queue, done_queue))
        p.start()
        jobs.append(p)
        work_queue.put('STOP')

    for p in jobs:
        p.join()

    done_queue.put('STOP')


#todo: deprecated (this should be removed)
def make_reference_files_and_alignment_args(working_folder, reference_sequence_string, reference_map,
                                            n_positions=None):
    # make paths for working txt files that contain this STEPs Ns
    forward_reference = working_folder.add_file_path("forward_reference.txt")
    backward_reference = working_folder.add_file_path("backward_reference.txt")

    # make N-ed reference sequence for this iteration, writes the strings to files

    check = make_degenerate_reference(reference_sequence_string, n_positions,
                                      forward_reference, backward_reference)
    assert check, "Problem making degenerate reference"

    # perform alignment for this step
    #instead of updating the alignment_args with these, we now update values in the map for this single contig
    reference_map["forward"] = forward_reference #key was forward_reference
    reference_map["backward"] = backward_reference #key was backward_reference
    return True


def scan_for_proposals(working_folder, step, reference_map, reference_sequence_string, list_of_fast5s, alignment_args,
                       workers, use_saved_alignments=True, save_alignments=True):
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
    proposals = []

    for s in xrange(step):
        print("\n[info] starting step %d" % s)
        saved_step_dir = os.path.join(working_folder, "step_{}".format(s))
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
            print("[info] Found %d alignment files here %s" % (alignment_count, working_folder.path))

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

        # get proposal sites
        prob_glob = os.path.join(working_folder.path, "{}*".format(marginal_probability_prefix))
        marginal_prob_files = [x for x in
                               glob.glob(prob_glob)
                               if os.stat(x).st_size != 0]
        print("[info] found {} marginal probability files matching {}".format(len(marginal_prob_files), prob_glob))
        for marginal_prob_file in marginal_prob_files:
            print("[info] calling sites with marginal probabilities from %s" % marginal_probability_prefix)
            proposals += call_sites_with_marginal_probs(marginal_prob_file, reference_sequence_string,
                                                        min_depth=0, get_sites=True)
        # remove or old alignments
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

    # proposals is a list of lists containing (position, delta_prob) where position in the position in the
    # reference sequence that is being proposed to be edited, and delta_prob is the difference in probability
    # of the reference base to the proposed base
    print("[info] returning %d proposals" % len(proposals))
    return proposals


def update_reference_with_marginal_probs(working_folder, proposals, reference_sequence_string, list_of_fast5s,
                                         alignment_args, workers):
    check = make_reference_files_and_alignment_args(working_folder, reference_sequence_string, alignment_args,
                                                    n_positions=proposals)
    assert check, "[update_reference_with_marginal_probs]: problem making reference files and args dict"
    run_service(aligner, list_of_fast5s, alignment_args, workers, "in_fast5")

    alignments = [x for x in glob.glob(working_folder.path + "*.tsv") if os.stat(x).st_size != 0]

    marginal_probability_file = working_folder.add_file_path("proposals.calls")

    proposal_args = {
        "sequence": None,
        "out_file": marginal_probability_file,
        "positions": {"forward": proposals, "backward": proposals},
        "degenerate_type": alignment_args["degenerate"]
    }
    #for alignment in alignments:
    #    a = dict({"alignment_file": alignment}, **proposal_args)
    #    c = CallMethylation(**a)
    #    c.write()
    run_service(variant_caller, alignments, proposal_args, workers, "alignment_file")

    # get proposal sites
    updated_reference_sequence = call_sites_with_marginal_probs(marginal_probability_file, reference_sequence_string,
                                                                min_depth=0, get_sites=True)
    # clean up
    working_folder.remove_file(marginal_probability_file)

    # remove old alignments
    for f in glob.glob(working_folder.path + "*.tsv"):
        os.remove(f)

    return updated_reference_sequence








