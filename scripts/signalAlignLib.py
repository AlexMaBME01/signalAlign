"""Small library for working with MinION data
"""
from __future__ import print_function, division
import os
import h5py
import sys
sys.path.append("../")
import subprocess
import re
import numpy as np
from itertools import islice, izip
from random import shuffle
from serviceCourse.sequenceTools import reverse_complement
from serviceCourse.parsers import read_fasta
from serviceCourse.file_handlers import FolderHandler

# Globals
NORM_DIST_PARAMS = 2
NB_MODEL_PARAMS = 5

def kmer_iterator(dna, k):
    for i in xrange(len(dna)):
        kmer = dna[i:(i+k)]
        if len(kmer) == k:
            yield kmer


def write_fasta(id, sequence, destination):
    print(">", id, sep="", end="\n", file=destination)
    print(sequence, end="\n", file=destination)
    destination.close()
    return


def cull_fast5_files(path_to_files, maximum_files):
    # list of alignment files
    fast5s = [x for x in os.listdir(path_to_files) if x.endswith(".fast5")]
    fast5s = [path_to_files + x for x in fast5s]

    if len(fast5s) == 0 or fast5s is None:
        print("[cull_fast5_files] : error culling .fast5 files")
        sys.exit(1)

    # take only some
    if maximum_files < len(fast5s):
        shuffle(fast5s)
        fast5s = fast5s[:maximum_files]

    return fast5s


def get_bwa_index(reference, dest):
    bwa = Bwa(reference)
    bwa.build_index(dest)
    bwa_ref_index = dest + "temp_bwaIndex"
    return bwa_ref_index


def get_npRead_2dseq_and_models(fast5, npRead_path, twod_read_path):
    """process a MinION .fast5 file into a npRead file for use with signalAlign also extracts
    the 2D read into fasta format
    """
    # setup
    out_file = open(npRead_path, 'w')
    temp_fasta = open(twod_read_path, "w")

    # load MinION read
    npRead = NanoporeRead(fast5)

    # only working with 2D reads right now
    if npRead.has2D_alignment_table is False:
        npRead.close()
        return False, None, False

    proceed = npRead.write_npRead(out_file=out_file)

    if proceed:
        # make the 2d read
        write_fasta(id=fast5, sequence=npRead.alignment_table_sequence, destination=temp_fasta)

        if npRead.complement_model_id == "complement_median68pA_pop1.model":
            pop1_complement = True
        else:
            pop1_complement = False
        version = npRead.version
        npRead.close()
        return True, version, pop1_complement
    else:
        npRead.close()
        print("problem making npRead for {fast5}".format(fast5=fast5), file=sys.stderr)
        return False, None, False


def parse_substitution_file(substitution_file):
    fH = open(substitution_file, 'r')
    line = fH.readline().split()
    forward_sub = line[0]
    forward_pos = map(np.int64, line[1:])
    line = fH.readline().split()
    backward_sub = line[0]
    backward_pos = map(np.int64, line[1:])
    return (forward_sub, forward_pos), (backward_sub, backward_pos)


def make_temp_sequence(fasta, sequence_outfile, rc_sequence_outfile):
    """extract the sequence from a fasta and put into a simple file that is used by signalAlign
    """
    assert (not os.path.isfile(sequence_outfile)), "ERROR: forward file already exists"
    assert (not os.path.isfile(rc_sequence_outfile)), "ERROR: backward file already exists"

    for header, comment, sequence in read_fasta(fasta):
        print(sequence, end='\n', file=open(sequence_outfile, 'w'))
        complement_sequence = reverse_complement(sequence, reverse=False, complement=True)
        print(complement_sequence, end='\n', file=open(rc_sequence_outfile, 'w'))
        break


def add_ambiguity_chars_to_reference(input_fasta, substitution_file, sequence_outfile, rc_sequence_outfile,
                                     degenerate_type, sub_out="C", ambig_char="X"):
    assert os.path.isfile(input_fasta), "ERROR: Didn't find reference FASTA {}".format(input_fasta)
    assert os.path.isfile(substitution_file), "ERROR: Didn't find substitution file {}".format(substitution_file)
    assert (not os.path.isfile(sequence_outfile)), "ERROR: forward file already exists"
    assert (not os.path.isfile(rc_sequence_outfile)), "ERROR: forward file already exists"

    # get the first sequence from the FASTA
    seq = ""
    for header, comment, sequence in read_fasta(input_fasta):
        seq += sequence
        break

    # we want the complement, not the reverse complement, we actually flip it around later
    r_seq = reverse_complement(dna=seq, reverse=False, complement=True)

    # turn the sequence into a list so we can change the nucleotides
    seq = list(seq)
    r_seq = list(r_seq)

    # parse the substitution file
    f, b = parse_substitution_file(substitution_file=substitution_file)
    forward_pos = f[1]
    backward_pos = b[1]

    for position in forward_pos:
        if degenerate_type in ["twoWay", "threeWay"]:
            assert seq[position] in list(sub_out), "ERROR: trying to sub {seq_pos} not allowed"\
                .format(seq_pos=seq[position])
        seq[position] = ambig_char
    for position in backward_pos:
        if degenerate_type in ["twoWay", "threeWay"]:
            assert r_seq[position] in list(sub_out), "ERROR: trying to sub {seq_pos} not allowed"\
                .format(seq_pos=seq[position])
        r_seq[position] = ambig_char

    # make them back into strings
    seq = ''.join(seq)
    r_seq = ''.join(r_seq)

    # write to files
    print(seq, end='\n', file=open(sequence_outfile, "w"))
    print(r_seq, end='\n', file=open(rc_sequence_outfile, "w"))
    return


def parse_cigar(cigar_string, ref_start):
    # use a regular expression to parse the string into operations and lengths
    cigar_tuples = re.findall(r'([0-9]+)([MIDNSHPX=])', cigar_string)

    clipping = {"S", "H"}
    alignment_operations = {"M", "I", "D"}

    # make some containers
    query_start = 0
    past_start = False
    query_end = 0
    reference_start = ref_start - 1  # fence posts adjustment
    reference_end = 0

    exonerated_cigar = " ".join(["%s %i" % (operation, int(length)) for length, operation in
                                 cigar_tuples if operation in alignment_operations])

    # this is how you calculate the reference map region
    for length, op in cigar_tuples:
        if op in clipping and past_start is False:
            query_start += int(length)
        if op == "M" or op == "D":
            reference_end += int(length)
            if past_start is False:
                past_start = True
        if op == "M" or op == "I":
            query_end += int(length)
            if past_start is False:
                past_start = True

    query_end = query_end + query_start
    reference_end = reference_end + reference_start

    return query_start, query_end, reference_start, reference_end, exonerated_cigar


def exonerated_bwa(bwa_index, query, target_regions=None):
    # align with bwa
    command = "bwa mem -x ont2d {index} {query}".format(index=bwa_index, query=query)

    # this is a small SAM file that comes from bwa
    aln = subprocess.check_output(command.split())
    aln = aln.split("\t")  # split

    query_start, query_end, reference_start, reference_end, cigar_string = parse_cigar(aln[11], int(aln[9]))

    strand = ""
    if int(aln[7]) == 16:
        # todo redo this swap
        strand = "-"
        temp = reference_start
        reference_start = reference_end
        reference_end = temp
    if int(aln[7]) == 0:
        strand = "+"
    elif int(aln[7]) != 0 and int(aln[7]) != 16:
        print("unknown alignment flag, exiting", file=sys.stderr)
        return False, False

    completeCigarString = "cigar: %s %i %i + %s %i %i %s 1 %s" % (
    aln[6].split()[-1], query_start, query_end, aln[8], reference_start, reference_end, strand, cigar_string)

    if target_regions is not None:
        keep = target_regions.check_aligned_region(reference_start, reference_end)
        if keep is False:
            return False, False
        else:
            pass

    return completeCigarString, strand


def default_template_model_from_version(version):
    supported_versions = ["1.15.0", "1.19.0", "1.20.0", "1.22.2", "1.22.4"]
    assert version in supported_versions
    version_index = supported_versions.index(version)

    if version_index <= 2:
        r7_3_default_template_model = "../models/testModel_template.model"
        assert os.path.exists(r7_3_default_template_model), "Didn't find default template R7.3 model"
        return r7_3_default_template_model
    else:
        r9_default_template_model = "../models/testModelR9_template.model"
        assert os.path.exists(r9_default_template_model), "Didn't find default template R9 model"
        return r9_default_template_model


def default_complement_model_from_version(version, pop1_complement=False):
    supported_versions = ["1.15.0", "1.19.0", "1.20.0", "1.22.2", "1.22.4"]
    assert version in supported_versions
    version_index = supported_versions.index(version)

    if version_index <= 2:
        r7_3_default_complement_model = "../models/testModel_complement.model" if not pop1_complement \
            else "../models/testModelR9_complement_pop2.model"
        assert os.path.exists(r7_3_default_complement_model), "Didn't find default complement R7.3 model"
        return r7_3_default_complement_model
    else:
        r9_default_complement_model = "../models/testModelR9_complement.model"
        assert os.path.exists(r9_default_complement_model), "Didn't find default complement R9 model"
        return r9_default_complement_model



def degenerate_enum(degenerate_request_string):
    degenerate_type = {
        "twoWay": 0,
        "threeWay": 1,
        "variant": 3,
    }

    assert (degenerate_request_string in degenerate_type.keys()), "Requested degenerate nucleotide set not recognized."
    return degenerate_type[degenerate_request_string]


class TargetRegions(object):
    def __init__(self, tsv, already_sorted=False):
        assert(os.stat(tsv).st_size != 0), "Empty regions file"

        self.region_array = np.loadtxt(tsv, usecols=(0, 1), dtype=np.int32)

        if len(self.region_array.shape) == 1:
            a = np.empty([1, 2], dtype=np.int32)
            a[0] = self.region_array
            self.region_array = a

        if not already_sorted:
            self.region_array = np.sort(self.region_array, axis=1)

    def check_aligned_region(self, left, right):
        if right < left:
            left, right = right, left
        for region in self.region_array:
            if (region[0] >= left) and (region[1] <= right):
                return True
            else:
                continue
        return False


class Bwa(object):
    """run BWA, mostly used to make index files
    """
    def __init__(self, target):
        self.target = target
        self.db_handle = ''

    def build_index(self, destination):
        # make a place to put the database
        path_to_bwa_index = destination

        # build database
        self.db_handle = path_to_bwa_index + '/temp_bwaIndex'
        os.system("bwa index -p {0} {1}".format(self.db_handle, self.target))


class NanoporeRead(object):
    def __init__(self, fast_five_file):
        # load the fast5
        self.filename = fast_five_file
        self.is_open = self.open()
        self.template_event_map = []
        self.complement_event_map = []
        self.stay_prob = 0
        self.skip_prob_bins = []
        self.template_model_name = ""
        self.complement_model_name = ""
        self.initialize_twoD()

    def open(self):
        try:
            self.fastFive = h5py.File(self.filename, 'r')
            return True
        except Exception, e:
            self.close()
            print("Error opening file {filename}".format(filename=self.filename), file=sys.stderr)
            return False

    def initialize_twoD(self, get_sequence=False):
        # init
        self.has2D = False
        self.has2D_alignment_table = False

        twoD_alignment_table_address = "/Analyses/Basecall_2D_000/BaseCalled_2D/Alignment"
        if twoD_alignment_table_address in self.fastFive:
            self.twoD_alignment_table = self.fastFive[twoD_alignment_table_address]
            if len(self.twoD_alignment_table) > 0:
                self.has2D_alignment_table = True
            self.kmer_length = len(self.twoD_alignment_table[0][2])

        if get_sequence is True:
            twoD_read_sequence_address = "/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq"
            if twoD_read_sequence_address in self.fastFive:
                self.has2D = True
                self.twoD_read_sequence = self.fastFive[twoD_read_sequence_address][()].split()[2]
                self.twoD_id = self.fastFive[twoD_read_sequence_address][()].split()[0:2][0][1:]

        supported_versions = ["1.15.0", "1.19.0", "1.20.0", "1.22.2", "1.22.4"]
        self.version = self.fastFive["/Analyses/Basecall_2D_000"].attrs["dragonet version"]

        if self.version not in supported_versions:
            print("Unsupported Version (1.15.0, 1.19.0, 1.20.0, 1.22.2, 1.22.4 supported)", file=sys.stdout)
            return False

        # initialize version-specific paths
        if self.version == "1.15.0":
            self.template_event_table_address = '/Analyses/Basecall_2D_000/BaseCalled_template/Events'
            self.template_model_address = "/Analyses/Basecall_2D_000/BaseCalled_template/Model"
            self.template_model_id = self.get_model_id("/Analyses/Basecall_2D_000/Summary/basecall_1d_template")
            self.template_read = self.fastFive["/Analyses/Basecall_2D_000/BaseCalled_template/Fastq"][()].split()[2]

            self.complement_event_table_address = '/Analyses/Basecall_2D_000/BaseCalled_complement/Events'
            self.complement_model_address = "/Analyses/Basecall_2D_000/BaseCalled_complement/Model"
            self.complement_model_id = self.get_model_id("/Analyses/Basecall_2D_000/Summary/basecall_1d_complement")
            self.complement_read = self.fastFive["/Analyses/Basecall_2D_000/BaseCalled_complement/Fastq"][()].split()[2]
            return True

        elif self.version == "1.19.0" or self.version == "1.20.0":
            self.template_event_table_address = '/Analyses/Basecall_1D_000/BaseCalled_template/Events'
            self.template_model_address = "/Analyses/Basecall_1D_000/BaseCalled_template/Model"
            self.template_model_id = self.get_model_id("/Analyses/Basecall_1D_000/Summary/basecall_1d_template")
            self.template_read = self.fastFive["/Analyses/Basecall_1D_000/BaseCalled_template/Fastq"][()].split()[2]

            self.complement_event_table_address = '/Analyses/Basecall_1D_000/BaseCalled_complement/Events'
            self.complement_model_address = "/Analyses/Basecall_1D_000/BaseCalled_complement/Model"
            self.complement_model_id = self.get_model_id("/Analyses/Basecall_1D_000/Summary/basecall_1d_complement")
            self.complement_read = self.fastFive["/Analyses/Basecall_1D_000/BaseCalled_complement/Fastq"][()].split()[2]
            return True

        elif self.version == "1.22.2" or self.version == "1.22.4":
            self.template_event_table_address = '/Analyses/Basecall_1D_000/BaseCalled_template/Events'
            self.template_model_address = ""
            self.template_model_id = None
            self.template_read = self.fastFive["/Analyses/Basecall_1D_000/BaseCalled_template/Fastq"][()].split()[2]

            self.complement_event_table_address = '/Analyses/Basecall_1D_000/BaseCalled_complement/Events'
            self.complement_model_address = ""
            self.complement_model_id = None
            self.complement_read = self.fastFive["/Analyses/Basecall_1D_000/BaseCalled_complement/Fastq"][()].split()[2]
            return True
        else:
            print("Unsupported Version (1.15.0, 1.19.0, 1.20.0, 1.22.2, 1.22.4 supported)", file=sys.stdout)
            return False

    def assemble_2d_sequence_from_table(self):
        """The 2D read sequence contains kmers that may not map to a template or complement event, which can make
        mapping difficult downstream. This function makes a sequence from the 2D alignment table, which is usually
        pretty similar to the 2D read, except it is guaranteed to have an event map to every position.

        returns: sequence made from alignment table
        """
        def find_kmer_overlap(k_i, k_j):
            """ finds the overlap between two non-identical kmers.
            k_i: one kmer
            k_j: another kmer
            returns: The number of positions not matching
            """
            for i in xrange(1, len(k_i)):
                sk_i = k_i[i:]
                sk_j = k_j[:-i]
                if sk_i == sk_j:
                    return i
            return len(k_i)

        self.alignment_table_sequence = ''
        self.alignment_table_sequence = self.twoD_alignment_table[0][2]
        p_kmer = self.twoD_alignment_table[0][2]

        # iterate through the 6-mers in the alignment table
        for t, c, kmer in self.twoD_alignment_table:
            # if we're at a new 6-mer
            if kmer != p_kmer:
                # find overlap, could move up to len(6-mer) - 1 bases
                i = find_kmer_overlap(p_kmer, kmer)

                # append the suffix of the new 6-mer to the sequence
                self.alignment_table_sequence += kmer[-i:]

                # update
                p_kmer = kmer
            else:
                continue
        return

    def init_1d_event_maps(self):
        """Maps the events from the template and complement strands to their base called kmers the map
        generated by this function is called the "strand_event_map" because it only works for mapping the
        strand read (1D read) to to it's events. Uses the same fields as 'get_twoD_event_map' below.
        """
        def make_map(events):
            event_map = [0]
            previous_prob = 0
            for i, line in islice(enumerate(events), 1, None):
                move = line['move']
                this_prob = line['p_model_state']
                if move == 1:
                    event_map.append(i)
                if move > 1:
                    for skip in xrange(move - 1):
                        event_map.append(i - 1)
                    event_map.append(i)
                if move == 0:
                    if this_prob > previous_prob:
                        event_map[-1] = i
                previous_prob = this_prob
            final_event_index = [event_map[-1]]
            padding = final_event_index * (self.kmer_length - 1)
            event_map = event_map + padding
            return event_map
        self.template_strand_event_map = make_map(self.template_events)
        self.complement_strand_event_map = make_map(self.complement_events)
        assert len(self.template_strand_event_map) == len(self.template_read)
        assert len(self.complement_strand_event_map) == len(self.complement_read)
        return True

    def get_twoD_event_map(self):
        """Maps the kmers in the alignment table sequence read to events in the template and complement strand reads
        """
        # initialize
        alignment_row = 0
        prev_alignment_kmer = ''
        nb_template_gaps = 0
        previous_complement_event = None
        previous_template_event = None

        #twoD_init = self.initialize_twoD()
        #if twoD_init is False:
        #    return False

        if not self.has2D_alignment_table:
            print("{file} doesn't have 2D alignment table".format(file=self.filename))
            return False

        self.assemble_2d_sequence_from_table()

        # go thought the kmers in the read sequence and match up the events
        for i, seq_kmer in enumerate(kmer_iterator(self.alignment_table_sequence, self.kmer_length)):
            # assign the current row's kmer
            current_alignment_kmer = self.twoD_alignment_table[alignment_row][2]

            # in the situation where there is a repeat kmer in the alignment then
            # we want to pick the best event to kmer alignment, TODO implement this
            # right now we just use the first alignment
            while current_alignment_kmer == prev_alignment_kmer:
                alignment_row += 1
                current_alignment_kmer = self.twoD_alignment_table[alignment_row][2]

            # a match
            if seq_kmer == current_alignment_kmer:
                template_event = self.twoD_alignment_table[alignment_row][0]
                complement_event = self.twoD_alignment_table[alignment_row][1]

                # handle template event
                # if there is a gap, count it and don't add anything to the map
                if template_event == -1:
                    nb_template_gaps += 1

                # if there is an aligned event
                if template_event != -1:
                    # if it is an aligned event and there are no gaps, add it to the map
                    if nb_template_gaps == 0:
                        self.template_event_map.append(template_event)
                        # update
                        previous_template_event = template_event
                    # if there were gaps in the alignment we have to add 'best guess'
                    # event alignments to the map which is the current aligned event
                    if nb_template_gaps > 0:
                        self.template_event_map += [template_event] * (nb_template_gaps + 1)
                        # reset template gaps
                        nb_template_gaps = 0
                        # update
                        previous_template_event = template_event

                # handle complement event
                # if there is a gap, add the last aligned complement event to the map
                if complement_event == -1:
                    self.complement_event_map.append(previous_complement_event)

                # if there is an aligned complement event add it to the map
                if complement_event != -1:
                    self.complement_event_map.append(complement_event)
                    # update the most recent aligned complement event
                    previous_complement_event = complement_event

                # update previous alignment kmer and increment alignment row
                prev_alignment_kmer = current_alignment_kmer
                alignment_row += 1
                continue

            # not a match, meaning that this kmer in the read sequence is not
            # in the event alignment but we need to assign an event to it so
            # we use the heuristic that we use the alignment of the most
            # recent aligned events to this base
            if seq_kmer != current_alignment_kmer:
                self.template_event_map.append(previous_template_event)
                self.complement_event_map.append(previous_complement_event)
                continue

        # fill in the final events for the partial last kmer
        for _ in xrange(self.kmer_length - 1):
            self.template_event_map += [previous_template_event] * (nb_template_gaps + 1)
            self.complement_event_map.append(previous_complement_event)
            nb_template_gaps = 0

        # check that we have mapped all of the bases in the 2D read
        assert(len(self.template_event_map) == len(self.alignment_table_sequence))
        assert(len(self.complement_event_map) == len(self.alignment_table_sequence))
        return True

    def adjust_events_for_drift(self, events, drift):
        """Adjust event means by drift
        """
        if events is None or drift is None:
            return False

        # transform events by time
        # events have format [[mean], [start_time], [std_dev], [length]]
        # get the start time of the first event
        start_time = events[0][1]
        for event in events:
            # time since first event
            delta_time = event[1] - start_time
            # drift adjust
            # TODO change adjustment here
            event[0] -= (delta_time * drift)
        return True

    def get_template_events(self):
        if self.template_event_table_address in self.fastFive:
            self.template_events = self.fastFive[self.template_event_table_address]
            return True

        if self.template_event_table_address not in self.fastFive:
            return False

    def get_complement_events(self):
        if self.complement_event_table_address in self.fastFive:
            self.complement_events = self.fastFive[self.complement_event_table_address]
            return True

        if self.complement_event_table_address not in self.fastFive:
            return False

    def get_template_model_adjustments(self):
        if self.template_model_address in self.fastFive:
            self.has_template_model = True
            self.template_scale = self.fastFive[self.template_model_address].attrs["scale"]
            self.template_shift = self.fastFive[self.template_model_address].attrs["shift"]
            self.template_drift = self.fastFive[self.template_model_address].attrs["drift"]
            self.template_var = self.fastFive[self.template_model_address].attrs["var"]
            self.template_scale_sd = self.fastFive[self.template_model_address].attrs["scale_sd"]
            self.template_var_sd = self.fastFive[self.template_model_address].attrs["var_sd"]

        if self.template_model_address not in self.fastFive:
            self.has_template_model = False
            self.template_scale = 1
            self.template_shift = 1
            self.template_drift = 0
            self.template_var = 1
            self.template_scale_sd = 1
            self.template_var_sd = 1
        return

    def get_complement_model_adjustments(self):
        if self.complement_model_address in self.fastFive:
            self.has_complement_model = True
            self.complement_scale = self.fastFive[self.complement_model_address].attrs["scale"]
            self.complement_shift = self.fastFive[self.complement_model_address].attrs["shift"]
            self.complement_drift = self.fastFive[self.complement_model_address].attrs["drift"]
            self.complement_var = self.fastFive[self.complement_model_address].attrs["var"]
            self.complement_scale_sd = self.fastFive[self.complement_model_address].attrs["scale_sd"]
            self.complement_var_sd = self.fastFive[self.complement_model_address].attrs["var_sd"]

        if self.complement_model_address not in self.fastFive:
            self.has_complement_model = False
            self.complement_scale = 1
            self.complement_shift = 1
            self.complement_drift = 0
            self.complement_var = 1
            self.complement_scale_sd = 1
            self.complement_var_sd = 1
        return

    @staticmethod
    def calculate_lambda(noise_mean, noise_stdev):
        return (np.power(noise_mean, 3)) / (np.power(noise_stdev, 2))

    def export_model(self, skip_bins, model_address, destination):
        """Exports the model to a file. Format:
        line 1: [correlation coefficient] [level_mean] [level_sd] [noise_mean]
                    [noise_sd] [noise_lambda ] (.../kmer) \n
        line 2: skip bins \n
        line 3: [correlation coefficient] [level_mean] [level_sd, scaled]
                    [noise_mean] [noise_sd] [noise_lambda ] (.../kmer) \n
        """

        assert self.is_open, "ERROR: Fast5 file is not open"

        lambdas = []

        if model_address in self.fastFive:
            model = self.fastFive[model_address]
            # line 1
            print("0", end=' ', file=destination)  # placeholder for correlation parameter
            for kmer, level_mean, level_sd, noise_mean, noise_sd, weight in model:
                lam = self.calculate_lambda(noise_mean, noise_sd)
                lambdas.append(lam)
                print(level_mean, level_sd, noise_mean, noise_sd, lam, end=' ', file=destination)
            print("", end="\n", file=destination)
            # line 2
            for p in skip_bins:
                print(p, end=' ', file=destination)
            print("", end="\n", file=destination)
            # line 3
            print("0", end=' ', file=destination) # placeholder for correlation parameter
            i = 0
            for kmer, level_mean, level_sd, noise_mean, noise_sd, weight in model:
                lam = lambdas[i]
                print(level_mean, (level_sd * 1.75), noise_mean, noise_sd, lam, end=' ', file=destination)
                i += 1
            print("", end="\n", file=destination)
            return True
        else:
            return False

    def export_template_model(self, destination):
        # for conditional HMM (as per JTS)
        t_skip_prob_bins = [0.487, 0.412, 0.311, 0.229, 0.174, 0.134, 0.115, 0.103, 0.096, 0.092,
                            0.088, 0.087, 0.084, 0.085, 0.083, 0.082, 0.085, 0.083, 0.084, 0.082,
                            0.080, 0.085, 0.088, 0.086, 0.087, 0.089, 0.085, 0.090, 0.087, 0.096]

        got_model = self.export_model(t_skip_prob_bins, self.template_model_address, destination)

        return got_model

    def export_complement_model(self, destination):
        c_skip_prob_bins = [0.531, 0.478, 0.405, 0.327, 0.257, 0.207, 0.172, 0.154, 0.138, 0.132,
                            0.127, 0.123, 0.117, 0.115, 0.113, 0.113, 0.115, 0.109, 0.109, 0.107,
                            0.104, 0.105, 0.108, 0.106, 0.111, 0.114, 0.118, 0.119, 0.110, 0.119]

        got_model = self.export_model(c_skip_prob_bins, self.complement_model_address, destination)

        return got_model

    def get_model_id(self, address):
        if address in self.fastFive:
            model_name = self.fastFive[address].attrs["model_file"]
            model_name = model_name.split('/')[-1]
            return model_name
        else:
            return None

    def write_npRead(self, out_file):
        if self.is_open is False:
            print("problem opeining file {filename}".format(filename=self.filename), file=sys.stderr)
            self.close()
            return False

        twoD_map_check = self.get_twoD_event_map()
        template_events_check = self.get_template_events()
        complement_events_check = self.get_complement_events()
        oneD_event_map_check = self.init_1d_event_maps()

        proceed = False not in [twoD_map_check, template_events_check, complement_events_check, oneD_event_map_check]

        if proceed:
            # get model params
            self.get_template_model_adjustments()
            self.get_complement_model_adjustments()

            # transform events
            # drift adjustment happens within signalMachine now
            #if self.version in ["1.15.0", "1.19.0"]:
            #    t_transformed = self.adjust_events_for_drift(self.template_events, self.template_drift)
            #    c_transformed = self.adjust_events_for_drift(self.complement_events, self.complement_drift)
                # check if that worked
            #    if t_transformed is False or c_transformed is False:
            #        return False

            # Make the npRead

            # line 1 parameters
            print(len(self.alignment_table_sequence), end=' ', file=out_file)  # 0alignment read length
            print(len(self.template_events), end=' ', file=out_file)           # 1nb of template events
            print(len(self.complement_events), end=' ', file=out_file)         # 2nb of complement events
            print(len(self.template_read), end=' ', file=out_file)             # 3length of template read
            print(len(self.complement_read), end=' ', file=out_file)           # 4length of template read
            print(self.template_scale, end=' ', file=out_file)                 # 5template scale
            print(self.template_shift, end=' ', file=out_file)                 # 67template shift
            print(self.template_var, end=' ', file=out_file)                   # 7template var
            print(self.template_scale_sd, end=' ', file=out_file)              # 8template scale_sd
            print(self.template_var_sd, end=' ', file=out_file)                # 9template var_sd
            print(self.template_drift, end=' ', file=out_file)                 # 0template_drift
            print(self.complement_scale, end=' ', file=out_file)               # 1complement scale
            print(self.complement_shift, end=' ', file=out_file)               # 2complement shift
            print(self.complement_var, end=' ', file=out_file)                 # 3complement var
            print(self.complement_scale_sd, end=' ', file=out_file)            # 4complement scale_sd
            print(self.complement_var_sd, end=' ', file=out_file)              # 5complement var_sd
            print(self.complement_drift, end='\n', file=out_file)              # 6complement_drift

            # line 2 alignment table sequence
            print(self.alignment_table_sequence, end='\n', file=out_file)

            # line 3 template read
            print(self.template_read, end='\n', file=out_file)

            # line 4 tempalte strand map
            for _ in self.template_strand_event_map:
                print(_, end=' ', file=out_file)
            print("", end="\n", file=out_file)

            # line 5 complement read
            print(self.complement_read, end='\n', file=out_file)

            # line 6 complement strand map
            for _ in self.complement_strand_event_map:
                print(_, end=' ', file=out_file)
            print("", end="\n", file=out_file)

            # line 7 template 2D event map
            for _ in self.template_event_map:
                print(_, end=' ', file=out_file)
            print("", end="\n", file=out_file)

            # line 8 template events
            template_start_time = self.template_events[0]['start']
            for mean, stdev, length, start in self.template_events['mean', 'stdv', 'length', 'start']:
                print(mean, stdev, length, (start - template_start_time), sep=' ', end=' ', file=out_file)
            print("", end="\n", file=out_file)

            # line 9 complement 2D event map
            for _ in self.complement_event_map[::-1]:
                print(_, end=' ', file=out_file)
            print("", end="\n", file=out_file)

            # line 10 complement events
            complement_start_time = self.complement_events[0]['start']
            for mean, stdev, length, start in self.complement_events['mean', 'stdv', 'length', 'start']:
                print(mean, stdev, length, (start - complement_start_time), sep=' ', end=' ', file=out_file)
            print("", end="\n", file=out_file)

            # line 11 model_state (template)
            for _ in self.template_events['model_state']:
                print(_, sep=' ', end=' ', file=out_file)
            print("", end="\n", file=out_file)

            # line 12 p(model) (template)
            for _ in self.template_events['p_model_state']:
                print(_, sep=' ', end=' ', file=out_file)
            print("", end="\n", file=out_file)

            # line 13 model_state (complement)
            for _ in self.complement_events['model_state']:
                print(_, sep=' ', end=' ', file=out_file)
            print("", end="\n", file=out_file)

            # line 14 p(model) (complement)
            for _ in self.complement_events['p_model_state']:
                print(_, sep=' ', end=' ', file=out_file)
            print("", end="\n", file=out_file)

            return True
        else:
            print("write_npRead: proceed was False", file=sys.stderr)
            return False

    def close(self):
        self.fastFive.close()


class SignalAlignment(object):
    def __init__(self, in_fast5, forward_reference, backward_reference, path_to_EC_refs, destination, stateMachineType,
                 banded, bwa_index, in_templateHmm, in_complementHmm, in_templateHdp, in_complementHdp,
                 threshold, diagonal_expansion, constraint_trim, degenerate,
                 target_regions=None, sparse_output=False):
        self.in_fast5 = in_fast5  # fast5 file to align
        self.forward_reference = forward_reference  # forward 'FASTA-oriented' reference
        self.backward_reference = backward_reference  # complement of the forward reference
        self.path_to_EC_refs = path_to_EC_refs  # place where the reference sequence with ambiguous characters is
        self.destination = destination  # place where the alignments go, should already exist
        self.stateMachineType = stateMachineType  # flag for signalMachine
        self.banded = banded  # use banded or not
        self.bwa_index = bwa_index  # index of reference sequence
        self.threshold = threshold  # min posterior probability to keep
        self.diagonal_expansion = diagonal_expansion  # alignment algorithm param
        self.constraint_trim = constraint_trim  # alignment algorithm param
        self.target_regions = target_regions  # only signal-align reads that map to these positions
        self.sparse_output = sparse_output  # smaller output files
        self.degenerate = degenerate  # set of nucleotides for degenerate characters

        # if we're using an input hmm, make sure it exists
        if (in_templateHmm is not None) and os.path.isfile(in_templateHmm):
            self.in_templateHmm = in_templateHmm
        else:
            self.in_templateHmm = None
        if (in_complementHmm is not None) and os.path.isfile(in_complementHmm):
            self.in_complementHmm = in_complementHmm
        else:
            self.in_complementHmm = None

        # similarly for HDPs
        if (in_templateHdp is not None) and os.path.isfile(in_templateHdp):
            self.in_templateHdp = in_templateHdp
        else:
            self.in_templateHdp = None
        if (in_complementHdp is not None) and os.path.isfile(in_complementHdp):
            self.in_complementHdp = in_complementHdp
        else:
            self.in_complementHdp = None

    def run(self, get_expectations=False):
        if get_expectations:
            assert self.in_templateHmm is not None and self.in_complementHmm is not None, "Need HMM files for model " \
                                                                                          "training"
        # file checks
        if os.path.isfile(self.in_fast5) is False:
            print("signalAlign - problem with file path {file}".format(file=self.in_fast5))
            return False

        # containers and defaults
        read_label = self.in_fast5.split("/")[-1]      # used in the posteriors file as identifier
        read_name = self.in_fast5.split("/")[-1][:-6]  # get the name without the '.fast5'

        # object for handling temporary files
        temp_folder = FolderHandler()
        temp_dir_path = temp_folder.open_folder(self.destination + "tempFiles_{readLabel}".format(readLabel=read_label))

        # read-specific files, could be removed later but are kept right now to make it easier to rerun commands
        temp_np_read = temp_folder.add_file_path("temp_{read}.npRead".format(read=read_label))
        temp_2d_read = temp_folder.add_file_path("temp_2Dseq_{read}.fa".format(read=read_label))

        # make the npRead and fasta
        success, version, pop1_complement = get_npRead_2dseq_and_models(fast5=self.in_fast5,
                                                                        npRead_path=temp_np_read,
                                                                        twod_read_path=temp_2d_read)

        if success is False:
            print("file {file} does not have 2D or is corrupt".format(file=read_label), file=sys.stderr)
            return False

        # add an indicator for the model being used
        if self.stateMachineType == "threeState":
            model_label = ".sm"
            stateMachineType_flag = ""
        elif self.stateMachineType == "threeStateHdp":
            model_label = ".sm3Hdp"
            stateMachineType_flag = "--sm3Hdp "
            assert (self.in_templateHdp is not None) and (self.in_complementHdp is not None), "Need to provide HDPs"
        else:  # make invalid stateMachine control?
            model_label = ".sm"
            stateMachineType_flag = ""

        # get orientation and cigar from BWA this serves as the guide alignment
        cigar_string, strand = exonerated_bwa(bwa_index=self.bwa_index, query=temp_2d_read,
                                              target_regions=self.target_regions)

        # this gives the format: /directory/for/files/file.model.orientation.tsv
        posteriors_file_path = ''

        # forward strand
        if strand == "+":
            forward = True
            posteriors_file_path = self.destination + read_name + model_label + ".forward.tsv"

        # backward strand
        if strand == "-":
            forward = False
            posteriors_file_path = self.destination + read_name + model_label + ".backward.tsv"

        # didn't map
        elif (strand != "+") and (strand != "-"):
            print("signalAlign - {read} didn't map got flag: {flag}".format(read=read_label, flag=strand),
                  file=sys.stderr)
            temp_folder.remove_folder()
            return False

        # Alignment/Expectations routine

        # containers and defaults
        path_to_signalAlign = "./signalMachine"

        # flags

        # input (match) models
        if self.in_templateHmm is None:
            self.in_templateHmm = default_template_model_from_version(version=version)
        if self.in_complementHmm is None:
            self.in_complementHmm = default_complement_model_from_version(version=version,
                                                                          pop1_complement=pop1_complement)

        assert self.in_templateHmm is not None and self.in_complementHmm is not None
        template_model_flag = "-T {} ".format(self.in_templateHmm)
        complement_model_flag = "-C {} ".format(self.in_complementHmm)
        print("signalAlign - NOTICE: template model {t} complement model {c}\n"
              "".format(t=self.in_templateHmm, c=self.in_complementHmm), file=sys.stderr)

        # reference sequences
        if self.forward_reference is not None or self.backward_reference is not None:
            assert self.forward_reference is not None and self.backward_reference is not None, \
                "Need forward and backward reference sequences"
            forward_ref_flag = "-f {f_ref} ".format(f_ref=self.forward_reference)
            backward_ref_flag = "-b {b_ref} ".format(b_ref=self.backward_reference)
            error_correct_ref_path = ""
        else:
            assert self.forward_reference is None and self.backward_reference is None, \
                "Erroneously gave forward and backward references when trying to do error correction"
            assert self.path_to_EC_refs is not None, "Need to provide path to ambiguous reference sequences"
            forward_ref_flag = ""
            backward_ref_flag = ""
            error_correct_ref_path = "-p {path}".format(path=self.path_to_EC_refs)

        # input HDPs
        if (self.in_templateHdp is not None) and (self.in_complementHdp is not None):
            hdp_flags = "-v {tHdp_loc} -w {cHdp_loc} ".format(tHdp_loc=self.in_templateHdp,
                                                              cHdp_loc=self.in_complementHdp)
        else:
            hdp_flags = ""

        # threshold
        if self.threshold is not None:
            threshold_flag = "-D {threshold} ".format(threshold=self.threshold)
        else:
            threshold_flag = ""

        # diagonal expansion
        if self.diagonal_expansion is not None:
            diag_expansion_flag = "-x {expansion} ".format(expansion=self.diagonal_expansion)
        else:
            diag_expansion_flag = ""

        # constraint trim
        if self.constraint_trim is not None:
            trim_flag = "-m {trim} ".format(trim=self.constraint_trim)
        else:
            trim_flag = ""

        # banded alignment
        if self.banded is True:
            pass
        else:
            trim_flag = "-m 9999"

        # sparse output
        if self.sparse_output is True:
            sparse_flag = "--sparse_output "
        else:
            sparse_flag = ""

        # degenerate nucleotide information
        if self.degenerate is not None:
            degenerate_flag = "-o {} ".format(self.degenerate)
        else:
            degenerate_flag = ""
        # commands
        if get_expectations:
            template_expectations_file_path = self.destination + read_name + ".template.expectations"
            complement_expectations_file_path = self.destination + read_name + ".complement.expectations"

            command = \
                "echo {cigar} | {vA} {degen}{sparse}{model}{f_ref}{b_ref}{eC} -q {npRead} " \
                "{t_model}{c_model}{thresh}{expansion}{trim} {hdp}-L {readLabel} " \
                "-t {templateExpectations} -c {complementExpectations}"\
                .format(cigar=cigar_string, vA=path_to_signalAlign, model=stateMachineType_flag,
                        f_ref=forward_ref_flag, b_ref=backward_ref_flag,
                        npRead=temp_np_read, readLabel=read_label,
                        templateExpectations=template_expectations_file_path, hdp=hdp_flags,
                        complementExpectations=complement_expectations_file_path, t_model=template_model_flag,
                        c_model=complement_model_flag, thresh=threshold_flag, expansion=diag_expansion_flag,
                        trim=trim_flag, degen=degenerate_flag, sparse=sparse_flag, eC=error_correct_ref_path)
        else:
            command = \
                "echo {cigar} | {vA} {degen}{sparse}{model}{f_ref}{b_ref}{eC} -q {npRead} " \
                "{t_model}{c_model}{thresh}{expansion}{trim} " \
                "-u {posteriors} {hdp}-L {readLabel}"\
                .format(cigar=cigar_string, vA=path_to_signalAlign, model=stateMachineType_flag, sparse=sparse_flag,
                        f_ref=forward_ref_flag, b_ref=backward_ref_flag, eC=error_correct_ref_path,
                        readLabel=read_label, npRead=temp_np_read,
                        t_model=template_model_flag, c_model=complement_model_flag,
                        posteriors=posteriors_file_path, thresh=threshold_flag, expansion=diag_expansion_flag,
                        trim=trim_flag, hdp=hdp_flags, degen=degenerate_flag)

        # run
        print("signalAlign - running command: ", command, end="\n", file=sys.stderr)
        os.system(command)
        temp_folder.remove_folder()
        return True


class SignalHmm(object):
    def __init__(self, model_type):
        self.match_model_params = 5  # level_mean, level_sd, noise_mean, noise_sd, noise_lambda
        self.model_type = model_type  # ID of model type
        self.state_number = {"threeState": 3, "threeStateHdp": 3}[model_type]
        self.symbol_set_size = 0
        self.transitions = np.zeros(self.state_number**2)
        self.transitions_expectations = np.zeros(self.state_number**2)
        self.likelihood = 0.0
        self.running_likelihoods = []
        self.alphabet_size = 0
        self.alphabet = ""
        self.kmer_length = 0
        self.has_model = False
        self.normalized = False

        # event model for describing normal distributions for each kmer
        self.event_model = {"means": np.zeros(self.symbol_set_size),
                            "SDs": np.zeros(self.symbol_set_size),
                            "noise_means": np.zeros(self.symbol_set_size),
                            "noise_SDs": np.zeros(self.symbol_set_size),
                            "noise_lambdas": np.zeros(self.symbol_set_size)}

    def normalize_transitions_expectations(self):
        # normalize transitions
        for from_state in xrange(self.state_number):
            i = self.state_number * from_state
            j = sum(self.transitions_expectations[i:i+self.state_number])
            for to_state in xrange(self.state_number):
                self.transitions_expectations[i + to_state] = self.transitions_expectations[i + to_state] / j

    def set_default_transitions(self):
        MATCH_CONTINUE = np.exp(-0.23552123624314988)     # stride
        MATCH_FROM_GAP_X = np.exp(-0.21880828092192281)   # 1 - skip'
        MATCH_FROM_GAP_Y = np.exp(-0.013406326748077823)  # 1 - (skip + stay)
        GAP_OPEN_X = np.exp(-1.6269694202638481)          # skip
        GAP_OPEN_Y = np.exp(-4.3187242127300092)          # 1 - (skip + stride)
        GAP_EXTEND_X = np.exp(-1.6269694202638481)        # skip'
        GAP_EXTEND_Y = np.exp(-4.3187242127239411)        # stay (1 - (skip + stay))
        GAP_SWITCH_TO_X = 0.000000001
        GAP_SWITCH_TO_Y = 0.0
        self.transitions = [
            MATCH_CONTINUE, GAP_OPEN_X, GAP_OPEN_Y,
            MATCH_FROM_GAP_X, GAP_EXTEND_X, GAP_SWITCH_TO_Y,
            MATCH_FROM_GAP_Y, GAP_SWITCH_TO_X, GAP_EXTEND_Y
        ]
        return

    def load_model(self, model_file):
        # the model file has the format:
        # line 0: stateNumber \t alphabetSize \t alphabet \t kmerLength
        # line 1: match->match \t match->gapX \t match->gapY \t
        #         gapX->match \t gapX->gapX \t gapX->gapY \t
        #         gapY->match \t gapY->gapX \t gapY->gapY \n
        # line 2: [level_mean] [level_sd] [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n
        assert os.path.exists(model_file), "signalHmm.load_model - didn't find model here{}?".format(model_file)

        fH = open(model_file, 'r')

        line = fH.readline().split()
        # check for correct header length
        assert len(line) == 4, "signalHmm.load_model - incorrect line length line:{}".format(''.join(line))
        # check stateNumber
        assert int(line[0]) == self.state_number, "signalHmm.load_model - incorrect stateNumber got {got} should be {exp}" \
                                                  "".format(got=int(line[0]), exp=self.state_number)
        # load model parameters
        self.alphabet_size = int(line[1])
        self.alphabet = line[2]
        self.kmer_length = int(line[3])
        self.symbol_set_size = self.alphabet_size**self.kmer_length
        assert self.symbol_set_size > 0, "signalHmm.load_model - Got 0 for symbol_set_size"
        assert self.symbol_set_size <= 6**6, "signalHmm.load_model - Got more than 6^6 for symbol_set_size got {}" \
                                             "".format(self.symbol_set_size)

        line = map(float, fH.readline().split())
        assert len(line) == len(self.transitions) + 1, "signalHmm.load_model incorrect transitions line"
        self.transitions = line[:-1]
        self.likelihood = line[-1]

        line = map(float, fH.readline().split())
        assert len(line) == self.symbol_set_size * NB_MODEL_PARAMS, \
            "signalHmm.load_model incorrect event model line"
        self.event_model["means"] = line[::NB_MODEL_PARAMS]
        self.event_model["SDs"] = line[1::NB_MODEL_PARAMS]
        self.event_model["noise_means"] = line[2::NB_MODEL_PARAMS]
        self.event_model["noise_SDs"] = line[3::NB_MODEL_PARAMS]
        self.event_model["noise_lambdas"] = line[4::NB_MODEL_PARAMS]

        assert not np.any(self.event_model["means"] == 0.0), "signalHmm.load_model, this model has 0 E_means"
        assert not np.any(self.event_model["SDs"] == 0.0), "signalHmm.load_model, this model has 0 E_means"
        assert not np.any(self.event_model["noise_means"] == 0.0), "signalHmm.load_model, this model has 0 E_noise_means"
        assert not np.any(self.event_model["noise_SDs"] == 0.0), "signalHmm.load_model, this model has 0 E_noise_SDs"
        self.has_model = True

    def write(self, out_file):
        # the model file has the format:
        # line 0: stateNumber \t alphabetSize \t alphabet \t kmerLength
        # line 1: match->match \t match->gapX \t match->gapY \t
        #         gapX->match \t gapX->gapX \t gapX->gapY \t
        #         gapY->match \t gapY->gapX \t gapY->gapY \n
        # line 2: [level_mean] [level_sd] [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n
        assert self.has_model, "Shouldn't be writing down a Hmm that has no Model"
        assert self.normalized, "Shouldn't be writing down a not normalized HMM"

        f = open(out_file, 'w')

        # line 0
        f.write("{stateNumber}\t{alphabetSize}\t{alphabet}\t{kmerLength}\n"
                "".format(stateNumber=self.state_number, alphabetSize=self.alphabet_size,
                          alphabet=self.alphabet, kmerLength=self.kmer_length))
        # line 1 transitions
        for i in xrange(self.state_number * self.state_number):
            f.write("{transition}\t".format(transition=str(self.transitions[i])))
        # likelihood
        f.write("{}\n".format(str(self.likelihood)))

        # line 2 Event Model
        for k in xrange(self.symbol_set_size):
            f.write("{level_mean}\t{level_sd}\t{noise_mean}\t{noise_sd}\t{noise_lambda}\t"
                    "".format(level_mean=self.event_model["means"][k], level_sd=self.event_model["SDs"][k],
                              noise_mean=self.event_model["noise_means"][k], noise_sd=self.event_model["noise_SDs"][k],
                              noise_lambda=self.event_model["noise_lambdas"][k]))
        f.write("\n")

        f.close()


class ContinuousPairHmm(SignalHmm):
    def __init__(self, model_type):
        super(ContinuousPairHmm, self).__init__(model_type=model_type)
        self.set_default_transitions()

        # bins for expectations
        self.mean_expectations = np.zeros(self.symbol_set_size)
        self.sd_expectations = np.zeros(self.symbol_set_size)
        self.posteriors = np.zeros(self.symbol_set_size)
        self.observed = np.zeros(self.symbol_set_size, dtype=bool)
        self.has_model = False
        self.normalized = False

    def add_expectations_file(self, expectations_file):
        # expectations files have the format:
        # line 0: stateNumber \t alphabetSize \t alphabet \t kmerLength
        # line 1: match->match \t match->gapX \t match->gapY \t
        #         gapX->match \t gapX->gapX \t gapX->gapY \t
        #         gapY->match \t gapY->gapX \t gapY->gapY \n
        # line 2: [level_mean] [level_sd] [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n
        # line 3: event expectations [mean] [sd] / kmer \n
        # line 4: posteriors 1 per kmer \n
        # line 5: observed 1 per kmer \n
        if not os.path.exists(expectations_file) or os.stat(expectations_file).st_size == 0:
            print("Empty or missing file {}".format(expectations_file))
            return False

        fH = open(expectations_file, 'r')

        # line 0
        line = fH.readline().split()
        if len(line) != 4:
            print("cpHMM: check_file - incorrect header (param line): {}".format(expectations_file), file=sys.stderr)
            fH.close()
            return False
        if int(line[0]) != self.state_number:
            print("cpHMM: state number error should be {exp} got {obs}"
                  "".format(exp=self.state_number, obs=line[0]), file=sys.stderr)
            fH.close()
            return False
        if int(line[1]) != self.alphabet_size:
            print("cpHMM: alphabet size error - incorrect parameters: {file}, line {line}"
                  "".format(file=expectations_file, line=''.join(line)), file=sys.stderr)
            fH.close()
            return False
        if line[2] != self.alphabet:
            print("cpHMM: alphabet error - incorrect parameters: {file}, line {line}"
                  "".format(file=expectations_file, line=''.join(line)), file=sys.stderr)
            fH.close()
            return False
        if int(line[3]) != self.kmer_length:
            print("cpHMM: kmer length error - incorrect parameters: {file}, line {line}"
                  "".format(file=expectations_file, line=''.join(line)), file=sys.stderr)
            fH.close()
            return False

        # line 1: transitions, likelihood
        line = map(float, fH.readline().split())
        # check if valid
        if len(line) != (len(self.transitions) + 1):
            print("cpHMM: check_file - bad file (transitions expectations): {}".format(expectations_file),
                  file=sys.stderr)
            fH.close()
            return False

        self.likelihood += line[-1]
        self.transitions_expectations = map(lambda x: sum(x), zip(self.transitions_expectations, line[0:-1]))

        # line 2: event model
        line = map(float, fH.readline().split())
        if len(line) != self.symbol_set_size * NB_MODEL_PARAMS:
            print("cpHMM: check_file - bad file (event model): {}".format(expectations_file), file=sys.stderr)
            fH.close()
            return False

        # line 3 event expectations [E_mean, E_sd]
        line = map(float, fH.readline().split())
        if len(line) != self.symbol_set_size * NORM_DIST_PARAMS:
            print("cpHMM: check_file - bad file (event expectations): {}".format(expectations_file), file=sys.stderr)
            fH.close()
            return False

        self.mean_expectations = [i + j for i, j in izip(self.mean_expectations, line[::NORM_DIST_PARAMS])]
        self.sd_expectations = [i + j for i, j in izip(self.sd_expectations, line[1::NORM_DIST_PARAMS])]

        # line 4, posteriors
        line = map(float, fH.readline().split())
        if len(line) != self.symbol_set_size:
            print("cpHMM: check_file - bad file (posteriors): {}".format(expectations_file), file=sys.stderr)
            fH.close()
            return False

        self.posteriors = map(lambda x: sum(x), zip(self.posteriors, line))

        line = map(bool, fH.readline().split())
        if len(line) != self.symbol_set_size:
            print("cpHMM: check_file - bad file (observations): {}".format(expectations_file), file=sys.stderr)
            fH.close()
            return False

        self.observed = [any(b) for b in zip(self.observed, line)]

        fH.close()
        return True

    def normalize(self, update_transitions, update_emissions):
        # normalize transitions expectations
        self.normalize_transitions_expectations()

        # update
        if update_transitions is True:
            for i in xrange(self.state_number**2):
                self.transitions[i] = self.transitions_expectations[i]

        # calculate the new expected mean and standard deviation for the kmer normal distributions
        if update_emissions:
            for k in xrange(self.symbol_set_size):  # TODO implement learning rate
                if self.observed[k] is True:
                    u_k = self.mean_expectations[k] / self.posteriors[k]
                    o_k = np.sqrt(self.sd_expectations[k] / self.posteriors[k])
                    if u_k > 0:
                        self.event_model["means"][k] = u_k
                        self.event_model["SDs"][k] = o_k
                else:
                    continue
        self.normalized = True


class HdpSignalHmm(SignalHmm):
    def __init__(self, model_type, threshold):
        super(HdpSignalHmm, self).__init__(model_type=model_type)
        self.set_default_transitions()
        self.threshold = threshold
        self.kmer_assignments = []
        self.event_assignments = []
        self.assignments_record = []

    def add_expectations_file(self, expectations_file):
        # expectations files have the format:
        # line 0: stateNumber \t alphabetSize \t alphabet \t kmerLength
        # line 1: match->match \t match->gapX \t match->gapY \t
        #         gapX->match \t gapX->gapX \t gapX->gapY \t
        #         gapY->match \t gapY->gapX \t gapY->gapY \n
        # line 2: [level_mean] [level_sd] [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n
        # line 3: event assignments
        # line 4: kmer assignments
        if not os.path.exists(expectations_file) or os.stat(expectations_file).st_size == 0:
            print("Empty or missing file {}".format(expectations_file))
            return

        fH = open(expectations_file, 'r')

        # line 0: smType stateNumber, symbolSetSize
        # line 0
        line = fH.readline().split()
        if len(line) != 4:
            print("cpHMM: check_file - bad file (param line): {}".format(expectations_file), file=sys.stderr)
            fH.close()
            return
        if line[0] != self.state_number or \
                        int(line[1]) != self.alphabet_size or \
                        line[2] != self.alphabet or \
                        int(line[3]) != self.kmer_length:
            print("cpHMM: check_file - bad file (Hmm params): {}".format(expectations_file), file=sys.stderr)
            fH.close()
            return

        # line 1: transitions, likelihood
        line = map(float, fH.readline().split())

        # check if valid file
        if len(line) != (len(self.transitions) + 1):
            print("PYSENTINAL - problem with file {}".format(expectations_file), file=sys.stdout)
            fH.close()
            return

        self.likelihood += line[-1]
        self.transitions_expectations = map(lambda x: sum(x), zip(self.transitions_expectations, line[0:-1]))

        # line 3: event assignments
        line = map(float, fH.readline().split())
        self.event_assignments += line

        # line 4: kmer assignments
        line = map(str, fH.readline().split())
        self.kmer_assignments += line

        fH.close()

        #assert (len(self.kmer_assignments) == self.number_of_assignments) and \
        #       (len(self.event_assignments) == self.number_of_assignments), \
        #    "trainModels - add_expectations_file: invalid number of assignments"

    def reset_assignments(self):
        self.assignments_record.append(len(self.event_assignments))
        self.event_assignments = []
        self.kmer_assignments = []

    def normalize(self, update_transitions, update_emissions=None):
        self.normalize_transitions_expectations()
        if update_transitions is True:
            for i in xrange(self.state_number**2):
                self.transitions[i] = self.transitions_expectations[i]
        self.normalized = True
