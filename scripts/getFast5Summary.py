
from __future__ import print_function
import sys
import os
from argparse import ArgumentParser
import glob
from signalAlignLib import NanoporeRead


def parse_args():
    parser = ArgumentParser(description=__doc__)

    parser.add_argument('--fast5_glob', '-g', action='store',
                        dest='fast5_glob', required=False, type=str, default="*.fast5",
                        help="glob for matching fast5 files")
    parser.add_argument('--desired_read_ids', '-r', action='store',
                        dest='desired_read_ids', required=False, type=str, default=None,
                        help="file with read ids (one per line). if set, only these will be output")
    parser.add_argument('--output_read_id', '-R', action='store_true', dest='output_read_id', default=False,
                        help='output read id')
    parser.add_argument('--output_fasta', '-A', action='store_true', dest='output_fasta', default=False,
                        help='output fasta from read')
    parser.add_argument('--no_fast5_path', '-5', action='store_false', dest="output_fast5", default=True,
                        help='remove fast5 path from output')

    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    if not (args.output_read_id or args.output_fasta):
        print("No output parameters specified; will output read_id", file=sys.stderr)
        args.output_read_id = True

    if args.desired_read_ids is not None:
        read_ids = set()
        with open(args.desired_read_ids, 'r') as file:
            for line in file:
                read_ids.add(line.strip())
        if len(read_ids) == 0:
            raise Exception("No read ids in file {}".format(args.desired_read_ids))
        print("Found {} read ids for inclusion".format(len(read_ids)), file=sys.stderr)

    files = glob.glob(args.fast5_glob)
    print("Found {} files matching {}".format(len(files), args.fast5_glob), file=sys.stderr)
    if len(files) == 0:
        return

    header = []
    if args.output_fast5: header.append("FAST5")
    if args.output_read_id: header.append("READ_ID")
    if args.output_fasta: header.append("FASTA")
    header = "\t".join(header)
    print("#" + header)

    output_count = 0
    for file in files:
        file = os.path.abspath(file)
        try:
            fast5 = NanoporeRead(file)
            if args.desired_read_ids is not None and fast5.read_label not in read_ids:
                    continue
            row = []
            if args.output_fast5: row.append(file)
            if args.output_read_id: row.append(fast5.read_label)
            if args.output_fasta: row.append(fast5.template_read if len(fast5.template_read) != 0 else fast5.complement_read)
            print("\t".join(row))
            output_count += 1
        except Exception, e:
            print("Error reading file {}: {}".format(file, e), file=sys.stderr)

    if args.desired_read_ids is not None:
        print("Found {} / {} desired read ids".format(output_count, len(read_ids)), file=sys.stderr)
    else:
        print("Found {} read ids".format(output_count), file=sys.stderr)


if __name__ == "__main__":
    main()


