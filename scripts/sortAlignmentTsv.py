
from __future__ import print_function
import sys
import os
from argparse import ArgumentParser
import glob
from signalAlignLib import NanoporeRead


def parse_args():
    parser = ArgumentParser(description=__doc__)

    parser.add_argument('--alignment_glob', '-a', action='store',
                        dest='alignment_glob', required=False, type=str, default="*.tsv",
                        help="glob for matching alignment files")
    parser.add_argument('--sort_column_idx', '-c', action='store',
                        dest='sort_column_idx', required=False, type=int, default=1,
                        help="0-based index to sort off of (default is 1, the position column)")
    parser.add_argument('--sort_column_string', required=False, default=False, action='store_true', dest="sort_column_string",
                        help="the sort column should be treated as a string (float is default)")
    parser.add_argument('--output_fmt', '-o', action='store', dest='output_fmt', default="{basename}.sorted.{ext}",
                        required=False, help='string describing output fmt: {basename} is name without extension, '
                                             '{filename} is whole filename, {ext} is extension, {dir} is directory.  '
                                             'default: "{dir}/{basename}.sorted.{ext}"' )

    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    files = glob.glob(args.alignment_glob)
    print("Found {} files matching {}".format(len(files), args.alignment_glob), file=sys.stderr)
    if len(files) == 0:
        return

    for file in files:
        file = os.path.abspath(file)
        filename = os.path.basename(file)
        dir = os.path.dirname(file)
        basename = os.path.splitext(filename)[0]
        ext = os.path.splitext(filename)[1][1:]

        lines = list()
        with open(file, 'r') as input:
            for line in input:
                if len(line.strip()) == 0 or line.startswith("#"): continue
                lines.append(line.strip().split())

        type_function = str if args.sort_column_string else float
        lines.sort(key=lambda x: type_function(x[args.sort_column_idx]))

        output_file = args.output_fmt.format(basename=basename, filename=filename, ext=ext, dir=dir)
        with open(output_file, 'w') as output:
            for line in lines:
                output.write("\t".join(line) + "\n")

        print("Sorted {} lines from {} into {}".format(len(lines), file, output_file), file=sys.stderr)


if __name__ == "__main__":
    main()


