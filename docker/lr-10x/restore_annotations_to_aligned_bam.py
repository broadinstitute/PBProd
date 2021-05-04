import argparse
import sys
import time
import os

import tqdm
import pysam


def main(bam_name, aligned_bam_name, out_bam_name):

    t_start = time.time()

    tag_names = []
    name_taglist_dict = dict()
    pysam.set_verbosity(0)  # silence message about the .bai file not being found

    print("Verifying input files exist...")
    files_ok = True
    for f in [bam_name, aligned_bam_name]:
        if not os.path.exists(f):
            print(f"ERROR: Input file does not exist: {f}")
            files_ok = False
    if not files_ok:
        sys.exit(1)
    print("Input files verified.")

    print(f"Reading in tags from: {bam_name}")
    with pysam.AlignmentFile(
            bam_name, "rb", check_sq=False, require_index=False
    ) as bam_file, tqdm.tqdm(
        desc="Progress",
        unit=" read",
        file=sys.stderr
    ) as pbar:
        for read in bam_file:
            tags_to_keep = []
            if len(tag_names) == 0:
                # Keep all tags:
                for t in read.get_tags():
                    tags_to_keep.append((t[0], t[1]))
            else:
                # Keep only the tags we specified:
                for tag_name in tag_names:
                    try:
                        tags_to_keep.append((tag_name, read.get_tag(tag_name)))
                    except KeyError:
                        pass
            name_taglist_dict[read.query_name] = tags_to_keep
            pbar.update(1)

    print("Tag list assembled.")

    print(f"Applying tags to reads in {aligned_bam_name} and writing to {out_bam_name}")
    with pysam.AlignmentFile(
            aligned_bam_name, "rb", check_sq=False, require_index=False
    ) as aligned_bam_file, tqdm.tqdm(
        desc="Progress",
        unit=" read",
        file=sys.stderr
    ) as pbar:

        # Get our header from the input bam file:
        out_bam_header_dict = aligned_bam_file.header.to_dict()

        # Add our program group to it:
        pg_dict = {
            "ID": f"restore-annotations-0.0.1",
            "PN": "restore-annotations",
            "VN": f"0.0.1",
            "DS": 'Reads an unaligned bam with annotations and an aligned bam containing the same reads.  '
                  'Copies tags on the unaligned bam to reads in the aligned file with the same name.',
            "CL": " ".join(sys.argv),
        }
        if "PG" in out_bam_header_dict:
            out_bam_header_dict["PG"].append(pg_dict)
        else:
            out_bam_header_dict["PG"] = [pg_dict]
        out_header = pysam.AlignmentHeader.from_dict(out_bam_header_dict)

        with pysam.AlignmentFile(out_bam_name, "wb", header=out_header) as out_bam_file:
            for read in aligned_bam_file:
                if read.query_name in name_taglist_dict:
                    for tag in name_taglist_dict[read.query_name]:
                        read.set_tag(tag[0], tag[1])

                out_bam_file.write(read)

                pbar.update(1)
    t_end = time.time()
    print("Done")
    print(f"Elapsed time: {t_end - t_start:2.2f}s")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Reads an unaligned bam with annotations and an aligned bam containing the same reads.'
                    'Copies tags on the unaligned bam to reads in the aligned file with the same name.')
    required_named_args = parser.add_argument_group('required named arguments')
    required_named_args.add_argument('-b', '--bam', help='Bam file from which to read annotations.', required=True)
    required_named_args.add_argument('-a', '--aligned-bam', help='Aligned bam to which to apply annotations.',
                                     required=True)
    required_named_args.add_argument('-o', '--out-name',
                                     help='Output bam file name', required=True)

    args = parser.parse_args()

    main(args.bam, args.aligned_bam, args.out_name)
