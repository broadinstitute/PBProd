#!/usr/bin/env python

import os
import sys
import argparse
import csv
import pickle

import scipy
import anndata
import numpy as np
import pandas as pd

from tqdm import tqdm

TX_ENTRY_STRING = "transcript"

CONTIG_FIELD = "contig"
START_FIELD = "start"
END_FIELD = "end"

GENE_ID_FIELD = "gene_id"
TX_ID_FIELD = "transcript_id"

GENCODE_GENE_NAME_FIELD = "gene_name"
GENCODE_TX_NAME_FIELD = "transcript_name"

STRINGTIE_GENE_ID_FIELD = "ref_gene_id"
STRINGTIE_TX_ID_FIELD = "reference_id"
STRINGTIE_GENE_NAME_FIELD = "ref_gene_name"

ALT_NAME_SUFFIX = "_PAR_Y"


def get_gtf_field_val_dict(gtf_file, force_rebuild=False):

    global TX_ENTRY_STRING
    global GENCODE_GENE_NAME_FIELD
    global GENE_ID_FIELD
    global GENCODE_TX_NAME_FIELD
    global TX_ID_FIELD
    global ALT_NAME_SUFFIX

    global CONTIG_FIELD
    global START_FIELD
    global END_FIELD

    gtf_dict = dict()

    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"GTF file does not exist: {gtf_file}")

    pickle_file_name = os.path.splitext(os.path.basename(gtf_file))[0] + ".big_map.pickle"

    if not force_rebuild and os.path.exists(pickle_file_name):
        print(f"Loading GTF field map from {pickle_file_name}...", end="\t", file=sys.stderr)
        gtf_dict = pickle.load(open(pickle_file_name, "rb"))
        print("Done!", file=sys.stderr)
    else:
        print(f"Creating master GTF field map keyed by transcript ID using {gtf_file}...",
              flush=True, file=sys.stderr)

        with open(gtf_file, "r") as f, tqdm(desc="Processing GTF File", unit=" line") as pbar:
            tsv_file = csv.reader(f, delimiter="\t")
            for row in tsv_file:
                # Ignore comments and make sure we only process transcript entries:
                if not row[0].startswith("#") and row[2] == TX_ENTRY_STRING:

                    # Parse the row data into a dict and make sure we're only looking at protein coding rows:
                    row_data_dict = {
                        field.strip().split(" ")[0].replace('"', ""): field.strip().split(" ")[1].replace('"', "") for
                        field in row[8].split(";") if len(field) != 0
                    }
                    row_data_dict[CONTIG_FIELD] = row[0]
                    row_data_dict[START_FIELD] = int(row[4])
                    row_data_dict[END_FIELD] = int(row[5])

                    # Make sure our names are unique:
                    if row_data_dict[TX_ID_FIELD].endswith(ALT_NAME_SUFFIX):
                        row_data_dict[GENCODE_TX_NAME_FIELD] = row_data_dict[GENCODE_TX_NAME_FIELD] + ALT_NAME_SUFFIX
                        row_data_dict[GENCODE_GENE_NAME_FIELD] = row_data_dict[GENCODE_GENE_NAME_FIELD] + ALT_NAME_SUFFIX

                    # Add this row to our dict keyed by transcript ID:
                    gtf_dict[row_data_dict[TX_ID_FIELD]] = row_data_dict

                pbar.update(1)

        print("Pickling data...", file=sys.stderr)
        pickle.dump(gtf_dict, open(pickle_file_name, "wb"))
        print("Done!", file=sys.stderr)

    return gtf_dict


def intervals_overlap(contig, start, end, contig2, start2, end2):
    if contig == contig2:
        if (start <= start2 <= end) or (start <= end2 <= end) or (start2 <= start and end2 >= end):
            return True
    return False


def interval_overlaps_any_in_interval_list(contig, start, end, interval_list):
    for contig2, start2, end2 in interval_list:
        if intervals_overlap(contig, start, end, contig2, start2, end2):
            return True
    return False


def create_combined_anndata(input_tsv, gtf_field_dict, overlap_intervals=None,
                            overlap_intervals_label="overlaps_intervals_of_interest", force_recount=False):

    """Create an anndata object holding the given gene/transcript information.
    NOTE: This MUST be a sparse matrix - we got lots of data here.

    First convert the counts to a matrix that looks like what scanpy expects:
    columns = Transcripts name / Transcript ID / Gene Name / Gene ID (variables)
    Rows = cell barcodes (observations)
    data = counts"""

    global GENCODE_GENE_NAME_FIELD
    global GENE_ID_FIELD
    global GENCODE_TX_NAME_FIELD
    global TX_ID_FIELD

    global STRINGTIE_GENE_ID_FIELD
    global STRINGTIE_TX_ID_FIELD
    global STRINGTIE_GENE_NAME_FIELD

    global CONTIG_FIELD
    global START_FIELD
    global END_FIELD

    cell_barcode_to_tx_to_umi_dict = dict()
    cell_barcode_to_tx_count_dict = dict()

    pickle_file_name = os.path.splitext(os.path.basename(input_tsv))[0] + ".tx_raw_count_matrix.pickle"

    if not force_recount and os.path.exists(pickle_file_name):
        print(f"Loading count map from {pickle_file_name}...", end="\t", file=sys.stderr)
        cell_barcode_to_tx_count_dict = pickle.load(open(pickle_file_name, "rb"))
        print("Done!", file=sys.stderr)
    else:
        # Get our cell tx counts:
        with open(input_tsv, "r") as f, tqdm(desc="Processing Raw Cell Counts", unit=" count") as pbar:
            tsv_file = csv.reader(f, delimiter="\t")
            next(tsv_file)
            for row in tsv_file:
                tx_col = row[0]
                pipe_pos = tx_col.find("|")
                if pipe_pos >= 0:
                    tx_id = tx_col[:tx_col.find("|")]
                else:
                    tx_id = tx_col

                cell_barcode = row[1]
                umi = row[2]

                # Handle the Transcript names:
                if cell_barcode not in cell_barcode_to_tx_to_umi_dict:
                    cell_barcode_to_tx_to_umi_dict[cell_barcode] = dict()
                    cell_barcode_to_tx_to_umi_dict[cell_barcode][tx_id] = set()
                    cell_barcode_to_tx_to_umi_dict[cell_barcode][tx_id].add(umi)

                    cell_barcode_to_tx_count_dict[cell_barcode] = dict()
                    cell_barcode_to_tx_count_dict[cell_barcode][tx_id] = 1

                elif tx_id not in cell_barcode_to_tx_to_umi_dict[cell_barcode]:
                    cell_barcode_to_tx_to_umi_dict[cell_barcode][tx_id] = set()
                    cell_barcode_to_tx_to_umi_dict[cell_barcode][tx_id].add(umi)

                    cell_barcode_to_tx_count_dict[cell_barcode][tx_id] = 1

                elif umi not in cell_barcode_to_tx_to_umi_dict[cell_barcode][tx_id]:
                    cell_barcode_to_tx_to_umi_dict[cell_barcode][tx_id].add(umi)

                    cell_barcode_to_tx_count_dict[cell_barcode][tx_id] += 1

                pbar.update(1)

        print("Pickling data...", file=sys.stderr)
        pickle.dump(cell_barcode_to_tx_count_dict, open(pickle_file_name, "wb"))
        print("Done!", file=sys.stderr)

    ##################################################################################################################
    # Write our cell TX counts to our adata object:

    # Create unique row / column identifiers into which to aggregate data:
    cell_barcodes = np.array(list(cell_barcode_to_tx_count_dict.keys()))
    tx_ids = np.unique(np.array(list(gtf_field_dict.keys())))

    # Populate the count matrix:
    pickle_file_name = os.path.splitext(os.path.basename(input_tsv))[0] + ".cell_transcript_count_matrix.pickle"
    if not force_recount and os.path.exists(pickle_file_name):
        print(f"Loading count map from {pickle_file_name}...", end="\t", file=sys.stderr)
        count_mat = pickle.load(open(pickle_file_name, "rb"))
        print("Done!", file=sys.stderr)
    else:
        count_mat = scipy.sparse.lil_matrix((len(cell_barcodes), len(tx_ids)), dtype=np.uint32)

        tx_id_index_dict = {name: i for i, name in enumerate(tx_ids)}
        with tqdm(desc=f"Creating cell transcript count matrix", unit=" cell",
                  total=len(cell_barcode_to_tx_count_dict)) as pbar:

            for i, (cb, counts_dict) in enumerate(cell_barcode_to_tx_count_dict.items()):
                # Put the counts for each transcript in the right indices:
                for tx_id, count in counts_dict.items():
                    count_mat[i, tx_id_index_dict[tx_id]] = count
                pbar.update(1)

        print("Pickling data...", file=sys.stderr)
        pickle.dump(count_mat, open(pickle_file_name, "wb"))
        print("Done!", file=sys.stderr)

    # Now we set up the variables that we're going to apply to each observation:
    # We'll try two different ways - one for Gencode GTFs and one for stringtie GTFs:

    # First we handle our overlapping tests:
    # Mark the transcripts that overlap our intervals:
    if overlap_intervals:
        print(f"Determining overlap intervals for the given interval list... ", end="", file=sys.stderr)
        tx_overlap_flags = [
            interval_overlaps_any_in_interval_list(
                gtf_field_dict[tx_id][CONTIG_FIELD],
                gtf_field_dict[tx_id][START_FIELD],
                gtf_field_dict[tx_id][END_FIELD],
                overlap_intervals
            ) for tx_id in tx_ids
        ]
        print("Done!", file=sys.stderr)

    is_gencode = True
    try:
        # Try Gencode first:
        transcript_names = [gtf_field_dict[tx_id][GENCODE_TX_NAME_FIELD] for tx_id in tx_ids]
        gene_names = [gtf_field_dict[tx_id][GENCODE_GENE_NAME_FIELD] for tx_id in tx_ids]
        gene_ids = [gtf_field_dict[tx_id][GENE_ID_FIELD] for tx_id in tx_ids]

        de_novo_gene_ids = ["N/A"] * len(transcript_names)
        de_novo_transcript_ids = ["N/A"] * len(transcript_names)
        is_de_novo = [False] * len(transcript_names)
        is_gene_id_ambiguous = [False] * len(transcript_names)

    except KeyError:
        # We must have stringtie data:

        de_novo_transcript_ids = tx_ids
        de_novo_gene_ids = [gtf_field_dict[tx_id][GENE_ID_FIELD] for tx_id in tx_ids]

        # TX Names for stringtie are the same as TX IDs:
        transcript_names = tx_ids

        # Get our fields for existing transcripts:
        gene_names = [gtf_field_dict[tx_id][STRINGTIE_GENE_NAME_FIELD] if STRINGTIE_GENE_NAME_FIELD in gtf_field_dict[tx_id] else gtf_field_dict[tx_id][GENE_ID_FIELD] for tx_id in tx_ids]
        gene_ids = [gtf_field_dict[tx_id][STRINGTIE_GENE_ID_FIELD] if STRINGTIE_GENE_ID_FIELD in gtf_field_dict[tx_id] else gtf_field_dict[tx_id][GENE_ID_FIELD] for tx_id in tx_ids]

        # Get our new transcript ids last so we can cue off of them earlier:
        tx_ids = [gtf_field_dict[tx_id][STRINGTIE_TX_ID_FIELD] if STRINGTIE_GENE_NAME_FIELD in gtf_field_dict[tx_id] else tx_id for tx_id in tx_ids]

        # Mark or de novo transcripts:
        is_de_novo = [tx_ids[i] == de_novo_transcript_ids[i] for i in range(len(tx_ids))]

        # Now we have to clean up the gene names.
        # If a stringtie gene overlaps a known gene at any point, we must use the known gene name.
        de_novo_gene_id_to_canonical_gene_name_dict = dict()
        for de_novo_gene_id, canonical_gene_id, canonical_gene_name in zip(de_novo_gene_ids, gene_ids, gene_names):
            if de_novo_gene_id != canonical_gene_id:
                try:
                    de_novo_gene_id_to_canonical_gene_name_dict[de_novo_gene_id].append((canonical_gene_id, canonical_gene_name))
                except KeyError:
                    de_novo_gene_id_to_canonical_gene_name_dict[de_novo_gene_id] = \
                        [(canonical_gene_id, canonical_gene_name)]

        # Make sure we don't have any ambiguous mappings:
        print("Ambiguously / multi-mapped:", file=sys.stderr)
        multimapped = 0
        ambiguous_de_novo_gene_set = set()
        for de_novo_gene_id, gene_info_tuple_list in de_novo_gene_id_to_canonical_gene_name_dict.items():
            if len(gene_info_tuple_list) > 1:
                conflicting = False
                for canonical_gene_id, _ in gene_info_tuple_list:
                    if canonical_gene_id != gene_info_tuple_list[0][0]:
                        conflicting = True
                        break
                if conflicting:
                    multimapped += 1
                    print(f"    {de_novo_gene_id}")
                    gene_info_tuple_count_dict = dict()
                    for canonical_gene_id, canonical_gene_name in gene_info_tuple_list:
                        try:
                            gene_info_tuple_count_dict[canonical_gene_id] += 1
                        except KeyError:
                            gene_info_tuple_count_dict[canonical_gene_id] = 1
                    for canonical_gene_id, canonical_gene_name in set(gene_info_tuple_list):
                        print(f"        {canonical_gene_id} -> {canonical_gene_name} ({gene_info_tuple_count_dict[canonical_gene_id]})")
                    ambiguous_de_novo_gene_set.add(de_novo_gene_id)
        print(f"Num ambiguously / multi-mapped = {multimapped}", file=sys.stderr)

        is_gene_id_ambiguous = [True if dngid in ambiguous_de_novo_gene_set else False for dngid in de_novo_gene_ids]

        # Now let's relabel the genes we know are not ambiguous:
        print("Relabeling unambiguous gene IDs and names...", end="", file=sys.stderr)
        num_relabeled = 0
        for i in range(len(is_de_novo)):
            if is_de_novo[i] and not is_gene_id_ambiguous[i]:
                de_novo_gene_id = de_novo_gene_ids[i]

                fixed_gene_id = None
                fixed_gene_name = None
                try:
                    fixed_gene_id = de_novo_gene_id_to_canonical_gene_name_dict[de_novo_gene_id][0][0]
                    fixed_gene_name = de_novo_gene_id_to_canonical_gene_name_dict[de_novo_gene_id][0][1]
                except KeyError:
                    pass

                if fixed_gene_id and fixed_gene_name:
                    gene_ids[i] = fixed_gene_id
                    gene_names[i] = fixed_gene_name
                    num_relabeled += 1
        print("Done!", file=sys.stderr)
        print(f"Successfully relabeled {num_relabeled} unambiguous gene names / IDs", file=sys.stderr)

        # Let downstream processing know we're not gencode:
        is_gencode = False

    # Create our anndata object now:
    count_adata = anndata.AnnData(count_mat.tocsr())

    # Add our variables:
    col_df = pd.DataFrame()
    col_df["transcript_ids"] = tx_ids
    col_df["gene_ids"] = gene_ids
    col_df["gene_names"] = gene_names
    col_df["transcript_names"] = transcript_names

    col_df["de_novo_gene_ids"] = de_novo_gene_ids
    col_df["de_novo_transcript_ids"] = de_novo_transcript_ids
    col_df["is_de_novo"] = is_de_novo
    col_df["is_gene_id_ambiguous"] = is_gene_id_ambiguous

    # If we're doing interval overlaps, add our label:
    if overlap_intervals:
        print(f"Adding {overlap_intervals_label} column for overlapping intervals...", file=sys.stderr)
        col_df[f"{overlap_intervals_label}"] = tx_overlap_flags
        print("Done!", file=sys.stderr)

    count_adata.var = col_df

    if is_gencode:
        count_adata.var_names = transcript_names
    else:
        count_adata.var_names = tx_ids

    # Add our observations:
    row_df = pd.DataFrame()
    row_df["Cell Barcode"] = cell_barcodes

    count_adata.obs = row_df
    count_adata.obs_names = cell_barcodes

    return count_adata


def read_intervals_from_tsv(filename):
    with open(filename, "r") as f, tqdm(desc="Processing intervals file", unit=" line") as pbar:
        intervals = []
        tsv_file = csv.reader(f, delimiter="\t")
        for row in tsv_file:
            if (not row[0].startswith("#")) and (not row[0].startswith("@")):
                try:
                    contig = row[0]
                    start = row[1]
                    end = row[2]
                    intervals.append((contig, start, end))
                except IndexError:
                    print("ERROR: Input interval file does not seem to be a properly formatted TSV.  "
                          "It must have the format: CONTIG    START    END", file=sys.stderr)
                    sys.exit(1)
            pbar.update(1)
    return intervals


def main(input_tsv, gtf_file, out_prefix, overlap_interval_filename=None, overlap_intervals_label=None):

    print("Verifying input file(s) exist...", file=sys.stderr)
    files_ok = True
    files_to_check = [input_tsv, gtf_file]
    if overlap_interval_filename:
        files_to_check.append(overlap_interval_filename)
    for f in files_to_check:
        if not os.path.exists(f):
            print(f"ERROR: Input file does not exist: {f}", file=sys.stderr)
            files_ok = False
    if not files_ok:
        sys.exit(1)

    overlap_intervals = None
    if overlap_interval_filename:
        print("Verifying contents of {overlap_interval_filename}...")
        overlap_intervals = read_intervals_from_tsv(overlap_interval_filename)

    print("Input files verified.", file=sys.stderr)

    # Create our gtf field map:
    gtf_field_dict = get_gtf_field_val_dict(gtf_file)

    # Create our anndata objects from the given data:
    print("Creating master anndata objects from transcripts counts data...", file=sys.stderr)
    master_adata = create_combined_anndata(input_tsv, gtf_field_dict, overlap_intervals, overlap_intervals_label)

    # Write our data out as pickles:
    print("Pickling data...", file=sys.stderr)
    pickle.dump(master_adata, open(f"{out_prefix}_tx_gene_counts_adata.pickle", "wb"))

    # Write our data as h5ad files:
    print("Writing data to h5ad file...", file=sys.stderr)
    master_adata.write(f"{out_prefix}_tx_gene_counts_adata.h5ad")

    print("Done!", file=sys.stderr)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=f"Creates anndata objects from the given count TSV and a GTF file.",
        epilog="The input TSV should have been created by create_count_matrix_from_annotated_bam.py and therefore"
               "should be of the form:"
               "Gene/TX    CBC    UMI    Count"
    )

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-t', '--tsv',
                               help='TSV file containing gene/transcript counts.',
                               required=True)
    requiredNamed.add_argument('-g', '--gtf',
                               help='GTF file containing gene annotations.  Can be from either Gencode or stringtie.',
                               required=True)
    requiredNamed.add_argument('-o', '--out-base-name',
                               help='Base name for the output files',
                               required=True)

    parser.add_argument("--overlap-intervals", 
                        help="TSV/Interval list file containing intervals on which to mark transcripts as overlapping.",
                        type=str)
    parser.add_argument("--overlap-interval-label", 
                        help="Label to add to all overlapping intervals.", 
                        type=str, 
                        default="")

    args = parser.parse_args()
    main(args.tsv, args.gtf, args.out_base_name, args.overlap_intervals, args.overlap_interval_label)
