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

GENE_NAME_FIELD = "gene_name"
GENE_ID_FIELD = "gene_id"
TX_NAME_FIELD = "transcript_name"
TX_ID_FIELD = "transcript_id"

ALT_NAME_SUFFIX = "_PAR_Y"


def get_gtf_field_val_dict(gencode_gtf_file, force_rebuild=False):

    global TX_ENTRY_STRING
    global GENE_NAME_FIELD
    global GENE_ID_FIELD
    global TX_NAME_FIELD
    global TX_ID_FIELD
    global ALT_NAME_SUFFIX

    gtf_dict = dict()

    if not os.path.exists(gencode_gtf_file):
        raise FileNotFoundError(f"Gencode GTF file does not exist: {gencode_gtf_file}")

    pickle_file_name = os.path.splitext(os.path.basename(gencode_gtf_file))[0] + ".big_map.pickle"

    if not force_rebuild and os.path.exists(pickle_file_name):
        print(f"Loading gencode map from {pickle_file_name}...", end="\t", file=sys.stderr)
        gtf_dict = pickle.load(open(pickle_file_name, "rb"))
        print("Done!", file=sys.stderr)
    else:
        print(f"Creating master GTF field map keyed by transcript ID using {gencode_gtf_file}...",
              flush=True, file=sys.stderr)

        with open(gencode_gtf_file, "r") as f, tqdm(desc="Processing Gencode File", unit="line") as pbar:
            tsv_file = csv.reader(f, delimiter="\t")
            for row in tsv_file:
                # Ignore comments and make sure we only process transcript entries:
                if not row[0].startswith("#") and row[2] == TX_ENTRY_STRING:

                    # Parse the row data into a dict and make sure we're only looking at protein coding rows:
                    row_data_dict = {
                        field.strip().split(" ")[0].replace('"', ""): field.strip().split(" ")[1].replace('"', "") for
                        field in row[8].split(";") if len(field) != 0}

                    # Make sure our names are unique:
                    if row_data_dict[TX_ID_FIELD].endswith(ALT_NAME_SUFFIX):
                        row_data_dict[TX_NAME_FIELD] = row_data_dict[TX_NAME_FIELD] + ALT_NAME_SUFFIX
                        row_data_dict[GENE_NAME_FIELD] = row_data_dict[GENE_NAME_FIELD] + ALT_NAME_SUFFIX

                    # Add this row to our dict keyed by transcript ID:
                    gtf_dict[row_data_dict[TX_ID_FIELD]] = row_data_dict

                pbar.update(1)

        print("Pickling data...", file=sys.stderr)
        pickle.dump(gtf_dict, open(pickle_file_name, "wb"))
        print("Done!", file=sys.stderr)

    return gtf_dict


def create_combined_anndata(input_tsv, gtf_field_dict, force_recount=False):
    """Create an anndata object holding the given gene/transcript information.
    NOTE: This MUST be a sparse matrix - we got lots of data here.

    First convert the counts to a matrix that looks like what scanpy expects:
    columns = Transcripts name / Transcript ID / Gene Name / Gene ID (variables)
    Rows = cell barcodes (observations)
    data = counts"""

    global GENE_NAME_FIELD
    global GENE_ID_FIELD
    global TX_NAME_FIELD
    global TX_ID_FIELD

    cell_barcode_to_tx_to_umi_dict = dict()
    cell_barcode_to_tx_count_dict = dict()

    pickle_file_name = os.path.splitext(os.path.basename(input_tsv))[0] + ".tx_raw_count_matrix.pickle"

    if not force_recount and os.path.exists(pickle_file_name):
        print(f"Loading count map from {pickle_file_name}...", end="\t", file=sys.stderr)
        cell_barcode_to_tx_count_dict = pickle.load(open(pickle_file_name, "rb"))
        print("Done!", file=sys.stderr)
    else:
        # Get our cell tx counts:
        with open(input_tsv, "r") as f, tqdm(desc="Processing Raw Cell Counts", unit="count") as pbar:
            tsv_file = csv.reader(f, delimiter="\t")
            next(tsv_file)
            for row in tsv_file:
                tx_col = row[0]
                tx_id = tx_col[:tx_col.find("|")]

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
    count_mat = scipy.sparse.lil_matrix((len(cell_barcodes), len(tx_ids)), dtype=np.uint32)

    tx_id_index_dict = {name: i for i, name in enumerate(tx_ids)}
    with tqdm(desc=f"Creating cell transcript count matrix", unit="cell",
              total=len(cell_barcode_to_tx_count_dict)) as pbar:

        for i, (cb, counts_dict) in enumerate(cell_barcode_to_tx_count_dict.items()):
            # Put the counts for each transcript in the right indices:
            for tx_id, count in counts_dict.items():
                count_mat[i, tx_id_index_dict[tx_id]] = count
            pbar.update(1)

    # Now we set up the variables that we're going to apply to each observation:
    transcript_names = [gtf_field_dict[tx_id][TX_NAME_FIELD] for tx_id in tx_ids]
    gene_ids = [gtf_field_dict[tx_id][GENE_ID_FIELD] for tx_id in tx_ids]
    gene_names = [gtf_field_dict[tx_id][GENE_NAME_FIELD] for tx_id in tx_ids]

    # Create our anndata object now:
    count_adata = anndata.AnnData(count_mat.tocsr())

    # Add our variables:
    col_df = pd.DataFrame()
    col_df["transcript_names"] = transcript_names
    col_df["transcript_ids"] = tx_ids
    col_df["gene_names"] = gene_names
    col_df["gene_ids"] = gene_ids

    count_adata.var = col_df
    count_adata.var_names = transcript_names

    # Add our observations:
    row_df = pd.DataFrame()
    row_df["Cell Barcode"] = cell_barcodes

    count_adata.obs = row_df
    count_adata.obs_names = cell_barcodes

    return count_adata


def main(input_tsv, gencode_gtf, out_prefix):

    print("Verifying input file(s) exist...", file=sys.stderr)
    files_ok = True
    for f in [input_tsv, gencode_gtf]:
        if not os.path.exists(f):
            print(f"ERROR: Input file does not exist: {f}", file=sys.stderr)
            files_ok = False
    if not files_ok:
        sys.exit(1)
    print("Input files verified.", file=sys.stderr)

    # Create our gencode map:
    gtf_field_dict = get_gtf_field_val_dict(gencode_gtf)

    # Create our anndata objects from the given data:
    print("Creating master anndata objects from transcripts counts data...", file=sys.stderr)
    master_adata = create_combined_anndata(input_tsv, gtf_field_dict)

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
                               help='Gencode GTF file containing gene annotations.',
                               required=True)
    requiredNamed.add_argument('-o', '--out-base-name',
                               help='Base name for the output files',
                               required=True)

    args = parser.parse_args()
    main(args.tsv, args.gtf, args.out_base_name)
