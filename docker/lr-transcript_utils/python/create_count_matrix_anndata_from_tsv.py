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


def create_transcript_and_gene_name_maps(gencode_gtf_file="gencode.v34.primary_assembly.annotation.gtf",
                                         force_rebuild=False):
    """Create transcript and gene name maps from the given gencode GTF file."""

    if not os.path.exists(gencode_gtf_file):
        raise FileNotFoundError(f"Gencode GTF file does not exist: {gencode_gtf_file}")

    tx_id_name_map_pickle_file_name = os.path.splitext(os.path.basename(gencode_gtf_file))[0] + ".tx_id_name_map.pickle"
    tx_id_gene_map_pickle_file_name = os.path.splitext(os.path.basename(gencode_gtf_file))[0] + ".tx_id_gene_map.pickle"

    if not force_rebuild and os.path.exists(tx_id_name_map_pickle_file_name) and os.path.exists(
            tx_id_gene_map_pickle_file_name):
        print(f"Loading tx id name map from {tx_id_name_map_pickle_file_name}...", end="\t", file=sys.stderr)
        tx_id_name_map = pickle.load(open(tx_id_name_map_pickle_file_name, "rb"))
        print("Done!", file=sys.stderr)

        print(f"Loading tx id gene map from {tx_id_gene_map_pickle_file_name}...", end="\t", file=sys.stderr)
        tx_id_gene_map = pickle.load(open(tx_id_gene_map_pickle_file_name, "rb"))
        print("Done!", file=sys.stderr)
    else:
        print(f"Creating transcript and gene name maps from {gencode_gtf_file}...", flush=True, file=sys.stderr)
        tx_id_name_map = dict()
        tx_id_gene_map = dict()

        TX_ENTRY_STRING = "transcript"

        GENE_NAME_FIELD = "gene_name"
        TX_NAME_FIELD = "transcript_name"
        TX_ID_FIELD = "transcript_id"
        GENE_TYPE_FIELD = "gene_type"
        PROTEIN_CODING_VALUE = "protein_coding"

        ALT_NAME_SUFFIX = "_PAR_Y"

        with open(gencode_gtf_file, "r") as f, tqdm(desc="Processing Gencode File", unit="line") as pbar:
            tsv_file = csv.reader(f, delimiter="\t")
            for row in tsv_file:
                # Ignore comments and make sure we only process transcript entries:
                if not row[0].startswith("#") and row[2] == TX_ENTRY_STRING:

                    # Parse the row data into a dict and make sure we're only looking at protein coding rows:
                    row_data_dict = {
                        field.strip().split(" ")[0].replace('"', ""): field.strip().split(" ")[1].replace('"', "") for
                        field in row[8].split(";") if len(field) != 0}

                    # Disable the protein coding check for now:
                    #                     if row_data_dict[GENE_TYPE_FIELD] != PROTEIN_CODING_VALUE:
                    #                         pbar.update(1)
                    #                         continue

                    # Add our data to our maps:
                    if row_data_dict[TX_ID_FIELD].endswith(ALT_NAME_SUFFIX):
                        tx_id_name_map[row_data_dict[TX_ID_FIELD]] = row_data_dict[TX_NAME_FIELD] + ALT_NAME_SUFFIX
                        tx_id_gene_map[row_data_dict[TX_ID_FIELD]] = row_data_dict[GENE_NAME_FIELD] + ALT_NAME_SUFFIX
                    else:
                        tx_id_name_map[row_data_dict[TX_ID_FIELD]] = row_data_dict[TX_NAME_FIELD]
                        tx_id_gene_map[row_data_dict[TX_ID_FIELD]] = row_data_dict[GENE_NAME_FIELD]

                pbar.update(1)

        print("Pickling data...", file=sys.stderr)
        pickle.dump(tx_id_name_map, open(tx_id_name_map_pickle_file_name, "wb"))
        pickle.dump(tx_id_gene_map, open(tx_id_gene_map_pickle_file_name, "wb"))
        print("Done!", file=sys.stderr)

    return tx_id_name_map, tx_id_gene_map


def get_cell_transcript_counts(filename, tx_id_name_map=None, tx_id_gene_map=None):
    """Collapse the given counts into two nested dictionaries:
    {Cell barcode : {transcript : count}}
    {Cell barcode : {gene : count}}
    """

    print(f"Extracting TX and Gene counts for {filename}", flush=True, file=sys.stderr)

    cell_transcript_umi_store = dict()

    cell_transcript_counts = dict()
    cell_gene_counts = dict()

    with open(filename, "r") as f, tqdm(desc="Processing Raw Cell Counts", unit="count") as pbar:
        tsv_file = csv.reader(f, delimiter="\t")
        next(tsv_file)
        for row in tsv_file:

            gene_name = None

            tx_id = row[0]
            tx_name = tx_id

            if tx_id_name_map:
                # If we have a name map, we should rename the transcript to its real name:
                tx_name = tx_id_name_map[tx_id[:tx_id.find("|")]]

            if tx_id_gene_map:
                # If we have a name map, we should rename the transcript to its real name:
                gene_name = tx_id_gene_map[tx_id[:tx_id.find("|")]]

            # Handle the Transcript names:
            if row[1] not in cell_transcript_umi_store:
                cell_transcript_umi_store[row[1]] = dict()
                cell_transcript_umi_store[row[1]][tx_name] = set()
                cell_transcript_umi_store[row[1]][tx_name].add(row[2])

                cell_transcript_counts[row[1]] = dict()
                cell_transcript_counts[row[1]][tx_name] = 1

                if gene_name:
                    cell_gene_counts[row[1]] = dict()
                    cell_gene_counts[row[1]][gene_name] = 1

            elif tx_name not in cell_transcript_umi_store[row[1]]:
                cell_transcript_umi_store[row[1]][tx_name] = set()
                cell_transcript_umi_store[row[1]][tx_name].add(row[2])

                cell_transcript_counts[row[1]][tx_name] = 1

                if gene_name:
                    cell_gene_counts[row[1]][gene_name] = 1

            elif row[2] not in cell_transcript_umi_store[row[1]][tx_name]:
                cell_transcript_umi_store[row[1]][tx_name].add(row[2])

                cell_transcript_counts[row[1]][tx_name] += 1

                if gene_name:
                    cell_gene_counts[row[1]][gene_name] += 1

            pbar.update(1)

    print("Done!", file=sys.stderr)
    return cell_transcript_counts, cell_gene_counts


def create_anndata(counts_nested_dict, name_map, col_field_name="Gene"):
    """Create an anndata object holding the given gene/transcript information.
    NOTE: This MUST be a sparse matrix - we got lots of data here.

    First convert the counts to a matrix that looks like what scanpy expects:
    columns = Transcripts/Genes (variables)
    Rows = cell barcodes (observations)
    data = counts"""
    # Do Genes first:
    row_labels = np.array([cb for cb, _ in counts_nested_dict.items()])
    col_labels = np.unique(np.array([name for _, name in name_map.items()]))
    count_mat = scipy.sparse.lil_matrix((len(row_labels), len(col_labels)), dtype=np.uint32)

    # Populate the count matrix:
    name_index_dict = {name: i for i, name in enumerate(col_labels)}
    with tqdm(desc=f"Creating cell {col_field_name.lower()} count matrix", unit="cell",
              total=len(counts_nested_dict)) as pbar:
        for i, (cb, counts_dict) in enumerate(counts_nested_dict.items()):
            # Put the counts for each gene in the right indices:
            for name, count in counts_dict.items():
                count_mat[i, name_index_dict[name]] = count
            pbar.update(1)

    count_adata = anndata.AnnData(count_mat.tocsr())
    col_df = pd.DataFrame()
    col_df[col_field_name] = col_labels

    row_df = pd.DataFrame()
    row_df["Cell Barcode"] = row_labels

    count_adata.var = col_df
    count_adata.var_names = col_labels

    count_adata.obs = row_df
    count_adata.obs_names = row_labels

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

    # Create our maps for TX and Gene names:
    print("Creating transcript / gene ID -> name maps...", file=sys.stderr)
    tx_id_name_map, tx_id_gene_map = create_transcript_and_gene_name_maps(gencode_gtf)

    # Get the counts themselves for the maps we created:
    print("Tallying counts from input data...", file=sys.stderr)
    tx_counts, gene_counts = get_cell_transcript_counts(input_tsv, tx_id_name_map, tx_id_gene_map)

    # Create our anndata objects from the given data:
    print("Creating anndata objects from counts data...", file=sys.stderr)
    gene_count_adata = create_anndata(gene_counts, tx_id_gene_map)
    tx_count_adata = create_anndata(tx_counts, tx_id_name_map, col_field_name="Transcript")

    # Write our data out as pickles:
    print("Pickling data...", file=sys.stderr)
    pickle.dump(gene_count_adata, open(f"{out_prefix}_gene_counts_adata.pickle", "wb"))
    pickle.dump(tx_count_adata, open(f"{out_prefix}_tx_counts_adata.pickle", "wb"))

    # Write our data as h5ad files:
    print("Writing data to h5ad files...", file=sys.stderr)
    gene_count_adata.write(f"{out_prefix}_gene_counts_adata.h5ad")
    tx_count_adata.write(f"{out_prefix}_tx_counts_adata.h5ad")

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
