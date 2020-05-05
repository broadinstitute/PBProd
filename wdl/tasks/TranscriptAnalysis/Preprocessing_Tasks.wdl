version 1.0

task SplitBamBySampleAndCellBarcodeTask {

    meta {
        description : "Convert a single annotated (via the 10x tool), aligned bam file into individual FASTA files named by sample name and cell barcode.  Also produces a manifest file for FLAIR to easily quantify output."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File aligned_annotated_bam
        String output_base_name = "reads"
    }

    parameter_meta {
        aligned_annotated_bam : "Bam file containing aligned reads that have been annotated with the 10x tool."
        output_base_name : "[optional] base name to give to every output file.  Should correspond to some unique identifier from this dataset."
    }

    # 10x the total size of the input bam (uncompressed reads)
    # 1x for the file itself
    # 1x for wiggle-room
    # 2x for tar/gz-ing the output:
    Int disk_size = ((10+1+1)*2) * ceil(size(aligned_annotated_bam, "GB"))

    String fasta_tar_gz_name = "fasta_files_by_sample_and_barcode.tar.gz"

    command {
        /python_scripts/split_annotated_reads_by_sample_and_cell_barcode.py -b ~{aligned_annotated_bam} -o ~{output_base_name}
        tar -zcf ~{fasta_tar_gz_name} *.fasta
    }

    output {
        File flair_manifest = "${output_base_name}_flair_reads_manifest.tsv"
        Array[File] sample_cell_barcode_fasta_files = glob("*.fasta")
        File fasta_tar_gz_out = "~{fasta_tar_gz_name}"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.1"
        memory: 16 + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        boot_disk_gb: 10
        preemptible: 0
        cpu: 8
    }
}
