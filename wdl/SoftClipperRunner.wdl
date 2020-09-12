version 1.0

import "tasks/SoftClipper.wdl" as SoftClipper
import "tasks/Finalize.wdl" as FF

workflow SoftClipperRunner {
    input {
        String reads_fastq
        File reference_fasta
        Int rounds
        Int clipping_threshold

        File output_dir
    }

    call SoftClipper.SplitSoftClippedReads {
        input:
            reads_fastq = reads_fastq,
            reference_fasta = reference_fasta,
            rounds = rounds,
            clipping_threshold = clipping_threshold
    }

    call FF.FinalizeToDir {
        input:
            files = SplitSoftClippedReads.split_reads,
            outdir = output_dir
    }
}