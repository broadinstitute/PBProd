version 1.0

import "tasks/Canu.wdl" as Canu
import "tasks/Quast.wdl" as Quast
import "tasks/Finalize.wdl" as FF

workflow MalariaAssembly {
    # needs evaluate

    input {
        String gcs_fast5_dir
        String medaka_model

        String out_dir
    }

    call Guppy.Guppy {
        input:
            gcs_fast5_dir = gcs_fast5_dir
    }

    call Guppy.MergeFastq {
        input:
            guppy_output_files = Guppy.output_files
    }

    call Canu.CorrectTrimAssemble {
        input:
            output_file_prefix = "malaria",
            genome_size = "22.9m",
            reads = MergeFastq.merged_fastq,
            correct_corrected_error_rate = 0.15,
            trim_corrected_error_rate = 0.15,
            assemble_corrected_error_rate = 0.15
    }

    call Medaka.MedakaPolish {
        input:
            basecalled_reads = MergeFastq.merged_fastq,
            draft_assembly = CorrectTrimAssemble.canu_contigs_fasta,
            model = medaka_model,
            n_rounds = 3
    }

    call FF.FinalizeToDir {
        input:
            files = [MedakaPolish.polished_assembly],
            outdir = out_dir
    }
}
