version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/AlignReads.wdl" as AR
import "tasks/Finalize.wdl" as FF

workflow PBCCSOnlySingleFlowcell {
    input {
        String gcs_input_dir

        String? sample_name
        Int num_reads_per_split = 100000

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams { input: gcs_input_dir = gcs_input_dir }

    scatter (subread_bam in FindBams.subread_bams) {
        call PB.GetRunInfo { input: subread_bam = subread_bam }

        String SM  = select_first([sample_name, GetRunInfo.run_info["SM"]])
        String PU  = GetRunInfo.run_info["PU"]
        String ID  = PU
        String DIR = SM + "." + ID

        call Utils.ShardLongReads { input: unmapped_files = [ subread_bam ], num_reads_per_split = num_reads_per_split }

        scatter (subreads in ShardLongReads.unmapped_shards) {
            call PB.CCS { input: subreads = subreads }
        }

        call AR.MergeBams as MergeChunks { input: bams = CCS.consensus }
        call PB.MergeCCSReports as MergeCCSReports { input: reports = CCS.report }
    }

    call AR.MergeBams as MergeRuns { input: bams = MergeChunks.merged_bam, prefix = "~{SM[0]}.~{ID[0]}" }
    call PB.MergeCCSReports as MergeAllCCSReports { input: reports = MergeCCSReports.report }

    ##########
    # Finalize
    ##########

    call FF.FinalizeToDir as FinalizeMergedRuns {
        input:
            files = [ MergeRuns.merged_bam ],
            outdir = outdir + "/" + DIR[0] + "/alignments"
    }

    call FF.FinalizeToDir as FinalizeCCSMetrics {
        input:
            files = [ MergeAllCCSReports.report ],
            outdir = outdir + "/" + DIR[0] + "/metrics/ccs_metrics"
    }
}
