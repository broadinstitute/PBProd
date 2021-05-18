version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF
import "tasks/Longbow.wdl" as LONGBOW

workflow PB10xMasSeqArrayPreProcessing {

    meta {
        description : "This workflow is designed to pre-process data from the MASSeq v2 protocol.  It will first filter PacBio Sequel IIe reads to keep only the high quality reads.  Then it will create several bam files - one for each model - and populate those files with the reads that best fit those models."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        String gcs_input_dir
        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/PB10xMasSeqArrayPreProcessing"

        # Default here is 0 because ccs uncorrected reads all seem to have RQ = -1.
        # All pathologically long reads also have RQ = -1.
        # This way we preserve the vast majority of the data, even if it has low quality.
        # We can filter it out at later steps.
        Float min_read_quality = 0.0

        String? sample_name
    }

    parameter_meta {
        gcs_input_dir : "Input folder on GCS in which to search for BAM files to process."
        gcs_out_root_dir : "Root output GCS folder in which to place results of this workflow."

        min_read_quality : "[optional] Min quality necessary to include a read.  Sequel IIe CCS uncorrected reads all seem to have RQ = -1.  All pathologically long reads also seem to have RQ = -1 (this makes sense given that the pathologically long reads I've seen are all reads that failed CCS).  (default: 0.0)."

        sample_name : "[optional] The name of the sample to associate with the data in this workflow."
    }

    # Create an object to disable preemption.  This should only be used for testing.
    RuntimeAttr disable_preemption = object {
        preemptible_tries:  0
    }

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as WdlExecutionStartTimestamp { input: }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams { input: gcs_input_dir = gcs_input_dir }

    scatter (reads_bam in FindBams.ccs_bams) {

        call PB.GetRunInfo { input: subread_bam = reads_bam }

        String SM  = select_first([sample_name, GetRunInfo.run_info["SM"]])
        String PL  = "PACBIO"
        String PU  = GetRunInfo.run_info["PU"]
        String DT  = GetRunInfo.run_info["DT"]
        String ID  = PU
        String DS  = GetRunInfo.run_info["DS"]
        String DIR = SM + "." + ID

        String RG_subreads  = "@RG\\tID:~{ID}.subreads\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
        String RG_consensus = "@RG\\tID:~{ID}.consensus\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
        String RG_array_elements = "@RG\\tID:~{ID}.array_elements\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

        # Shard widely so we can go faster:
        File read_pbi = sub(reads_bam, ".bam$", ".bam.pbi")
        call PB.ShardLongReads {
            input:
                unaligned_bam = reads_bam,
                unaligned_pbi = read_pbi,
                prefix = SM + "_shard",
                runtime_attr_override = disable_preemption
        }

        scatter (sharded_reads in ShardLongReads.unmapped_shards) {

            # 1 - filter the reads by the minimum read quality:
            String fbmrq_prefix = basename(sharded_reads, ".bam")
            call Utils.Bamtools as FilterByMinReadQuality {
                input:
                    bamfile = sharded_reads,
                    prefix = fbmrq_prefix + "_good_reads",
                    cmd = "filter",
                    args = '-tag "rq":">=' + min_read_quality + '"',
                    runtime_attr_override = disable_preemption
            }

            # 2 - split the reads by the model:
            String adis_prefix = basename(FilterByMinReadQuality.bam_out, ".bam")
            call LONGBOW.Demultiplex as AssignReadsToModels {
                input:
                    bam = FilterByMinReadQuality.bam_out,
                    prefix = adis_prefix,
                    runtime_attr_override = disable_preemption
            }
        }
        # TODO: Fix this to allow for an arbitrary number of models easily:
        call Utils.MergeBams as MergeMas10Bams {
            input:
                bams = AssignReadsToModels.mas10_bam,
                prefix = SM + "_mas10.reads",
                runtime_attr_override = disable_preemption
        }
        call PB.PBIndex as PbIndexMas10Bam {
            input:
                bam = MergeMas10Bams.merged_bam,
                runtime_attr_override = disable_preemption
        }

        call Utils.MergeBams as MergeMas15Bams {
            input:
                bams = AssignReadsToModels.mas15_bam,
                prefix = SM + "_mas15.reads",
                runtime_attr_override = disable_preemption
        }
        call PB.PBIndex as PbIndexMas15Bam {
            input:
                bam = MergeMas15Bams.merged_bam,
                runtime_attr_override = disable_preemption
        }

        # Finalize our merged reads:
        String base_out_dir = outdir + "/" + DIR + "/" + WdlExecutionStartTimestamp.timestamp_string
        call FF.FinalizeToDir as FinalizeMas10SplitBam {
            input:
                files = [
                    MergeMas10Bams.merged_bam,
                    MergeMas10Bams.merged_bai,
                    PbIndexMas10Bam.pbindex,
                ],
                outdir = base_out_dir + "/" + SM + "_mas10",
                runtime_attr_override = disable_preemption
        }
        call FF.FinalizeToDir as FinalizeMas15SplitBam {
            input:
                files = [
                    MergeMas15Bams.merged_bam,
                    MergeMas15Bams.merged_bai,
                    PbIndexMas15Bam.pbindex,
                ],
                outdir = base_out_dir + "/" + SM + "_mas15",
                runtime_attr_override = disable_preemption
        }

        # Copy over the metadata to our finalized folders so we can just run our workflow on it:
        call PB.CopyMetadataFilesToNewDir as CopyMetadataForMas10Reads {
            input:
                input_gs_path = gcs_input_dir,
                dest_gs_path = base_out_dir + "/" + SM + "_mas10"
        }
        call PB.CopyMetadataFilesToNewDir as CopyMetadataForMas15Reads {
            input:
                input_gs_path = gcs_input_dir,
                dest_gs_path = base_out_dir + "/" + SM + "_mas15"
        }
    }
}