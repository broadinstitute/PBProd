version 1.0

##########################################################################################
# This pipeline calls small variants using DeepVariant on PacBio CCS BAM.
##########################################################################################

import "Structs.wdl"

task DeepVariant {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fai

        String model_class

        String output_prefix

        String? intervals

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam: "input BAM from which to call variants"
        bai: "index accompanying the BAM"

        ref_fasta: "reference to which the BAM was aligned to"
        ref_fai:   "index accompanying the reference"

        model_class: "class of model to be applied; currently only 'PACBIO' is accepted"

        intervals: "[optional] interval on which to call, '<chr>:<start>-<end>'; if ommited, the whole genome"

        output_prefix: "prefix to output files"
    }

    Int disk_size = ceil(size(bam, "GB")) + 50

    String extra_args = if defined(intervals) then "--regions=~{intervals}" else " "

    command <<<

        # example from https://github.com/google/deepvariant/blob/r0.8/docs/deepvariant-quick-start.md
        set -euo pipefail

        # this is just to future proof ourselves for ONT data
        if [[ ~{model_class} != "PACBIO" ]]; then
          echo "requesting a model type not supported by the model: ~{model_class}"
          exit 1
        fi

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        /opt/deepvariant/bin/run_deepvariant \
          --model_type=~{model_class} \
          --ref=~{ref_fasta} \
          --reads=~{bam} \
          --output_vcf=/cromwell_root/~{output_prefix}.deepvariant.vcf.gz \
          --output_gvcf=/cromwell_root/~{output_prefix}.deepvariant.g.vcf.gz \
          --num_shards=${num_core} \
          ~{extra_args}
    >>>

    output {
        # save both VCF and gVCF
        File gvcf = "~{output_prefix}.deepvariant.g.vcf.gz"
        File gvcf_tbi = "~{output_prefix}.deepvariant.g.vcf.gz.tbi"
        File vcf = "~{output_prefix}.deepvariant.vcf.gz"
        File vcf_tbi = "~{output_prefix}.deepvariant.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             120,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "gcr.io/deepvariant-docker/deepvariant:0.8.0-gpu"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        gpuType:                "nvidia-tesla-p100"
        gpuCount:               2
        nvidiaDriverVersion:    "418.87.00"
        zones:                  ["us-east1-b", "us-east1-c"]
        cpuPlatform:            "Intel Skylake"
    }
}
