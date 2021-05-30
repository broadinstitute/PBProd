version 1.0

#######################################################
# This pipeline calls small variants using DeepVariant.
#######################################################

import "Structs.wdl"

task Clair {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        File? sites_vcf
        File? sites_vcf_tbi

        String chr
        String preset

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:             "input BAM from which to call variants"
        bai:             "index accompanying the BAM"

        ref_fasta:       "reference to which the BAM was aligned to"
        ref_fasta_fai:   "index accompanying the reference"

        chr:             "chr on which to call variants"
        preset:          "calling preset (CCS, ONT)"
    }

    Int disk_size = 2*ceil(size(select_all([bam, bai, ref_fasta, ref_fasta_fai, sites_vcf]), "GB"))
    String platform = if preset == "CCS" then "hifi" else "ont"
    String vcf_fn = if defined(sites_vcf) then "--vcf_fn=" + select_first([sites_vcf]) else ""

    command <<<
        # example from https://github.com/HKU-BAL/Clair3#option-1--docker-pre-built-image
        set -euxo pipefail

        touch ~{bai}

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)
        SM=$(samtools view -H ~{bam} | grep -m1 '^@RG' | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g')

        /opt/bin/run_clair3.sh \
            --bam_fn=~{bam} --ref_fn=~{ref_fasta} ~{vcf_fn} \
            --threads=${num_core} \
            --platform=~{platform} \
            --model_path="/opt/models/~{platform}" \
            --ctg_name=~{chr} \
            --sample_name=$SM \
            --gvcf \
            --output="./"

        find . -type f -exec ls -lah {} \;
    >>>

    output {
#        # save both VCF and gVCF
#        File phased_vcf = "~{prefix}.phased.vcf.gz"
#        File phased_vcf_tbi = "~{prefix}.phased.vcf.gz.tbi"
#
#        File vcf = "~{prefix}.vcf.gz"
#        File vcf_tbi = "~{prefix}.vcf.gz.tbi"
#        File gvcf = "~{prefix}.g.vcf.gz"
#        File gvcf_tbi = "~{prefix}.g.vcf.gz.tbi"
#
#        File report = "~{prefix}.visual_report.html"
#        File phaseset_bed = "~{prefix}.phaseset.bed"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          12,
        mem_gb:             72,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "kishwars/pepper_deepvariant:r0.4.1"
        #docker:             "us.gcr.io/broad-dsp-lrma/lr-dvpepper:r0.4.1"
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
        gpuCount:               1
        nvidiaDriverVersion:    "418.152.00"
        zones:                  ["us-central1-c", "us-central1-f", "us-east1-b", "us-east1-c", "us-west1-a", "us-west1-b"]
        cpuPlatform:            "Intel Skylake"
    }
}
