version 1.0

import "../Structs.wdl"

task Run_Group {

    meta {
        description : "Run umi-tools group on a bam file."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File aligned_transcriptome_reads
        File aligned_transcriptome_reads_index

        String gene_tag = "XG"
        String cell_barcode_tag = "CB"
        String umi_tag = "ZU"

        Boolean do_per_cell = true

        String prefix = "umi_tools_group"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 4*ceil(size(aligned_transcriptome_reads, "GB"))

    String per_cell_args = if do_per_cell then " --per-cell --cell-tag " + cell_barcode_tag + " " else ""

    command <<<
        set -euxo pipefail

        # Run umi-tools group:
        umi_tools group \
          --buffer-whole-contig \
          --no-sort-output \
          --per-gene \
          ~{per_cell_args} \
          --gene-tag XG \
          --extract-umi-method tag \
          --umi-tag ZU \
          -I ~{aligned_transcriptome_reads} \
          --group-out=~{prefix}.tsv \
          --output-bam \
          --log=~{prefix}.log > ~{prefix}.bam
    >>>

    output {
        File output_bam = "~{prefix}.bam"
        File output_tsv = "~{prefix}.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.3"
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
    }
}