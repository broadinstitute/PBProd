version 1.0

import "Structs.wdl"

task Annotate
{
    input {
        File reads
        Boolean is_mas_seq_10_array = false
        String prefix = "annmas_annotated"

        RuntimeAttr? runtime_attr_override
    }

    String model_spec_arg = if is_mas_seq_10_array then " --m10 " else ""

    Int disk_size = 4*ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        source /annmas/venv/bin/activate
        annmas annotate -t8 -v INFO ~{reads} ~{model_spec_arg} -o ~{prefix}.bam
    >>>

    output {
        File annotated_bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-annmas:0.0.4"
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

task Segment
{
    input {
        File annotated_reads
        Boolean is_mas_seq_10_array = false
        String prefix = "annmas_segmented"

        RuntimeAttr? runtime_attr_override
    }

    String model_spec_arg = if is_mas_seq_10_array then " --m10 " else ""
    Int disk_size = 4*ceil(size(annotated_reads, "GB"))

    command <<<
        set -euxo pipefail

        source /annmas/venv/bin/activate
        annmas segment -v INFO -s ~{annotated_reads} ~{model_spec_arg} -o ~{prefix}.bam
    >>>

    output {
        File segmented_bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-annmas:0.0.4"
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

task ScSplit
{
    input {
        File reads_bam
        Boolean is_mas_seq_10_array = false

        String force_option = ""

        Int umi_length = 10
        String cbc_dummy = "CTGCCTAACCTGATCC"

        String prefix = "scsplit"

        RuntimeAttr? runtime_attr_override
    }

    String model_spec_arg = if is_mas_seq_10_array then " --m10 " else ""

    # On average bam compression is 10x and we have more data on output than that:
    Int disk_size = 15*ceil(size(reads_bam, "GB"))

    command <<<
        set -euxo pipefail

        source /annmas/venv/bin/activate
        annmas scsplit -t4 -v INFO ~{model_spec_arg} -b ~{force_option} -o ~{prefix} -u ~{umi_length} -c ~{cbc_dummy} ~{reads_bam}
    >>>

    output {
        File mates1 = "~{prefix}_mates1.fastq"
        File mates2 = "~{prefix}_mates2.fastq"
        File annotated_bam = "~{prefix}.cbc_umi_annotated.bam"
        File whitelist = "~{prefix}_whitelist.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-annmas:0.0.4"
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

task Inspect
{
    input {
        File reads
        File reads_pbi

        File read_names

        Boolean is_mas_seq_10_array = false
        String prefix = "annmas_inspected_reads"

        RuntimeAttr? runtime_attr_override
    }

    String model_spec_arg = if is_mas_seq_10_array then " --m10 " else ""
    Int disk_size = 4*ceil(size(reads, "GB")) + size(reads_pbi, "GB") + size(read_names, "GB")

    command <<<
        set -euxo pipefail

        source /annmas/venv/bin/activate
        annmas inspect ~{model_spec_arg} ~{reads} -p ~{reads_pbi} -r ~{read_names} -o ~{prefix}
        tar -zxf ~{prefix}.tar.gz ~{prefix}
    >>>

    output {
        File inspected_reads_tgz = "~{prefix}.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-annmas:0.0.4"
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
