version 1.0

import "Structs.wdl"

task IndexUnalignedBam {
    input {
        String input_bam

        RuntimeAttr? runtime_attr_override
    }

    String bri_path = basename(input_bam) + ".bri"

    command <<<
        set -euxo pipefail

        mkfifo /tmp/token_fifo
        ( while true ; do curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/token_fifo ; done ) &
        HTS_AUTH_LOCATION=/tmp/token_fifo bri index -v -i ~{bri_path} ~{input_bam}
    >>>

    output {
        File bri = bri_path
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            4,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-bri:0.1.19"
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

task MakeReadNameManifests {
    input {
        File input_bri
        Int N = 300

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        bri show ~{input_bri} | awk -F'/' 'BEGIN { getline; zmw=$2; line=$0 } { if ($2 != zmw) { print line; line = $0; } else { line = line "," $0; } zmw=$2; } END { print line; }' | sort -t'/' -n -k2 > read_list.txt
        split -a 5 -d --additional-suffix=".txt" -n l/~{N} read_list.txt chunk_
    >>>

    output {
        File manifest_full = "read_list.txt"
        Array[File] manifest_chunks = glob("chunk_*.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-bri:0.1.19"
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


task ExtractReadsInManifest {
    input {
        String input_bam
        File input_bri
        File read_name_manifest

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        sed 's/,/\n/g' ~{read_name_manifest} | sort -t'/' -n -k2 > readnames.txt

        mkfifo /tmp/token_fifo
        ( while true ; do curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/token_fifo ; done ) &
        ((HTS_AUTH_LOCATION=/tmp/token_fifo samtools view -H ~{input_bam}) && \
         (HTS_AUTH_LOCATION=/tmp/token_fifo cat readnames.txt | bri get -i ~{input_bri} ~{input_bam})) | samtools view -b > reads.bam

        wc -l readnames.txt
        samtools view reads.bam | wc -l
    >>>

    output {
        File reads = "reads.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-bri:0.1.19"
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
