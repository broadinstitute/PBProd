version 1.0

import "Structs.wdl"

# A wrapper to minimap2 for mapping & aligning (groups of) sequences to a reference
task Minimap2 {
    input {
        Array[File] reads
        File ref_fasta

        String RG
        String map_preset

        String prefix = "out"
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        reads:      "query sequences to be mapped and aligned"
        ref_fasta:  "reference fasta"
        RG:         "read group information to be supplied to parameter '-R' (note that tabs should be input as '\t')"
        map_preset: "preset to be used for minimap2 parameter '-x'"
        prefix:     "[default-valued] prefix for output BAM"
    }

    Int disk_size = 1 + 3*ceil(size(reads, "GB") + size(ref_fasta, "GB"))

    Int cpus = 4
    Int mem = 30

    command <<<
        set -euxo pipefail

        MAP_PARAMS="-ayYL --MD -x ~{map_preset} -R ~{RG} -t ~{cpus} ~{ref_fasta}"
        SORT_PARAMS="-@~{cpus} -m~{mem}G --no-PG -o ~{prefix}.bam"
        FILE="~{reads[0]}"

        if [[ "$FILE" =~ \.fastq$ ]] || [[ "$FILE" =~ \.fq$ ]]; then
            find . \( -name '*.fastq' -or -name '*.fq' \) -exec cat {} \; | minimap2 $MAP_PARAMS - | samtools sort $SORT_PARAMS -
        elif [[ "$FILE" =~ \.fastq.gz$ ]] || [[ "$FILE" =~ \.fq.gz$ ]]; then
            find . \( -name '*.fastq.gz' -or -name '*.fq.gz' \) -exec zcat {} \; | minimap2 $MAP_PARAMS - | samtools sort $SORT_PARAMS -
        elif [[ "$FILE" =~ \.fasta$ ]] || [[ "$FILE" =~ \.fa$ ]]; then
            find . \( -name '*.fasta' -or -name '*.fa' \) -not -name '~{basename(ref_fasta)}' -exec cat {} \; | python3 /usr/local/bin/cat_as_fastq.py | minimap2 $MAP_PARAMS - | samtools sort $SORT_PARAMS -
        elif [[ "$FILE" =~ \.fasta.gz$ ]] || [[ "$FILE" =~ \.fa.gz$ ]]; then
            find . \( -name '*.fasta.gz' -or -name '*.fa.gz' \) -not -name '~{basename(ref_fasta)}' -exec zcat {} \; | python3 /usr/local/bin/cat_as_fastq.py | minimap2 $MAP_PARAMS - | samtools sort $SORT_PARAMS -
        elif [[ "$FILE" =~ \.bam$ ]]; then
            samtools fastq $FILES | minimap2 $MAP_PARAMS - | samtools sort $SORT_PARAMS -
        else
            echo "Did not understand file format for '$FILE'"
            exit 1
        fi

        samtools index ~{prefix}.bam
    >>>

    output {
        File aligned_bam = "~{prefix}.bam"
        File aligned_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             mem,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.27"
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

# A simple task to covert SAM-formatted alignment to PAF format
task SAMtoPAF {
    input {
        File sam_formatted_file
        File? index

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        sam_formatted_file: "SAM-formated input file to be converted to PAF (note currently we only support SAM or BAM, not CRAM)"
        index:              "[optional] index for sam_formatted_file"
    }

    String prefix = basename(basename(sam_formatted_file, ".bam"), ".sam") # we have hack like this because WDL stdlib doesn't provide endsWith stuff

    Int disk_size = 2*ceil(size(sam_formatted_file, "GB"))

    command <<<
        set -eu

        filename=$(basename -- ~{sam_formatted_file})
        extension="${filename##*.}"
        if [[ "$extension" == "sam" ]]; then
            /minimap2-2.17_x64-linux/k8 \
                /minimap2-2.17_x64-linux/paftools.js \
                sam2paf \
                -L \
                ~{sam_formatted_file} \
                > ~{prefix}".paf"
        elif [[ "$extension" == "bam" ]]; then
            samtools view -h ~{sam_formatted_file} | \
                /minimap2-2.17_x64-linux/k8 \
                /minimap2-2.17_x64-linux/paftools.js \
                sam2paf \
                -L \
                - \
                > ~{prefix}".paf"
        else
            echo "Currently we only support SAM or BAM (not CRAM)." && exit 1;
        fi
    >>>

    output {
        File pat_formatted_file = "~{prefix}.paf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
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
