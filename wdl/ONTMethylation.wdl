version 1.0

import "tasks/Utils.wdl" as Utils
import "tasks/ONTUtils.wdl" as ONTUtils
import "tasks/Guppy.wdl" as Guppy
import "tasks/Finalize.wdl" as FF

workflow ONTMethylation {
    input {
        String gcs_fast5_dir

        File ref_map_file
        File variants

        #String prefix

        #String gcs_out_root_dir
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    #String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTMethylation/~{prefix}"

    call Utils.ListFilesOfType { input: gcs_dir = gcs_fast5_dir, suffixes = [ 'fast5' ] }

    scatter (fast5_file in ListFilesOfType.files) {
        call Megalodon {
            input:
                fast5_files = [ fast5_file ],
                ref_fasta   = ref_map['fasta'],
                variants    = variants
        }

        call Utils.SortBam as SortMappings { input: input_bam = Megalodon.mappings_bam }
        call Utils.SortBam as SortModMappings { input: input_bam = Megalodon.mod_mappings_bam }
        call Utils.SortBam as SortVarMappings { input: input_bam = Megalodon.variant_mappings_bam }
    }

    call Merge as MergeVariantDBs {
        input:
            dbs = Megalodon.per_read_variant_calls_db,
            merge_type = "variants",
            runtime_attr_override = { 'mem_gb': 48 }
    }

    call Merge as MergeModifiedBaseCallDBs {
        input:
            dbs = Megalodon.per_read_modified_base_calls_db,
            merge_type = "modified_bases"
    }

    call Utils.MergeFastqGzs { input: fastq_gzs = Megalodon.basecalls_fastq, prefix = "basecalls" }

    call Utils.MergeBams as MergeMappings { input: bams = SortMappings.sorted_bam }
    call Utils.MergeBams as MergeModMappings { input: bams = SortModMappings.sorted_bam }
    call Utils.MergeBams as MergeVarMappings { input: bams = SortVarMappings.sorted_bam }

    call Utils.Cat as CatModifiedBases5mC {
         input:
            files = Megalodon.modified_bases_5mC,
            has_header = false,
            out = "modified_bases.5mC.bed"
    }

    call Utils.Cat as CatMappingSummaries {
        input:
            files = Megalodon.mappings_summary,
            has_header = true,
            out = "mappings_summary.txt"
    }

    call Utils.Cat as CatSequencingSummaries {
        input:
            files = Megalodon.sequencing_summary,
            has_header = true,
            out = "sequencing_summary.txt"
    }

    output {
        #String gcs_basecall_dir = Guppy.gcs_dir
    }
}

task Megalodon {
    input {
        Array[File] fast5_files
        File ref_fasta
        File variants

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4 * ceil(size(fast5_files, "GB") + size([ref_fasta, variants], "GB"))

    command <<<
        set -x

        num_cores=$(grep -c '^processor' /proc/cpuinfo | awk '{ print $1 - 1 }')
        dir=$(dirname ~{fast5_files[0]})

        megalodon $dir \
            --guppy-params "-d /rerio/basecall_models" \
            --guppy-config res_dna_r941_prom_modbases_5mC_CpG_v001.cfg \
            --outputs basecalls mappings mod_mappings mods variant_mappings \
            --reference ~{ref_fasta} \
            --mod-motif m CG 0 \
            --variant-filename ~{variants} \
            --devices cuda:0 \
            --processes $num_cores \
            --guppy-server-path /usr/bin/guppy_basecall_server
    >>>

    output {
        File basecalls_fastq = "megalodon_results/basecalls.fastq"
        File log = "megalodon_results/log.txt"
        File mappings_bam = "megalodon_results/mappings.bam"
        File mappings_summary = "megalodon_results/mappings.summary.txt"
        File mod_mappings_bam = "megalodon_results/mod_mappings.bam"
        File variant_mappings_bam = "megalodon_results/variant_mappings.bam"
        File modified_bases_5mC = "megalodon_results/modified_bases.5mC.bed"
        File per_read_modified_base_calls_db = "megalodon_results/per_read_modified_base_calls.db"
        File per_read_variant_calls_db = "megalodon_results/per_read_variant_calls.db"
        File sequencing_summary = "megalodon_results/sequencing_summary.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             50,
        disk_gb:            disk_size,
        boot_disk_gb:       30,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-megalodon:2.3.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        gpuType:                "nvidia-tesla-p100"
        gpuCount:               1
        nvidiaDriverVersion:    "418.152.00"
        zones:                  ["us-central1-c", "us-central1-f", "us-east1-b", "us-east1-c", "us-west1-a", "us-west1-b"]
        cpuPlatform:            "Intel Haswell"
    }
}

task Merge {
    input {
        Array[File] dbs
        String merge_type

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(dbs, "GB")) + 1

    command <<<
        set -x

        DIRS=$(find /cromwell_root/ -name '*.db' -exec dirname {} \; | tr '\n' ' ')

        megalodon_extras merge ~{merge_type} $DIRS
    >>>

    output {
        File? per_read_variant_calls_db = "megalodon_merge_vars_results/per_read_variant_calls.db"
        File? per_read_modified_base_calls_db = "megalodon_merge_mods_results/per_read_modified_base_calls.db"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-megalodon:2.3.1"
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
