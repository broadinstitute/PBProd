version 1.0

import "Structs.wdl"

task PB10xMasSeqSingleFlowcellReport {

    meta {
        description : "Create a report for a given MASSeq run which summarizes the results  using a given Jupyter Notebook template."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File notebook_template

        File subreads_stats
        File ccs_reads_stats
        File array_elements_stats
        File ccs_report_file

        File ccs_bam_file
        File array_element_bam_file
        File ccs_rejected_bam_file

        File annotated_bam_file

        File zmw_subread_stats_file
        File polymerase_read_lengths_file
        File approx_raw_subread_array_lengths

        File? ten_x_metrics_file

        File workflow_dot_file

        String prefix = ""
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        notebook_template : "Jupyter notebook MASSeq template to run with the given data to produce a MASSeq report."

        subreads_stats : "Samtools stats file created from the raw subreads from the PacBio instrument."
        ccs_reads_stats : "Samtools raw stats file created from the aligned CCS corrected reads from the PacBio instrument."
        array_elements_stats : "Samtools raw stats file created from the aligned MASSeq array elements."

        ccs_report_file : "CCS report file from the CCS run for the data from the PacBio instrument."
        ccs_bam_file : "Unaligned reads file in BAM format from the CCS process (pre-array splitting)."

        array_element_bam_file : "Aligned reads file in BAM format containing aligned MASSeq array elements as individual reads."
        ccs_rejected_bam_file : "Bam file containing all subreads from zmws that were rejected by CCS."

        annotated_bam_file : "Bam file containing ccs corrected reads with annotated sections in the SG tag."

        zmw_subread_stats_file : "File containing statistics about the subreads from each ZMW (created by collect_zmw_subread_stats.py in the PBUtils docker container)."
        polymerase_read_lengths_file : "File containing the lengths of each polymerase read from the sequencer (as created by collect_polymerase_read_lengths.py)"
        approx_raw_subread_array_lengths : "File containing the approximate array length information from the raw (pre-ccs) subreads file  (created by get_approx_raw_subread_array_lengths.py in the Cartographer docker container)."

        ten_x_metrics_file : "Stats file from the 10x tool run for the data in this MASSeq run.  If not supplied stats will not be displayed in the resulting report."
        workflow_dot_file : "DOT file containing the representation of this workflow used to create and analyze the data.  This is included in the QC reports (the DOT file can be generated with womtool)."

        prefix : "[optional] Prefix to prepend to the name of the generated report."
        runtime_attr_override : "[optional] Runtime attributes struct with which to override the docker container runtime.."
    }

    String nb_name = prefix + "report.ipynb"
    String html_out = prefix + "report.html"
    String pdf_out = prefix + "report.pdf"

    Int disk_size = 20 + 8*ceil((
            size(notebook_template, "GB") +
            size(subreads_stats, "GB") +
            size(ccs_reads_stats, "GB") +
            size(ccs_report_file, "GB") +
            size(ccs_bam_file, "GB") +
            size(array_element_bam_file, "GB") +
            size(ccs_rejected_bam_file, "GB") +
            size(annotated_bam_file, "GB") +
            size(ccs_bam_file, "GB") +
            size(zmw_subread_stats_file, "GB") +
            size(polymerase_read_lengths_file, "GB") +
            size(ten_x_metrics_file, "GB") +
            size(workflow_dot_file, "GB")
        ))

    # Handle the optional 10x metrics file:
    String ten_x_metrics_file_name = if defined(ten_x_metrics_file) then ten_x_metrics_file else "NON-EXISTENT-PLACEHOLDER"

    command <<<
        set -euxo pipefail

        # Copy the notebook template to our current folder:
        cp ~{notebook_template} ~{nb_name}

        # Create a template to create the html report with collapsed code:
        echo "{%- extends 'full.tpl' -%}" > hidecode.tpl
        echo "" >> hidecode.tpl
        echo "{% block input_group %}" >> hidecode.tpl
        echo "    {%- if cell.metadata.get('nbconvert', {}).get('show_code', False) -%}" >> hidecode.tpl
        echo "        ((( super() )))" >> hidecode.tpl
        echo "    {%- endif -%}" >> hidecode.tpl
        echo "{% endblock input_group %}" >> hidecode.tpl

        # Set some environment variables for the notebook to read in:
        export DATE_RUN="$(date)"
        export WDL_NAME="PB10xMasSeqArraySingleFlowcell.wdl"
        export REPO_INFO="git@github.com:broadinstitute/long-read-pipelines.git"

        # Prepare the config file:
        rm -f mas-seq_qc_inputs.config

        echo "~{subreads_stats}" >> mas-seq_qc_inputs.config
        echo "~{ccs_reads_stats}" >> mas-seq_qc_inputs.config
        echo "~{array_elements_stats}" >> mas-seq_qc_inputs.config
        echo "~{ccs_report_file}" >> mas-seq_qc_inputs.config

        echo "~{ccs_bam_file}" >> mas-seq_qc_inputs.config
        echo "~{array_element_bam_file}" >> mas-seq_qc_inputs.config
        echo "~{ccs_rejected_bam_file}" >> mas-seq_qc_inputs.config

        echo "~{annotated_bam_file}" >> mas-seq_qc_inputs.config

        echo "~{zmw_subread_stats_file}" >> mas-seq_qc_inputs.config
        echo "~{polymerase_read_lengths_file}" >> mas-seq_qc_inputs.config
        echo "~{approx_raw_subread_array_lengths}" >> mas-seq_qc_inputs.config

        echo "~{ten_x_metrics_file_name}" >> mas-seq_qc_inputs.config

        echo "~{workflow_dot_file}" >> mas-seq_qc_inputs.config

        # Do the conversion:

        # Run the notebook and populate the notebook itself:
        jupyter nbconvert --execute ~{nb_name} --to notebook --inplace --no-prompt --no-input --clear-output --debug --ExecutePreprocessor.timeout=3600

        # Convert the notebook output we created just above here to the HTML report:
        jupyter nbconvert ~{nb_name} --to html --no-prompt --no-input --debug --ExecutePreprocessor.timeout=3600

        # One more for good measure - make a PDF so we don't need to wait for the browser all the time.
        jupyter nbconvert ~{nb_name} --to pdf --no-prompt --no-input --debug --ExecutePreprocessor.timeout=3600
    >>>

    output {
        File populated_notebook = nb_name
        File html_report = html_out
        File pdf_report = pdf_out
        File generated_config = "mas-seq_qc_inputs.config"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-jupyter_interactive:0.0.1"
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
