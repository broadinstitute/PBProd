version 1.0

workflow count_records {
    input{
        File bamfile
    }
    call count{
        input: bamfile = bamfile

    }
    output{
        File out_wf = count.out
    }
}

task count {
    input{
        File bamfile
    }
    command {
        set -euxo pipefail
        samtools view -c "~{bamfile}" > "counts_out.txt"
    }
    output {
        File out = "counts_out.txt"
    }
    runtime {
           docker: "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
           preemptible: 1
           memory: "4 GB"
    }

}