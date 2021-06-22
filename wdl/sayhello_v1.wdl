version 1.0

workflow myWorkflow {
    call task_print
}

task task_print {
    command {
        echo "hello world" > out_file.txt
    }
    output {
        File out = "out_file.txt"
    }
    runtime {
   docker: 'ubuntu:latest'
   preemptible: 1
   memory: "4 GB"
  }
}









