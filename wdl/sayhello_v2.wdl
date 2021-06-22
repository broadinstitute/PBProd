version 1.0

task printwords {
  input{
    String sayhello
  }

  command {
    echo “~{sayhello}”
  }
  output {
    String hello_out = read_string(stdout())
  }
  runtime {
   docker: 'ubuntu:latest'
   preemptible: 0
   memory: "4 GB"
  }
}

workflow helloworld {
  input{
    String sayhello
  }
  call printwords as printwords1{
    input: sayhello = sayhello

  }

  call printwords as printwords2{
    input: sayhello = sayhello + ",earth!"

  }
  output {
    String out1 = printwords1.hello_out
    String out2 = printwords2.hello_out


  }

}







