version 1.0

task greeting {

  command {
    echo 'hello world!'
  }
  output {
    File hello_out = stdout()
  }
  runtime {
   docker: 'ubuntu:latest'
  }
}

workflow test {
  call greeting
}








