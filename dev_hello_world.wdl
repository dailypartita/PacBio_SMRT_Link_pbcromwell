# The most minimal workflow possible, suitable for testing call caching

version 1.0

task hello_world {
  input {
    File eid_subread
    String hello_txt = "Hello, world!"
    Int exit_code = 0
  }
  command {
    echo "${hello_txt}" > hello.txt
    exit ${exit_code}
  }
  runtime {
    backend: "Local"
    cpu: 1
    memory: "1GB"
  }
  output {
    File txt_out = "hello.txt"
  }
}

workflow dev_hello_world {
  input {
    File eid_subread
    Int exit_code = 0
    String hello_txt = "Hello, world!"

    Int nproc = 1
    Int max_nchunks = 100
    String log_level = "INFO"
    String? tmp_dir
  }

  call hello_world {
    input:
      eid_subread = eid_subread,
      exit_code = exit_code
  }

  output {
    File hello_txt = hello_world.txt_out
  }
}
