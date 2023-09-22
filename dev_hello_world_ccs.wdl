# The most minimal workflow possible, suitable for testing call caching,
# using a ConsensusReadSet as input

version 1.0

task hello_world {
  input {
    File eid_ccs
    Int exit_code = 0
  }
  command {
    echo "Hello, world!" > hello.txt
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

workflow dev_hello_world_ccs {
  input {
    File eid_ccs
    Int exit_code = 0

    Int nproc = 1
    Int max_nchunks = 100
    String log_level = "INFO"
    String? tmp_dir
  }

  call hello_world {
    input:
      eid_ccs = eid_ccs,
      exit_code = exit_code
  }

  output {
    File hello_txt = hello_world.txt_out
  }
}
