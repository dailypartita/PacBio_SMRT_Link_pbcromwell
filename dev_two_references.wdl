version 1.0

task echo_references {
  input {
    File ref1
    File ref2
  }
  command {
    echo "Ref1 is ${ref1}" > echo.txt
    echo "Ref2 is ${ref2}" >> echo.txt
  }
  runtime {
    cpu: 1
    memory: "50MB"
  }
  output {
    File txt_out = "echo.txt"
  }
}

workflow dev_two_references {
  input {
    File eid_ccs
    File eid_ref_dataset
    File eid_ref_dataset_2

    String log_level
    Int nproc = 1
    Int max_nchunks = 1
    String? tmp_dir
  }

  call echo_references {
    input:
      ref1 = eid_ref_dataset,
      ref2 = eid_ref_dataset_2
  }

  output {
    File echo_txt = echo_references.txt_out
  }
}
