# Simple wrapper for the delete-datasets utility

version 1.0

task delete_datasets {
  input {
    File xml_or_fofn
    Boolean remove_files = true

    String log_level = "INFO"
    Int nproc = 1
  }
  command {
    delete-datasets --log2stdout --log-level ${log_level} \
      ${false="--no-remove" true="" remove_files} \
      --report-json delete_datasets.report.json \
      ${xml_or_fofn}
  }
  runtime {
    cpu: 1
    memory: "50MB"
  }
  output {
    File report = "delete_datasets.report.json"
  }
}

workflow sl_delete_datasets {
  input {
    File path

    String log_level = "INFO"
    Int nproc = 1
    File? tmp_dir
    Int max_nchunks = 1
  }

  call delete_datasets {
    input:
      xml_or_fofn = path,
      log_level = log_level
  }

  output {
    File report = delete_datasets.report
  }
}
