# Simple wrapper for the delete-job utility

version 1.0

task delete_job {
  input {
    String job_dir
    Boolean remove_files = true

    String log_level = "INFO"
    Int nproc = 1
  }
  command {
    delete-job --log2stdout --log-level ${log_level} \
      ${false="--no-remove" true="" remove_files} \
      --report-json delete_job.report.json \
      ${job_dir}
  }
  runtime {
    cpu: 1
    memory: "50MB"
  }
  output {
    File report = "delete_job.report.json"
  }
}

workflow sl_delete_job {
  input {
    String path
    Boolean remove_files = true

    String log_level = "INFO"
    Int nproc = 1
    File? tmp_dir
    Int max_nchunks = 1
  }

  call delete_job {
    input:
      job_dir = path,
      remove_files = remove_files,
      log_level = log_level
  }

  output {
    File report = delete_job.report
  }
}
