version 1.0

task pbreports_jobs_report {
  input {
    File jobs_json

    String log_level = "INFO"
  }
  command {
    python3 -m pbreports.tasks.jobs_meta_report \
      --log-level ${log_level} \
      --log-file pbreports.log \
      --datastore-json reports.datastore.json \
      ${jobs_json}
  }
  runtime {
    cpu: 1
    memory: "200MB"
  }
  output {
    File datastore = "reports.datastore.json"
  }
}

workflow sl_jobs_report {
  input {
    # this is its own special format
    File jobs_json

    Int nproc = 0
    Int max_nchunks = 0
    String log_level = "INFO"
    String tmp_dir = "/tmp"
  }

  call pbreports_jobs_report {
    input:
      jobs_json = jobs_json,
      log_level = log_level
  }

  output {
    File reports_datastore = pbreports_jobs_report.datastore
  }
}
