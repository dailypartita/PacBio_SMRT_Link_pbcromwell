# Very simple wrapper for the unzip-datasets CLI tool, so we can dump the
# blocking part of the job on Cromwell

version 1.0

task unzip_job {
  input {
    File zip_file
    String job_dir

    String log_level = "INFO"
    Int nproc = 1
  }
  command {
    set -e
    unzip-job --log2stdout --log-level ${log_level} ${zip_file} ${job_dir}
    ln -s ${job_dir}/workflow/datastore.json
  }
  runtime {
    memory: "50MB"
    cpu: 1
    backend: "Local"
  }
  output {
    # this should have been included in the zip file when unpacking to PWD
    File datastore = "datastore.json"
  }
}

workflow sl_unzip_job {
  input {
    File zip_file
    # this can't be a File because Cromwell will try to localize it
    String job_dir

    String log_level = "INFO"
    Int nproc = 1
    Int max_nchunks = 1
    File? tmp_dir
  }

  call unzip_job {
    input:
      zip_file = zip_file,
      job_dir = job_dir,
      log_level = log_level
  }

  output {
    File datastore = unzip_job.datastore
  }
}
