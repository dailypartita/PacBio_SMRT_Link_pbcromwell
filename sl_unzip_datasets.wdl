# Very simple wrapper for the unzip-datasets CLI tool, so we can dump the
# blocking part of the job on Cromwell

version 1.0

import "sl_download_update.wdl"

task unzip_datasets {
  input {
    File zip_file

    String log_level = "INFO"
    Int nproc = 1
  }
  command {
    unzip-datasets --log2stdout --log-level ${log_level} ${zip_file}
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

workflow sl_unzip_datasets {
  input {
    File zip_file
    Boolean is_aws = false
    Boolean is_http = false

    String log_level = "INFO"
    Int nproc = 1
    Int max_nchunks = 1
    File tmp_dir = "/tmp"
  }

  if (is_aws) {
    call sl_download_update.download_s3 {
      input:
        file_uri = zip_file,
        tmp_dir = tmp_dir
    }
  }

  if (is_http) {
    call sl_download_update.download_curl {
      input:
        file_uri = zip_file,
        tmp_dir = tmp_dir
    }
  }

  File local_zip_file = select_first([download_s3.downloaded,
                                      download_curl.downloaded,
                                      zip_file])

  call unzip_datasets {
    input:
      zip_file = local_zip_file,
      log_level = log_level
  }

  output {
    File datastore = unzip_datasets.datastore
  }
}
