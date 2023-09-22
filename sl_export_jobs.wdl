# workflow for export-job services job type

version 1.0

import "sl_export_datasets.wdl"

task export_job {
  input {
    File job_json
    String? zip_output_dir

    String log_level = "INFO"
    Int nproc = 1
  }
  command {
    set -e
    export-job \
      --log2stdout \
      --log-level ${log_level} \
      ${"--output-dir " + zip_output_dir} \
      --symlink exported_job.zip \
      ${job_json}
  }
  runtime {
    cpu: 1
    memory: "1GB"
  }
  output {
    File zip_out = "exported_job.zip"
  }
}

workflow sl_export_jobs {
  input {
    Array[File] jobs_json
    # this is a directory, not an actual zip file
    String? zip_out_path
    Array[File]? datasets

    String log_level = "INFO"
    File? tmp_dir
    Int nproc = 1
    Int max_nchunks = 1
  }

  scatter (job in jobs_json) {
    call export_job {
      input:
        job_json = job,
        zip_output_dir = zip_out_path,
        log_level = log_level
    }
  }

  if (defined(datasets)) {
    call sl_export_datasets.export_datasets {
      input:
        datasets = select_first([datasets]),
        zip_output_dir = zip_out_path,
        log_level = log_level
    }
  }

  output {
    Array[File] zip_out = export_job.zip_out
    File? datasets_zip_out = export_datasets.zip_out
  }
}
