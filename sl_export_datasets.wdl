# workflow for export-datasets services job type

version 1.0

task export_datasets {
  input {
    Array[File] datasets
    # the next two options are mutually exclusive
    String? zip_out_path
    String? zip_output_dir
    File? json_summaries
    Boolean? mock_uuid

    String log_level = "INFO"
    Int nproc = 1
  }
  command <<<
    set -e
    export-datasets \
      --log2stdout \
      --log-level ~{log_level} \
      ~{"-o " + zip_out_path} \
      ~{"--output-dir " + zip_output_dir} \
      ~{true="--mock-uuid" false="" mock_uuid} \
      ~{"--reports-json " + json_summaries} \
      --symlink exported_datasets.zip \
      ~{sep=" " datasets}
  >>>
  runtime {
    cpu: 1
    memory: "1GB"
  }
  output {
    File zip_out = "exported_datasets.zip"
  }
}

workflow sl_export_datasets {
  input {
    Array[File] datasets
    # optional serialized dataset summary reports, which will be turned into
    # PDFs and added to the zip archive
    File? json_summaries
    # unlike every other Cromwell workflow except sl_export_job, we allow
    # the user to select an external output directory and write to it from
    # within the task
    String? zip_out_path
    # for internal use
    Boolean? mock_uuid = false

    String log_level = "INFO"
    File? tmp_dir
    Int nproc = 1
    Int max_nchunks = 1
  }

  call export_datasets {
    input:
      datasets = datasets,
      zip_out_path = zip_out_path,
      json_summaries = json_summaries,
      mock_uuid = mock_uuid,
      log_level = log_level
  }

  output {
    File zip_out = export_datasets.zip_out
  }
}
