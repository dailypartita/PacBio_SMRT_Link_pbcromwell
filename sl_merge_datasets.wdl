# Workflow for running auto-merge

version 1.0

import "tasks/sl_tasks.wdl" as sl_tasks
import "sl_dataset_reports.wdl" as wf_dataset_reports

task dataset_merge {
  input {
    Array[File] datasets
    String dataset_ext
    String dataset_name = "Auto-merged input datasets"
    Boolean skip_counts = false

    Int nproc = 1
    String log_level = "DEBUG"
  }

  command {
    dataset \
      --log-level ${log_level} \
      ${true="--skipCounts" false="" skip_counts} \
      merge \
      --remove-parentage \
      --unique-collections \
      --no-sub-datasets \
      --name "${dataset_name}" \
      merged${dataset_ext} \
      ${sep=" " datasets} \
    && \
      dataset newuuid --random merged${dataset_ext}
  }
  runtime {
    cpu: 1
  }
  output {
    File merged = "merged${dataset_ext}"
  }
}

workflow sl_merge_datasets {
  input {
    Array[File] datasets
    String dataset_ext
    String dataset_name = "Auto-merged input datasets"
    # this will work for merge mode as long as the input datasets have the
    # correct counts already (which is true for anything that the instrument
    # or SMRT Link generates)
    Boolean skip_counts = true
    Boolean make_reports = true

    Int nproc = 1
    String log_level = "DEBUG"
    String? tmp_dir
  }

  call dataset_merge {
    input:
      datasets = datasets,
      dataset_ext = dataset_ext,
      dataset_name = dataset_name,
      skip_counts = skip_counts,
      nproc = nproc,
      log_level = log_level
  }

  if (make_reports) {
    call wf_dataset_reports.import_dataset_reports {
      input:
        dataset_xml = dataset_merge.merged,
        log_level = log_level,
        nproc = 1
    }
  }

  output {
    File merged = dataset_merge.merged
    File? report_raw_data = import_dataset_reports.report_raw_data
    File? report_adapters = import_dataset_reports.report_adapter
    File? report_loading = import_dataset_reports.report_loading
    File? report_control = import_dataset_reports.report_control
    File? report_subread_stats = import_dataset_reports.report_subread_stats
    File? report_ccs2 = import_dataset_reports.report_ccs2
  }
}
