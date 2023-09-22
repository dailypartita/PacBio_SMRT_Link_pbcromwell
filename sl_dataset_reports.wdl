# Standalone dataset report generation, for both subreads and CCS

version 1.0

import "tasks/pbcoretools.wdl" as pbcoretools

# This replaces the formerly separate tasks for individual reports
task import_dataset_reports {
  input {
    File dataset_xml

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 \
      -m pbreports.tasks.import_dataset_reports \
      --log-level ${log_level} \
      "`readlink -f '${dataset_xml}'`"
  }
  runtime {
    cpu: 1
    memory: "12GB"
  }
  output {
    File? report_adapter = "adapter.report.json"
    File? report_loading = "loading.report.json"
    File? report_raw_data = "raw_data.report.json"
    File? report_control = "control.report.json"
    File? report_subread_stats = "subread_stats.report.json"
    File? report_ccs2 = "ccs.report.json"
    File? report_detect_cpg_methyl = "detect_cpg_methyl.report.json"
    Array[File?] plot_pngs = glob("*.png")
  }
}

workflow sl_dataset_reports {
  input {
    File? eid_subread
    File? eid_ccs

    Int nproc = 1
    String log_level = "INFO"
    String? tmp_dir
    Int max_nchunks = 1
  }

  Array[File?] datasets_ = [eid_subread, eid_ccs]
  File dataset = select_first(datasets_)

  call import_dataset_reports {
    input:
      dataset_xml = dataset,
      nproc = nproc,
      log_level = log_level
  }

  output {
    File? report_raw_data = import_dataset_reports.report_raw_data
    File? report_adapters = import_dataset_reports.report_adapter
    File? report_loading = import_dataset_reports.report_loading
    File? report_control = import_dataset_reports.report_control
    File? report_subread_stats = import_dataset_reports.report_subread_stats
    File? report_ccs2 = import_dataset_reports.report_ccs2
    File? report_detect_cpg_methyl = import_dataset_reports.report_detect_cpg_methyl
  }
}
