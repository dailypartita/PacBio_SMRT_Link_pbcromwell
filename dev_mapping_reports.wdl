# workflow for generating mapping reports starting from existing dataset

version 1.0

import "tasks/pbreports.wdl" as pbreports
import "wf_coverage_reports.wdl"

workflow dev_mapping_reports {
  input {
    File? eid_alignment
    File? eid_ccs_alignment
    File eid_ref_dataset

    Int nproc = 1
    Int max_nchunks = 1
    String log_level = "INFO"
    String? tmp_dir
  }

  # workaround for input dataset type polymorphism
  Array[File?] datasets_ = [eid_alignment, eid_ccs_alignment]
  File mapped_dataset = select_first(datasets_)
  String mapping_stats_module = if (defined(eid_alignment)) then "mapping_stats_subreads" else "mapping_stats"

  call pbreports.mapping_stats as mapping_stats {
    input:
      mapped = mapped_dataset,
      report_module = mapping_stats_module,
      nproc = nproc,
      log_level = log_level
  }

  call wf_coverage_reports.coverage_reports {
    input:
      mapped = mapped_dataset,
      reference = eid_ref_dataset,
      nproc = nproc,
      log_level = log_level
  }

  output {
    File report_mapping_stats = mapping_stats.report
    File coverage_gff = coverage_reports.coverage_gff
    File report_coverage = coverage_reports.report_coverage
  }
}
