# CCS version of SAT workflow (no consensus step)

version 1.0

import "pb_align_ccs.wdl"

task report_site_acceptance_test_ccs {
  input {
    File mapped
    # see above
    File? coverage_report
    File? mapping_stats_report

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 \
      -m pbreports.report.sat_ccs \
      --log-level ${log_level} \
      ${mapped} \
      ${coverage_report} \
      ${mapping_stats_report} \
      site_acceptance_test_ccs.report.json
  }
  runtime {
    cpu: 1
    memory: "1GB"
  }
  output {
    File report = "site_acceptance_test_ccs.report.json"
  }
}

workflow pb_sat_ccs {
  input {
    File eid_ccs
    File eid_ref_dataset
    String dataset_filters = ""
    Int filter_min_qv = 20
    Int downsample_factor = 0

    String log_level = "INFO"
    Int nproc = 1
    Int max_nchunks = 1
    String? tmp_dir
  }

  call pb_align_ccs.pb_align_ccs {
    input:
      eid_ccs = eid_ccs,
      eid_ref_dataset = eid_ref_dataset,
      dataset_filters = dataset_filters,
      filter_min_qv = filter_min_qv,
      downsample_factor = downsample_factor,
      log_level = log_level,
      nproc = nproc,
      max_nchunks = max_nchunks,
      tmp_dir = tmp_dir
  }

  call report_site_acceptance_test_ccs {
    input:
      mapped = pb_align_ccs.mapped,
      mapping_stats_report = pb_align_ccs.report_mapping_stats,
      coverage_report = pb_align_ccs.report_coverage,
      log_level = log_level
  }

  output {
    File mapped = pb_align_ccs.mapped
    File? coverage_gff = pb_align_ccs.coverage_gff
    File? report_mapping_stats = pb_align_ccs.report_mapping_stats
    File? report_coverage = pb_align_ccs.report_coverage
    File? mapped_bam_datastore = pb_align_ccs.mapped_bam_datastore
    File report_sat_ccs = report_site_acceptance_test_ccs.report
  }
}
