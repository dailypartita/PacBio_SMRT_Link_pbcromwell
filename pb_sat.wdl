# SAT workflow using pbmm2

version 1.0

import "pb_resequencing.wdl" as resequencing_wdl

task report_site_acceptance_test {
  input {
    File mapped
    # XXX this is hacky, a side effect of making report generation optional
    # to accommodate sub-workflow re-usability
    File? variants_report
    File? mapping_stats_report

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 \
      -m pbreports.report.sat \
      --log-level ${log_level} \
      ${mapped} \
      ${variants_report} \
      ${mapping_stats_report} \
      site_acceptance_test.report.json
  }
  runtime {
    cpu: 1
    # only the basic XML metadata are loaded now
    memory: "1GB"
  }
  output {
    File report = "site_acceptance_test.report.json"
  }
}

workflow pb_sat {
  input {
    File eid_subread
    File eid_ref_dataset
    String dataset_filters = ""
    Int downsample_factor = 0

    Int nproc = 1
    Int max_nchunks = 100
    String log_level = "INFO"
    String? tmp_dir
  }

  String CONSENSUS_ALGORITHM = "arrow"

  call resequencing_wdl.pb_resequencing {
    input:
      eid_subread = eid_subread,
      eid_ref_dataset = eid_ref_dataset,
      dataset_filters = dataset_filters,
      consensus_algorithm = CONSENSUS_ALGORITHM,
      downsample_factor = downsample_factor,
      nproc = nproc,
      max_nchunks = max_nchunks,
      log_level = log_level,
      tmp_dir = tmp_dir
  }

  call report_site_acceptance_test as sat_report {
    input:
      mapped = pb_resequencing.mapped,
      variants_report = pb_resequencing.report_variants,
      mapping_stats_report = pb_resequencing.report_mapping_stats,
      log_level = log_level
  }

  output {
    File report_sat = sat_report.report
    File mapped = pb_resequencing.mapped
    File coverage_gff = pb_resequencing.coverage_gff
    File? report_mapping_stats = pb_resequencing.report_mapping_stats
    File report_coverage = pb_resequencing.report_coverage
    File consensus_fasta = pb_resequencing.consensus_fasta
    File consensus_fastq = pb_resequencing.consensus_fastq
    File variants_gff = pb_resequencing.variants_gff
    File variants_vcf = pb_resequencing.variants_vcf
    File consensus_gff = pb_resequencing.consensus_gff
    File report_variants = pb_resequencing.report_variants
    File report_top_variants = pb_resequencing.report_top_variants
  }
}
