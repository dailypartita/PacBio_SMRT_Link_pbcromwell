# elisgnment starting from SubreadSet

version 1.0

import "wf_mapping_subreads.wdl"
import "tasks/pbcoretools.wdl" as pbcoretools

workflow pb_align_subreads {
  input {
    File eid_subread
    File eid_ref_dataset
    File? target_regions_bed
    String dataset_filters = ""
    Int downsample_factor = 0
    Float mapping_min_concordance = 70
    Int mapping_min_length = 50
    Boolean mapping_hq_mode = false
    Boolean mapping_zmw_mode = false
    String mapping_preset_mode = "SUBREAD"
    String? mapping_biosample_name
    String? mapping_pbmm2_overrides
    Int add_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
    Int max_nchunks = 1
    String? tmp_dir
  }

  call pbcoretools.dataset_filter {
    input:
      dataset = eid_subread,
      filters = dataset_filters,
      downsample_factor = downsample_factor,
      nproc = 1,
      log_level = log_level
  }

  call wf_mapping_subreads.mapping_subreads as mapping {
    input:
      reads = dataset_filter.filtered,
      reference = eid_ref_dataset,
      target_regions_bed = target_regions_bed,
      min_concordance = mapping_min_concordance,
      min_length = mapping_min_length,
      hq_mode = mapping_hq_mode,
      zmw_mode = mapping_zmw_mode,
      biosample_name = mapping_biosample_name,
      pbmm2_overrides = mapping_pbmm2_overrides,
      preset_mode = mapping_preset_mode,
      add_memory_mb = add_memory_mb,
      target_size = 100000,
      max_nchunks = max_nchunks,
      nproc = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir
  }

  output {
    File mapped = mapping.mapped
    File? report_mapping_stats = mapping.report_mapping_stats
    File? coverage_gff = mapping.coverage_gff
    File? report_coverage = mapping.report_coverage
    File? report_target_coverage = mapping.report_target_coverage
    File? mapped_bam_datastore = mapping.mapped_bam_datastore
  }
}
