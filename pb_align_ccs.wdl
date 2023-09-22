# Alignment starting from ConsensusReadSet

version 1.0

import "wf_mapping.wdl"
import "wf_prepare_input.wdl"

workflow pb_align_ccs {
  input {
    File eid_ccs
    File eid_ref_dataset
    File? target_regions_bed
    Int filter_min_qv = 20
    String dataset_filters = ""
    Float mapping_min_concordance = 70
    Int mapping_min_length = 50
    String? mapping_biosample_name
    String? mapping_pbmm2_overrides
    Int downsample_factor = 0
    Boolean? report_show_calibration_plot
    Int add_memory_mb = 0

    Int nproc = 1
    # this is no longer used
    Int max_nchunks = 0
    String log_level = "INFO"
    String? tmp_dir
  }
  Int BASE_MEMORY_APP = 512
  Int base_memory_mb = BASE_MEMORY_APP + add_memory_mb

  call wf_prepare_input.prepare_input {
    input:
      dataset_xml = eid_ccs,
      dataset_ext = ".consensusreadset.xml",
      dataset_filters = dataset_filters,
      downsample_factor = downsample_factor,
      filter_min_qv = filter_min_qv,
      nproc = 1,
      log_level = log_level
  }

  call wf_mapping.mapping {
    input:
      reads = prepare_input.reads_file,
      reference = eid_ref_dataset,
      target_regions_bed = target_regions_bed,
      preset_mode = "HiFi",
      min_concordance = mapping_min_concordance,
      min_length = mapping_min_length,
      biosample_name = mapping_biosample_name,
      pbmm2_overrides = mapping_pbmm2_overrides,
      show_calibration_plot = report_show_calibration_plot,
      mem_scale_factor = 6,
      base_memory_mb = base_memory_mb,
      nproc = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir
  }

  output {
    File mapped = mapping.mapped
    File? coverage_gff = mapping.coverage_gff
    File? report_mapping_stats = mapping.report_mapping_stats
    File? report_coverage = mapping.report_coverage
    File? report_target_coverage = mapping.report_target_coverage
    File? mapped_bam_datastore = mapping.mapped_bam_datastore
    String? mapped_bam = mapping.mapped_bam
    String? mapped_bam_bai = mapping.mapped_bam_bai
  }
}
