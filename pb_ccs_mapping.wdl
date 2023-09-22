version 1.0

import "pb_ccs.wdl" as ccs_wdl
import "wf_mapping.wdl"
import "tasks/pbcoretools.wdl" as pbcoretools

workflow pb_ccs_mapping {
  input {
    File eid_subread
    File eid_ref_dataset
    File? target_regions_bed
    Float ccs_min_predicted_accuracy = 0.99
    Int ccs_min_passes = 3
    Float ccs_min_snr = 2.5
    Int ccs_min_length = 10
    Int ccs_max_length = 50000
    Boolean ccs_polish = true
    String? ccs_model_args
    String? ccs_override_args
    Int downsample_factor = 0
    Int filter_min_qv = 20
    Float mapping_min_concordance = 70
    Int mapping_min_length = 50
    String? mapping_biosample_name
    String? mapping_pbmm2_overrides
    Boolean? report_show_calibration_plot
    Boolean ccs_by_strand = false
    Boolean ccs_process_all = false
    Boolean ccs_include_kinetics = false

    Int nproc = 1
    Int max_nchunks = 300
    String log_level = "INFO"
    String? tmp_dir
  }

  call ccs_wdl.pb_ccs {
    input:
      eid_subread = eid_subread,
      ccs_min_passes = ccs_min_passes,
      ccs_min_predicted_accuracy = ccs_min_predicted_accuracy,
      ccs_min_snr = ccs_min_snr,
      ccs_min_length = ccs_min_length,
      ccs_max_length = ccs_max_length,
      ccs_polish = ccs_polish,
      ccs_model_args = ccs_model_args,
      ccs_process_all = ccs_process_all,
      ccs_include_kinetics = ccs_include_kinetics,
      ccs_override_args = ccs_override_args,
      ccs_by_strand = ccs_by_strand,
      nproc = nproc,
      max_nchunks = max_nchunks,
      log_level = log_level
  }

  call pbcoretools.dataset_filter {
    input:
      dataset = pb_ccs.ccsxml,
      dataset_ext = ".consensusreadset.xml",
      filters = "",
      downsample_factor = downsample_factor,
      min_qv = filter_min_qv,
      log_level = log_level
  }

  call wf_mapping.mapping as ccs_mapping {
    input:
      reads = dataset_filter.filtered,
      reference = eid_ref_dataset,
      target_regions_bed = target_regions_bed,
      preset_mode = "CCS",
      min_concordance = mapping_min_concordance,
      min_length = mapping_min_length,
      biosample_name = mapping_biosample_name,
      pbmm2_overrides = mapping_pbmm2_overrides,
      show_calibration_plot = report_show_calibration_plot,
      mem_scale_factor = 6,
      nproc = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir
  }

  output {
    File ccsxml = pb_ccs.ccsxml
    File report_ccs2 = pb_ccs.report_ccs2
    File report_ccs_processing = pb_ccs.report_ccs_processing
    File ccs_csv = pb_ccs.ccs_csv
    File ccs_zmws = pb_ccs.ccs_zmws
    File bam_datastore = pb_ccs.bam_datastore
    File fasta_datastore = pb_ccs.fasta_datastore
    File fastq_datastore = pb_ccs.fastq_datastore
    File mapped = ccs_mapping.mapped
    # these will always be generated here but need optional types because of
    # mapping workflow structure
    File? coverage_gff = ccs_mapping.coverage_gff
    File? report_mapping_stats = ccs_mapping.report_mapping_stats
    File? report_coverage = ccs_mapping.report_coverage
    File? report_target_coverage = ccs_mapping.report_target_coverage
    # these are truly optional
    File? mapped_bam_datastore = ccs_mapping.mapped_bam_datastore
    String? reads_bam = pb_ccs.reads_bam
    String? mapped_bam = ccs_mapping.mapped_bam
    String? mapped_bam_bai = ccs_mapping.mapped_bam_bai
  }
}
