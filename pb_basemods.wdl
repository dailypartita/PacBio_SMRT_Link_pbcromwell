# LEGACY WORKFLOW

version 1.0

import "wf_mapping_subreads.wdl"
import "wf_coverage_reports.wdl"
import "wf_basemods.wdl"
import "tasks/pbcoretools.wdl" as pbcoretools
import "tasks/chunking.wdl" as chunking

workflow pb_basemods {
  input {
    File eid_subread
    File eid_ref_dataset
    Boolean run_find_motifs = false

    String dataset_filters = ""
    Int downsample_factor = 0
    String kineticstools_identify_mods = "m4C,m6A"
    Float kineticstools_p_value = 0.001
    Int motif_min_score = 100
    Float motif_min_fraction = 0.30
    Float mapping_min_concordance = 70
    Int mapping_min_length = 50
    String? mapping_biosample_name
    String? mapping_pbmm2_overrides

    Int nproc = 1
    Int max_nchunks = 100
    Int target_size = 100000
    String log_level = "INFO"
    String? tmp_dir
  }

  call pbcoretools.dataset_filter {
    input:
      dataset = eid_subread,
      filters = dataset_filters,
      downsample_factor = downsample_factor,
      nproc = 1
  }

  call wf_mapping_subreads.mapping_subreads as mapping {
    input:
      reads = dataset_filter.filtered,
      reference = eid_ref_dataset,
      biosample_name = mapping_biosample_name,
      pbmm2_overrides = mapping_pbmm2_overrides,
      min_concordance = mapping_min_concordance,
      min_length = mapping_min_length,
      run_coverage = false,
      target_size = target_size,
      max_nchunks = max_nchunks,
      nproc = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir
  }

  call wf_coverage_reports.coverage_reports {
    input:
      mapped = mapping.mapped,
      reference = eid_ref_dataset,
      nproc = nproc,
      log_level = log_level
  }

  call wf_basemods.wf_basemods {
    input:
      alignments = mapping.mapped,
      reference = eid_ref_dataset,
      kineticstools_identify_mods = kineticstools_identify_mods,
      kineticstools_p_value = kineticstools_p_value,
      run_find_motifs = run_find_motifs,
      motif_min_score = motif_min_score,
      motif_min_fraction = motif_min_fraction,
      genome_length_mb = mapping.genome_length_mb,
      max_nchunks = max_nchunks,
      target_size = target_size,
      nproc = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir
  }

  output {
    File mapped = mapping.mapped
    File? report_mapping_stats = mapping.report_mapping_stats
    File coverage_gff = coverage_reports.coverage_gff
    File report_coverage = coverage_reports.report_coverage
    File basemods_gff = wf_basemods.basemods_gff
    File? basemods_csv = wf_basemods.basemods_csv
    File? report_modifications = wf_basemods.report_modifications
    File? bigwig_file = wf_basemods.bigwig_file
    File? motifs_csv = wf_basemods.motifs_csv
    File? motifs_gff = wf_basemods.motifs_gff
    File? report_motifs = wf_basemods.report_motifs
    File? mapped_bam_datastore = mapping.mapped_bam_datastore
  }
}
