# HiFi Basemods analysis using ccs-kinetics-bystrandify to generate mock
# subread-like input that can be processed by the CLR workflow

version 1.0

import "wf_prepare_input.wdl"
import "wf_basemods.wdl"
import "wf_coverage_reports.wdl"
import "tasks/pbmm2.wdl"
import "tasks/pbreports.wdl"

task ccs_kinetics_bystrandify {
  input {
    File ccs_xml
    Int base_memory_mb = 0

    String log_level = "INFO"
  }
  command {
    set -e
    ccs-kinetics-bystrandify \
      --log-level ${log_level} \
      --log-file ccs-kinetics-bystrandify.log \
      `readlink -f "${ccs_xml}"` \
      ccs_kinetics_bystrandify.subreads.bam
    # NOTE ccs_kinetics_bystrandify.subreads.bam is actually a dataset XML
    # ccs_kinetics_bystrandify.subreads.bystrand.bam is what we want
    dataset \
      --log-level ${log_level} \
      --log-file dataset-create.log \
      --strict \
      create \
      --type SubreadSet \
      --name "ccs-kinetics-bystrandify output" \
      ccs_kinetics_bystrandify.subreadset.xml \
      ccs_kinetics_bystrandify.subreads.bystrand.bam
  }
  Int total_mem_mb = 8192 + base_memory_mb
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File subreads_xml = "ccs_kinetics_bystrandify.subreadset.xml"
  }
}

workflow pb_basemods_hifi {
  input {
    File eid_ccs
    File eid_ref_dataset
    Boolean run_find_motifs = false

    Int filter_min_qv = 20
    String dataset_filters = ""
    String kineticstools_identify_mods = "m4C,m6A"
    Float kineticstools_p_value = 0.001
    Int motif_min_score = 35
    Float motif_min_fraction = 0.30
    Float mapping_min_concordance = 70
    Int mapping_min_length = 50
    String? mapping_biosample_name
    String? mapping_pbmm2_overrides
    Int add_memory_mb = 0

    Int nproc = 1
    Int max_nchunks = 100
    String log_level = "INFO"
    String? tmp_dir
  }

  call wf_prepare_input.prepare_input {
    input:
      dataset_xml = eid_ccs,
      dataset_filters = dataset_filters,
      filter_min_qv = filter_min_qv,
      nproc = 1,
      log_level = log_level
  }

  call ccs_kinetics_bystrandify {
    input:
      ccs_xml = prepare_input.reads_file,
      base_memory_mb = add_memory_mb,
      log_level = log_level
  }

  Int genome_length_mb_max = 20
  Float pbmm2_min_concordance = 98
  call pbmm2.pbmm2_align as mapping {
    input:
      unmapped = ccs_kinetics_bystrandify.subreads_xml,
      reference = eid_ref_dataset,
      genome_length_mb = genome_length_mb_max,
      min_concordance = pbmm2_min_concordance,
      base_memory_mb = add_memory_mb,
      nproc = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir
  }

  call pbreports.mapping_stats {
    input:
      mapped = mapping.mapped,
      all_reads = ccs_kinetics_bystrandify.subreads_xml,
      report_module = "mapping_stats_subreads",
      index_memory_gb = 4,
      base_memory_mb = add_memory_mb,
      nproc = 1,
      log_level = log_level
  }

  call wf_coverage_reports.coverage_reports {
    input:
      mapped = mapping.mapped,
      reference = eid_ref_dataset,
      index_memory_gb = 4,
      base_memory_mb = add_memory_mb,
      run_gc_coverage_plot = true,
      nproc = 1,
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
      genome_length_mb = genome_length_mb_max,
      base_memory_mb = add_memory_mb,
      max_nchunks = max_nchunks,
      nproc = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir
  }

  output {
    File mapped = mapping.mapped
    File? report_mapping_stats = mapping_stats.report
    File? coverage_gff = coverage_reports.coverage_gff
    File? report_coverage = coverage_reports.report_coverage
    File? basemods_gff = wf_basemods.basemods_gff
    File? basemods_csv = wf_basemods.basemods_csv
    File? report_modifications = wf_basemods.report_modifications
    File? bigwig_file = wf_basemods.bigwig_file
    File? motifs_csv = wf_basemods.motifs_csv
    File? motifs_gff = wf_basemods.motifs_gff
    File? report_motifs = wf_basemods.report_motifs
  }
}
