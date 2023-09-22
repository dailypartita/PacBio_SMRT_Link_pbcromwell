# Resequencing using pbmm2 and gcpp

version 1.0

import "wf_mapping_subreads.wdl"
import "wf_consensus.wdl"
import "wf_coverage_reports.wdl"
import "tasks/pbcoretools.wdl" as pbcoretools
import "tasks/chunking.wdl" as chunking
import "tasks/pbreports.wdl" as pbreports

workflow pb_resequencing {
  Int TARGET_SIZE = 100000
  input {
    File eid_subread
    File eid_ref_dataset
    File? target_regions_bed
    String consensus_algorithm = "arrow"
    String dataset_filters = ""
    Float mapping_min_concordance = 70
    Int mapping_min_length = 50
    Int downsample_factor = 0
    Boolean extract_unmapped_bam = false
    String? mapping_biosample_name
    String? mapping_pbmm2_overrides

    Int nproc = 8
    String log_level = "INFO"
    String? tmp_dir
    Int max_nchunks = 100
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
      biosample_name = mapping_biosample_name,
      pbmm2_overrides = mapping_pbmm2_overrides,
      run_coverage = false,
      nproc = nproc,
      target_size = TARGET_SIZE,
      max_nchunks = max_nchunks,
      log_level = log_level,
      tmp_dir = tmp_dir
  }

  call wf_coverage_reports.coverage_reports {
    input:
      mapped = mapping.mapped,
      reference = eid_ref_dataset,
      target_regions_bed = target_regions_bed,
      index_memory_gb = mapping.index_memory_gb,
      nproc = nproc,
      log_level = log_level
  }

  if (extract_unmapped_bam) {
    call pbcoretools.extract_unmapped_reads {
      input:
        mapped = mapping.mapped,
        unmapped = dataset_filter.filtered,
        nproc = nproc,
        log_level = log_level
    }
  }

  call wf_consensus.consensus {
    input:
      alignments = mapping.mapped,
      reference = eid_ref_dataset,
      consensus_algorithm = consensus_algorithm,
      genome_length_mb = mapping.genome_length_mb,
      log_level = log_level,
      nproc = nproc,
      max_nchunks = max_nchunks
  }

  call pbreports.consensus_reports {
    input:
      reference = eid_ref_dataset,
      coverage_gff = coverage_reports.coverage_gff,
      variants_gff = consensus.variants_gff,
      log_level = log_level
  }

  output {
    File mapped = mapping.mapped
    File? report_mapping_stats = mapping.report_mapping_stats
    File coverage_gff = coverage_reports.coverage_gff
    File report_coverage = coverage_reports.report_coverage
    File? report_target_coverage = coverage_reports.report_target_coverage
    File consensus_fasta = consensus.consensus_fasta
    File consensus_fastq = consensus.consensus_fastq
    File variants_gff = consensus.variants_gff
    File variants_vcf = consensus.variants_vcf
    File consensus_gff = consensus_reports.consensus_gff
    File report_variants = consensus_reports.report_variants
    File report_top_variants = consensus_reports.report_top_variants
    File? mapped_bam_datastore = mapping.mapped_bam_datastore
    File? unmapped_bam = extract_unmapped_reads.unmapped_bam
  }
}
