# Standalone genomic consensus starting from AlignmentSet

version 1.0

import "wf_consensus.wdl"
import "tasks/memory.wdl" as memory
import "tasks/pbreports.wdl" as pbreports

workflow pb_consensus {
  input {
    File eid_alignment
    File eid_ref_dataset
    # note that this needs to be handled separately from the other entry points
    File? coverage_gff

    String consensus_algorithm = "arrow"
    Int min_coverage = 5
    Int max_coverage = 100
    Int min_confidence = 40
    String? gcpp_advanced_args
    Boolean report_effective_coverage = false

    Int nproc = 8
    String log_level = "INFO"
    String? tmp_dir
    Int max_nchunks = 100
  }

  call memory.get_input_sizes {
    input:
      bam_dataset = eid_alignment,
      fasta_dataset = eid_ref_dataset,
      log_level = log_level
  }

  call wf_consensus.consensus {
    input:
      alignments = eid_alignment,
      reference = eid_ref_dataset,
      consensus_algorithm = consensus_algorithm,
      min_coverage = min_coverage,
      max_coverage = max_coverage,
      min_confidence = min_confidence,
      gcpp_advanced_args = gcpp_advanced_args,
      report_effective_coverage = report_effective_coverage,
      genome_length_mb = get_input_sizes.genome_length_mb,
      log_level = log_level,
      nproc = nproc,
      max_nchunks = max_nchunks
  }

  if (defined(coverage_gff)) {
    File coverage_gff_ = select_first([coverage_gff])

    if (coverage_gff != "") {
      call pbreports.consensus_reports {
        input:
          reference = eid_ref_dataset,
          coverage_gff = coverage_gff_,
          variants_gff = consensus.variants_gff
      }
    }
  }

  output {
    File consensus_fasta = consensus.consensus_fasta
    File consensus_fastq = consensus.consensus_fastq
    File variants_gff = consensus.variants_gff
    File variants_vcf = consensus.variants_vcf
    File? consensus_gff = consensus_reports.consensus_gff
    File? report_variants = consensus_reports.report_variants
    File? report_top_variants = consensus_reports.report_top_variants
  }
}
