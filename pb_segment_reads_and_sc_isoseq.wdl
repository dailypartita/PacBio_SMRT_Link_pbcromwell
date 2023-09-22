# MAS-Seq read segmentation and single-cell Iso-Seq workflow
# used internally by Run Design

version 1.0

import "pb_segment_reads.wdl"
import "pb_sc_isoseq.wdl"

workflow pb_segment_reads_and_sc_isoseq {
  input {
    File eid_ccs
    File? eid_barcode
    File eid_barcode_2
    File eid_ref_dataset
    File tenx_barcodes
    File? adapters_fasta

    # FIXME what is the default value?  is it an enumeration?
    String isoseq_design = "T-12U-16B"
    String output_prefix = "scisoseq"
    String cell_barcode_finding_method = "knee"
    Int cell_barcode_percentile_cutoff = 99
    Int add_memory_mb = 0

    String log_level = "INFO"
    Int nproc = 1
    # ignored
    Int max_nchunks = 1
    File? tmp_dir
  }

  call pb_segment_reads.pb_segment_reads {
    input:
      eid_ccs = eid_ccs,
      eid_barcode = eid_barcode,
      adapters_fasta = adapters_fasta,
      add_memory_mb = add_memory_mb,
      log_level = log_level,
      nproc = nproc
  }

  call pb_sc_isoseq.pb_sc_isoseq {
    input:
      eid_ccs = pb_segment_reads.segmented_reads,
      eid_barcode = eid_barcode_2,
      eid_ref_dataset = eid_ref_dataset,
      tenx_barcodes = tenx_barcodes,
      isoseq_design = isoseq_design,
      output_prefix = output_prefix,
      cell_barcode_finding_method = cell_barcode_finding_method,
      cell_barcode_percentile_cutoff = cell_barcode_percentile_cutoff,
      add_memory_mb = add_memory_mb,
      log_level = log_level,
      nproc = nproc
  }

  output {
    File segmented_reads = pb_segment_reads.segmented_reads
    Array[File] segmented_reads_bam = pb_segment_reads.segmented_reads_bam
    Array[File] non_passing_bam = pb_segment_reads.non_passing_bam
    File report_read_segmentation = pb_segment_reads.report_read_segmentation
    File report_sc_isoseq_read_statistics = pb_sc_isoseq.report_sc_isoseq_read_statistics
    File report_sc_isoseq_cell_statistics = pb_sc_isoseq.report_sc_isoseq_cell_statistics
    File report_sc_isoseq_transcript_statistics = pb_sc_isoseq.report_sc_isoseq_transcript_statistics
    File bcstats_report_tsv = pb_sc_isoseq.bcstats_report_tsv
    File dedup_fasta = pb_sc_isoseq.dedup_fasta
    File collapse_groups = pb_sc_isoseq.collapse_groups
    # pigeon outputs
    File mapped_transcripts_gff = pb_sc_isoseq.mapped_transcripts_gff
    File mapped_transcripts_filtered_gff = pb_sc_isoseq.mapped_transcripts_filtered_gff
    File classification_txt = pb_sc_isoseq.classification_txt
    File junctions_txt = pb_sc_isoseq.junctions_txt
    File classification_filtered_txt = pb_sc_isoseq.classification_filtered_txt
    File junctions_filtered_txt = pb_sc_isoseq.junctions_filtered_txt
    File seurat_info_tgz = pb_sc_isoseq.seurat_info_tgz
    # BAM file names
    String mapped_bam_fn = pb_sc_isoseq.mapped_bam_fn
    String mapped_bam_bai_fn = pb_sc_isoseq.mapped_bam_bai_fn
    String unmapped_bam_fn = pb_sc_isoseq.unmapped_bam_fn
  }
}
