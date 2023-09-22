# Reference-Guided Long Amplicon Analysis (LAAGC)

version 1.0

import "tasks/pbcoretools.wdl" as pbcoretools
import "tasks/chunking.wdl" as chunking
import "pb_laa.wdl" as pb_laa

task laagc {
  input {
    File subreads
    File reference

    Boolean chimera_filter = true
    Boolean clustering = true
    Boolean phasing = true
    Boolean full_length = false
    Int min_length = 3000
    Int max_length = 0
    Int max_reads = 2000
    Float min_snr = 3.75
    Float min_read_score = 0.75
    Int min_barcode_score = 26
    Int max_phasing_reads = 500
    Int max_clustering_reads = 400
    Int ignore_ends = 0
    Float min_predicted_accuracy = 0.95
    Int trim_ends = 0
    Int min_split_reads = 20
    Float min_split_fraction = 0.1
    String? other_args
    Int min_guide_score = 50
    Int min_guide_span = 500

    Int nproc
    String log_level = "DEBUG"
  }
  command {
    laagc \
      --log-level ${log_level} \
      -j ${nproc} \
      ${true="" false="--noChimeraFilter" chimera_filter} \
      ${true="" false="--noClustering" clustering} \
      ${true="" false="--noPhasing" phasing} \
      ${true="--fullLength" false="" full_length} \
      ${"--maxReads " + max_reads} \
      ${"--minLength " + min_length} \
      ${"--maxLength " + max_length} \
      ${"--minSnr " + min_snr} \
      ${"--minReadScore " + min_read_score} \
      ${"--minBarcodeScore " + min_barcode_score} \
      ${"--maxPhasingReads " + max_phasing_reads} \
      ${"--maxClusteringReads " + max_clustering_reads} \
      ${"--minPredictedAccuracy " + min_predicted_accuracy} \
      ${"--minSplitReads " + min_split_reads} \
      ${"--minSplitFraction " + min_split_fraction} \
      ${"--ignoreEnds " + ignore_ends} \
      ${"--trimEnds " + trim_ends} \
      ${"--minGuideScore " + min_guide_score} \
      ${"--minGuideSpan " + min_guide_span} \
      --alarms alarms.json \
      ${other_args} \
      `readlink -f ${reference}` \
      `readlink -f ${subreads}` \
      --resultFile amplicon_analysis.fastq \
      --junkFile amplicon_analysis_chimeras_noise.fastq \
      --reportFile amplicon_analysis_summary.csv \
      --inputReportFile amplicon_analysis_input.csv
  }
  runtime {
    cpu: nproc
    memory: "12GB"
  }
  output {
    File consensus_fastq = "amplicon_analysis.fastq"
    File chimeras_fastq = "amplicon_analysis_chimeras_noise.fastq"
    File summary_csv = "amplicon_analysis_summary.csv"
    File input_csv = "amplicon_analysis_input.csv"
    File locus_csv = "amplicon_analysis_per_locus.csv"
  }
}

task laagc_input {
  input {
    File input_csv
    File locus_csv
    File barcoded_subreads

    Int nproc = 1
    String log_level = "INFO"
  }

  command {
    python3 \
      -m pbreports.report.laagc_input \
      --log-level ${log_level} \
      ${input_csv} \
      amplicon_analysis_input.report.json \
      ${locus_csv} \
      ${barcoded_subreads}
  }
  runtime {
    cpu: 1
    memory: "4GB"
  }
  output {
    File report = "amplicon_analysis_input.report.json"
  }
}

workflow pb_laagc {
  input {
    File eid_subread
    File eid_ref_dataset
    String dataset_filters = ""
    Int downsample_factor = 0
    Boolean laa_chimera_filter = true
    Boolean laa_clustering = true
    Boolean laa_phasing = true
    Boolean laa_full_length = false
    Int laa_min_length = 3000
    Int laa_max_length = 0
    Int laa_max_reads = 2000
    Float laa_min_snr = 2.5
    Float laa_min_read_score = 0.75
    Int laa_min_barcode_score = 26
    Int laa_max_phasing_reads = 500
    Int laa_max_clustering_reads = 400
    Int laa_ignore_ends = 0
    Int laagc_min_guide_score = 50
    Int laagc_min_guide_span = 500
    Float laa_min_predicted_accuracy = 0.95
    Int laa_trim_ends = 0
    Int laa_min_split_reads = 20
    Float laa_min_split_fraction = 0.1
    String? laa_override_args

    Int nproc = 8
    Int max_nchunks = 100
    Int? target_size
    String log_level = "DEBUG"
    String? tmp_dir
  }

  call pbcoretools.dataset_filter {
    input:
      dataset = eid_subread,
      filters = dataset_filters,
      downsample_factor = downsample_factor,
      log_level = log_level
  }

  call chunking.split_dataset {
    input:
      ds_in = dataset_filter.filtered,
      max_nchunks = max_nchunks,
      target_size = target_size,
      split_mode = "barcodes",
      log_level = log_level
  }

  scatter (chunk_subreads in split_dataset.chunks) {
    call laagc as laa {
      input:
        subreads = chunk_subreads,
        reference = eid_ref_dataset,
        chimera_filter = laa_chimera_filter,
        clustering = laa_clustering,
        phasing = laa_phasing,
        full_length = laa_full_length,
        min_length = laa_min_length,
        max_length = laa_max_length,
        max_reads = laa_max_reads,
        min_snr = laa_min_snr,
        min_read_score = laa_min_read_score,
        min_barcode_score = laa_min_barcode_score,
        max_phasing_reads = laa_max_phasing_reads,
        max_clustering_reads = laa_max_clustering_reads,
        min_predicted_accuracy = laa_min_predicted_accuracy,
        min_split_reads = laa_min_split_reads,
        min_split_fraction = laa_min_split_fraction,
        trim_ends = laa_trim_ends,
        ignore_ends = laa_ignore_ends,
        other_args = laa_override_args,
        min_guide_score = laagc_min_guide_score,
        min_guide_span = laagc_min_guide_span,
        log_level = log_level,
        nproc = nproc
    }
  }

  call chunking.gather_generic as gather_consensus_fastq {
    input:
      chunks = laa.consensus_fastq,
      output_file_name = "amplicon_analysis.fastq"
  }

  call chunking.gather_generic as gather_chimeras_fastq {
    input:
      chunks = laa.chimeras_fastq,
      output_file_name = "amplicon_analysis_chimeras_noise.fastq"
  }

  call chunking.gather_generic as gather_summary_csv {
    input:
      chunks = laa.summary_csv,
      output_file_name = "amplicon_consensus_summary.csv"
  }

  call chunking.gather_generic as gather_input_csv {
    input:
      chunks = laa.input_csv,
      output_file_name = "amplicon_consensus_input.csv"
  }

  call chunking.gather_generic as gather_locus_csv {
    input:
      chunks = laa.locus_csv,
      output_file_name = "amplicon_consensus_per_locus.csv"
  }

  call pb_laa.report_amplicon_analysis_consensus {
    input:
      summary_csv = gather_summary_csv.gathered,
      log_level = log_level
  }

  # this replaces pb_laa.report_amplicon_analysis_input
  call laagc_input {
    input:
      input_csv = gather_input_csv.gathered,
      locus_csv = gather_locus_csv.gathered,
      barcoded_subreads = dataset_filter.filtered,
      log_level = log_level
  }

  call pb_laa.split_laa_fastq {
    input:
      consensus = gather_consensus_fastq.gathered,
      chimeras = gather_chimeras_fastq.gathered,
      subreads = dataset_filter.filtered,
      nproc = nproc,
      log_level = log_level
  }

  call pb_laa.make_combined_laa_zip {
    input:
      consensus = gather_consensus_fastq.gathered,
      summary_csv = gather_summary_csv.gathered,
      subreads = dataset_filter.filtered,
      nproc = nproc,
      log_level = log_level
  }

  output {
    File consensus_fastq = gather_consensus_fastq.gathered
    File chimeras_fastq = gather_chimeras_fastq.gathered
    File summary_csv = gather_summary_csv.gathered
    File locus_csv = gather_locus_csv.gathered
    File report_consensus = report_amplicon_analysis_consensus.report
    File report_input = laagc_input.report
    File consensus_fastq_split = split_laa_fastq.consensus_fastq_zip
    File chimeras_fastq_split = split_laa_fastq.chimeras_fastq_zip
    File combined_zip = make_combined_laa_zip.combined_zip
  }
}
