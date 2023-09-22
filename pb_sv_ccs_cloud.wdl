# Structural variants workflow that is independent of SMRT Link and various
# PacBio-specific file formats.  This is mostly copy-pasted from the customer
# facing pb_sv_ccs workflow, with several simplifications:
# 1. Reference FASTA and .fai index are passed explicitly to the workflow
#    and all dependent tasks
# 2. mapping outputs just BAM, no dataset
# 3. all tasks are defined inline (with some changes)
# 4. multi-sample support has been ripped out (we may add this back later)
# 5. scatter/gather is not used, all tasks run as single jobs

version 1.0

task pbmm2_align {
  input {
    File unmapped
    File reference_fasta
    File reference_fasta_fai

    String? preset_mode
    Float min_concordance = 70
    Int min_length = 50
    Boolean hq_mode = false
    Boolean zmw_mode = false
    String? pbmm2_overrides
    Boolean? median_filter = false
    Boolean? strip = false
    Boolean? split_by_sample = false

    Int nproc = 1
    String log_level = "DEBUG"
    String? tmp_dir
  }

  command {
    ln -s ${reference_fasta} reference.fasta
    ln -s ${reference_fasta_fai} reference.fasta.fai
    ${"TMPDIR=" + tmp_dir} \
    pbmm2 \
      align \
      reference.fasta \
      ${unmapped} \
      mapped.bam \
      --sort \
      --min-concordance-perc ${min_concordance} \
      --min-length ${min_length} \
      ${true="--median-filter" false="" median_filter} \
      ${true="--zmw" false="" zmw_mode} \
      ${true="--hqregion" false="" hq_mode} \
      ${"--preset " + preset_mode} \
      ${true="--strip" false="" strip} \
      ${pbmm2_overrides} \
      ${true="--split-by-sample" false="" split_by_sample} \
      -j ${nproc} \
      --log-level ${log_level}
  }
  runtime {
    cpu: nproc
  }
  output {
    File mapped = "mapped.bam"
  }
}

task pbsv_split_ref {
  input {
    File reference_fasta
    File reference_fasta_fai

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    ln -s ${reference_fasta} reference.fasta
    ln -s ${reference_fasta_fai} reference.fasta.fai
    python3 \
      -m pbsvtools.tasks.split_ref_to_chrs \
      --log-level ${log_level} \
      reference.fasta \
      pbsv.split_ref.csv
  }
  runtime {
    cpu: 1
  }
  output {
    File split_ref_csv = "pbsv.split_ref.csv"
  }
}

task pbsv_tandem_repeat_finder {
  input {
    File reference_fasta
    File reference_fasta_fai
    File split_ref_csv
    
    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    ln -s ${reference_fasta} reference.fasta
    ln -s ${reference_fasta_fai} reference.fasta.fai
    python3 \
      -m pbsvtools.tasks.tandem_repeat_finder \
      --log-level ${log_level} \
      reference.fasta \
      ${split_ref_csv} \
      pbsv.trf.bed
  }
  runtime {
    cpu: 1
  }
  output {
    File tandem_repeat_bed = "pbsv.trf.bed"
  }
}

task pbsv_discover {
  input {
    File mapped
    File tandem_repeat_bed

    Int nproc = 1
    String? log_level = "INFO"
  }
  command {
    pbsv discover -b ${tandem_repeat_bed} ${mapped} pbsv.svsig.gz
  }
  runtime {
    cpu: 1
  }
  output {
    File svsig_gz = "pbsv.svsig.gz"
  }
}

task pbsv_call {
  input {
    File reference_fasta
    File reference_fasta_fai
    File svsig_gz
    String chunk_length = "1M"
    Int min_sv_length = 20
    Int min_cnv_length = 1000
    Int min_percent_reads = 20
    Int min_reads_one_sample = 2
    Int min_reads_all_samples = 2
    Int min_reads_per_strand_all_samples = 1

    Int nproc
    String log_level = "INFO"
  }
  command {
    ln -s ${reference_fasta} reference.fasta
    ln -s ${reference_fasta_fai} reference.fasta.fai
    pbsv call \
      --log-level ${log_level} \
      -j ${nproc} \
      --chunk-length ${chunk_length} \
      --min-sv-length ${min_sv_length} \
      --min-cnv-length ${min_cnv_length} \
      --call-min-read-perc-one-sample ${min_percent_reads} \
      --call-min-reads-all-samples ${min_reads_all_samples} \
      --call-min-reads-one-sample ${min_reads_one_sample} \
      --call-min-reads-per-strand-all-samples ${min_reads_per_strand_all_samples} \
      reference.fasta \
      ${svsig_gz} \
      pbsv.vcf \
      structural_variants.lengths.json \
      structural_variants.summary.json
  }
  runtime {
    cpu: nproc
  }
  output {
    File variants = "pbsv.vcf"
    File lengths_report = "structural_variants.lengths.json"
    File summary_report = "structural_variants.summary.json"
  }
}

task structural_variants_2_report {
  input {
    File lengths_json
    File summary_json

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 \
      -m pbreports.report.structural_variants_2 \
      --log-level ${log_level} \
      ${summary_json} \
      ${lengths_json} \
      structural_variants.report.json
  }
  runtime {
    cpu: 1
  }
  output {
    File report = "structural_variants.report.json"
  }
}

workflow pb_sv_ccs_cloud {
  input {
    File reads
    File reference_fasta
    File reference_fasta_fai

    # Mapping parameters
    Int mapping_min_length = 50
    Float mapping_min_concordance = 70
    Boolean mapping_median_filter = false
    Boolean mapping_strip = true
    Boolean mapping_split_by_sample = false
    String? mapping_pbmm2_overrides

    # pbsv call parameters
    String chunk_length = "1M"
    Int min_sv_length = 20
    Int min_cnv_length = 1000
    # pbsv v2.2.2, CCS optimized parameters -A 1 -O 1 -S 0 -P 10
    Int min_percent_reads = 10                # -P 10
    Int min_reads_one_sample = 1              # -O 1
    Int min_reads_all_samples = 1             # -A 1
    Int min_reads_per_strand_all_samples = 0  # -S 0

    # Workflow resource configuration parameters
    Int nproc = 1
    String log_level = "INFO"
    Int max_nchunks = 100
  }

  call pbsv_split_ref {
    input:
      reference_fasta = reference_fasta,
      reference_fasta_fai = reference_fasta_fai,
      log_level = log_level,
      nproc = 1
  }

  call pbsv_tandem_repeat_finder {
    input:
      reference_fasta = reference_fasta,
      reference_fasta_fai = reference_fasta_fai,
      split_ref_csv = pbsv_split_ref.split_ref_csv,
      log_level = log_level,
      nproc = nproc
  }

  call pbmm2_align {
    input:
      unmapped = reads,
      reference_fasta = reference_fasta,
      reference_fasta_fai = reference_fasta_fai,
      min_concordance = mapping_min_concordance,
      median_filter = mapping_median_filter,
      strip = mapping_strip,
      split_by_sample = mapping_split_by_sample,
      preset_mode = "CCS",
      min_length = mapping_min_length,
      pbmm2_overrides = mapping_pbmm2_overrides,
      nproc = nproc,
      log_level = log_level
  }

  call pbsv_discover {
    input:
      mapped = pbmm2_align.mapped,
      tandem_repeat_bed = pbsv_tandem_repeat_finder.tandem_repeat_bed,
      log_level = log_level,
      nproc = nproc
  }

  call pbsv_call {  # svsig -> structural_variants.vcf, and json reports
    input:
      reference_fasta = reference_fasta,
      reference_fasta_fai = reference_fasta_fai,
      svsig_gz = pbsv_discover.svsig_gz,
      chunk_length = chunk_length,
      min_sv_length = min_sv_length,
      min_cnv_length = min_cnv_length,
      min_percent_reads = min_percent_reads,
      min_reads_one_sample = min_reads_one_sample,
      min_reads_all_samples = min_reads_all_samples,
      min_reads_per_strand_all_samples = min_reads_per_strand_all_samples,
      nproc = nproc,
      log_level = log_level
  }

  call structural_variants_2_report {
    input:
      lengths_json = pbsv_call.lengths_report,
      summary_json = pbsv_call.summary_report,
      nproc = 1,
      log_level = log_level
  }

  output {
    File variants = pbsv_call.variants
    File alignments_by_sample_datastore = pbmm2_align.mapped
    File report = structural_variants_2_report.report
  }
}
