# Single-cell Iso-Seq

version 1.0

task lima_sc_isoseq {
  input {
    File ccsxml
    File primers
    String prefix

    String log_level = "INFO"
    Int nproc = 8
  }
  String output_xml = "${prefix}.5p--3p.consensusreadset.xml"
  command {
    set -e
    # FIXME this is relatively inefficient, can we do this in lima instead?
    # this just emits alarms.json if sub-Q20 reads are present
    python3 -m pbcoretools.tasks.check_for_low_quality_reads \
      --log-level INFO \
      --ignore-segmented-reads \
      `readlink -f ${ccsxml}`

    lima \
      --log-level DEBUG \
      --log-file lima_isoseq.log \
      -j ${nproc} \
      --alarms alarms.json \
      --isoseq \
      --per-read \
      --ignore-xml-biosamples \
      `readlink -f ${ccsxml}` \
      `readlink -f ${primers}` \
      ${prefix}.bam

    # get the number of S-reads for future memory allocation
    python3 -m pbcoretools.tasks.memory.get_dataset_size \
      --skip-counts \
      ${output_xml}
  }
  runtime {
    cpu: nproc
    # lima itself takes very little memory; the rest is for the python task
    # which loads the .pbi, at 29 bytes per read, but we try to optimize this
    # by only loading raw CCS index, not segmented reads.
    memory: "4GB"
  }
  output {
    File xml = "${output_xml}"
    # this can be emitted by either task, so lima may overwrite the output
    # of the first task, but that usually means a fatal error anyway
    File? alarms = "alarms.json"
    File num_records_txt = "numrecords.txt"
  }
}

task isoseq_tag {
  input {
    File lima_reads
    String prefix
    String design
    Int number_of_sreads_m
    Int base_memory_mb

    Int nproc = 1
    String log_level = "INFO"
  }
  String output_bam = "${prefix}.5p--3p.tagged.bam"
  String output_xml = "${prefix}.5p--3p.tagged.consensusreadset.xml"
  command {
    isoseq3 tag \
      --log-level ${log_level} \
      --log-file isoseq_tag.log \
      -j ${nproc} \
      --alarms alarms.json \
      --design ${design} \
      `readlink -f ${lima_reads}` \
      ${output_xml}
  }
  # for 'dataset' (very generous, 4-cell job took 6.8GB)
  Int SCALE_FACTOR = 54
  Int total_mem_mb = base_memory_mb + (number_of_sreads_m * SCALE_FACTOR)
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File xml = "${output_xml}"
    File? alarms = "alarms.json"
  }
}

task isoseq_refine {
  input {
    File tagged_reads
    File primers
    String prefix
    Int base_memory_mb

    Int nproc = 1
    String log_level = "INFO"
  }
  String output_xml = "${prefix}.5p--3p.tagged.refined.consensusreadset.xml"
  String output_bam = "${prefix}.5p--3p.tagged.refined.bam"
  command {
    isoseq3 refine \
      --log-level ${log_level} \
      -j ${nproc} \
      --alarms alarms.json \
      --require-polya \
      `readlink -f ${tagged_reads}` \
      `readlink -f ${primers}` \
      ${output_xml}
  }
  runtime {
    cpu: nproc
    # even this is excessive
    memory: "${base_memory_mb}MB"
  }
  output {
    File xml = "${output_xml}"
    File report = "${prefix}.5p--3p.tagged.refined.filter_summary.report.json"
    File? alarms = "alarms.json"
  }
}

task isoseq_correct {
  input {
    File reads_in
    File tenx_barcodes
    String prefix
    String cell_barcode_finding_method = "knee"
    Int cell_barcode_percentile_cutoff = 99
    Int number_of_sreads_m
    Int base_memory_mb

    Int nproc = 1
    String log_level = "INFO"
  }
  String output_bam_1 = "${prefix}.5p--3p.tagged.refined.corrected.bam"
  String output_bam_2 = "${prefix}.5p--3p.tagged.refined.corrected.sorted.bam"
  String output_xml = "${prefix}.5p--3p.tagged.refined.corrected.sorted.transcriptset.xml"
  command {
    set -e

    isoseq3 correct \
      --log-level DEBUG \
      --verbose \
      --log-file isoseq_correct.log \
      -j ${nproc} \
      --alarms alarms.json \
      --method ${cell_barcode_finding_method} \
      --percentile ${cell_barcode_percentile_cutoff} \
      --barcodes ${tenx_barcodes} \
      `readlink -f ${reads_in}` \
      ${output_bam_1}

    samtools sort -t CB \
      -@ ${nproc} \
      -m ${memory_per_thread}M \
      ${output_bam_1} \
      -o ${output_bam_2}
  }
  # the correct step also apparently scales with nproc
  Int memory_per_thread = 512
  Int total_mem_mb = base_memory_mb + 8192 + memory_per_thread * nproc
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File bam = "${output_bam_2}"
    File report = "${prefix}.5p--3p.tagged.refined.corrected.report.json"
    File? alarms = "alarms.json"
  }
}

task isoseq_dedup {
  input {
    File sorted_transcripts
    String prefix
    Int number_of_sreads_m
    Int base_memory_mb

    Int nproc = 1
    String log_level = "INFO"
  }
  String output_prefix = "${prefix}.5p--3p.tagged.refined.corrected.sorted.dedup"
  command {
    set -e

    isoseq3 groupdedup \
      --log-level ${log_level} \
      --log-file isoseq_dedup.log \
      --alarms alarms.json \
      -j ${nproc} \
      `readlink -f ${sorted_transcripts}` \
      ${output_prefix}.bam

    echo -n "`pwd`/${output_prefix}.bam" > unmapped_bam.fofn
    echo -n "`pwd`/${output_prefix}.fasta" > dedup_fasta.fofn
  }
  # FIXME this is not a very good guess
  Int SCALE_FACTOR = 16 * nproc
  Int total_mem_mb = base_memory_mb + 1536 + (number_of_sreads_m * SCALE_FACTOR)
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    String fasta_fn = read_string("dedup_fasta.fofn")
    String unmapped_bam_fn = read_string("unmapped_bam.fofn")
    File? alarms = "alarms.json"
  }
}

task pbmm2_align {
  input {
    File transcripts
    File reference
    String prefix
    Int base_memory_mb

    Int nproc
    String log_level = "INFO"
  }
  String output_xml = "${prefix}.mapped.transcriptalignmentset.xml"
  command {
    set -e

    pbmm2 align \
      --log-level ${log_level} \
      --log-file pbmm2.log \
      -j ${nproc} \
      --alarms alarms.json \
      --preset ISOSEQ \
      --sort \
      --report-json mapping_stats.report.json \
      `readlink -f ${reference}` \
      `readlink -f ${transcripts}` \
      ${output_xml}

    echo -n "`pwd`/${prefix}.mapped.bam" > mapped_bam.fofn
    echo -n "`pwd`/${prefix}.mapped.bam.bai" > mapped_bam_bai.fofn
  }
  # XXX this is a human- and mouse-only application so I am hardcoding it
  Int genome_length_mb = 3300
  Int mem_scale_factor = 6
  Int total_mem_mb = base_memory_mb + 4096 + (genome_length_mb * mem_scale_factor) + (nproc * 256)
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File mapped = "${prefix}.mapped.transcriptalignmentset.xml"
    File report = "mapping_stats.report.json"
    # these are stored as String to avoid expensive call caching reads
    String mapped_bam_fn = read_string("mapped_bam.fofn")
    String mapped_bam_bai_fn = read_string("mapped_bam_bai.fofn")
    File? alarms = "alarms.json"
  }
}

task isoseq_collapse {
  input {
    File mapped_transcripts
    String prefix
    Int base_memory_mb

    Int nproc
    String log_level = "INFO"
  }
  String output_gff = "${prefix}.mapped_transcripts.collapse.gff"
  command {
    isoseq3 collapse \
      --log-level ${log_level} \
      --log-file isoseq_collapse.log \
      -j ${nproc} \
      --alarms alarms.json \
      `readlink -f ${mapped_transcripts}` \
      ${output_gff}
  }
  Int mem_total_mb = base_memory_mb + 16384 + nproc * 1024
  runtime {
    cpu: nproc
    memory: "${mem_total_mb}MB"
  }
  output {
    File gff = "${output_gff}"
    File abundance_txt = "${prefix}.mapped_transcripts.collapse.abundance.txt"
    File groups_txt = "${prefix}.mapped_transcripts.collapse.group.txt"
    File? alarms = "alarms.json"
  }
}

task pigeon {
  input {
    File gff_in
    File abundance
    File reference
    # FIXME this should be XML instead
    String dedup_fasta
    File groups
    File? poly_a
    File? cage
    String prefix
    Int number_of_sreads_m
    Int base_memory_mb

    String log_level = "INFO"
    Int nproc = 1
  }
  String sorted_gff = "${prefix}_transcripts.sorted.gff"
  String classify_out_1 = "${prefix}_classification.txt"
  String classify_out_2 = "${prefix}_classification.filtered_lite_classification.txt"
  String report_txt_out = "${prefix}_saturation.txt"
  command {
    set -e
    # FIXME pigeon should handle relative paths properly
    dataset \
      --log-level ${log_level} \
      --log-file dataset.log \
      absolutize \
      --output-xml absolutized.referenceset.xml \
      `readlink -f ${reference}`

    pigeon sort \
      --log-level ${log_level} \
      --log-file pigeon-sort.log \
      -o ${sorted_gff} \
      ${gff_in}

    pigeon classify \
      --log-level ${log_level} \
      --log-file pigeon-classify.log \
      -j ${nproc} \
      ${"--poly-a " + poly_a} \
      ${"--cage-peak " + cage} \
      --out-dir . \
      --out-prefix ${prefix} \
      --flnc ${abundance} \
      --ref absolutized.referenceset.xml \
      ${sorted_gff}

    pigeon filter \
      --log-level ${log_level} \
      --log-file pigeon-filter.log \
      --isoforms ${sorted_gff} \
      ${classify_out_1}

    pigeon report \
      --log-level ${log_level} \
      --log-file pigeon-report.log \
      ${classify_out_2} \
      ${report_txt_out}

    pigeon make-seurat \
      --log-level ${log_level} \
      --log-file pigeon-make-seurat.log \
      --num-threads ${nproc} \
      --annotations absolutized.referenceset.xml \
      --dedup ${dedup_fasta} \
      --group ${groups} \
      --out-dir . \
      --out-prefix ${prefix} \
      ${classify_out_2}

    tar cvzf ${prefix}.seurat_info.tar.gz \
      isoforms_seurat genes_seurat \
      ${prefix}.annotated.info.csv ${prefix}.info.csv
    echo "`pwd`/${prefix}.seurat_info.tar.gz" > seurat_info.fofn
  }
  # this is all `make-seurat`
  Int SCALE_FACTOR = 450
  # Cromwell integers are 32-bit, so we need to scale down to avoid overflow
  Int total_mem_mb = base_memory_mb + 4096 + (number_of_sreads_m * SCALE_FACTOR)
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File classification_txt = "${classify_out_1}"
    File junctions_txt = "${prefix}_junctions.txt"
    File classification_filtered_txt = "${classify_out_2}"
    File junctions_filtered_txt = "${prefix}_classification.filtered_lite_junctions.txt"
    File filter_reasons_txt = "${prefix}_classification.filtered_lite_reasons.txt"
    File report_classify = "${prefix}.report.json"
    File report_filter = "${prefix}_classification.filtered.report.json"
    File saturation_txt = "${report_txt_out}"
    File sorted_transcripts_gff = "${sorted_gff}"
    File filtered_gff = "${prefix}_transcripts.sorted.filtered_lite.gff"
    String seurat_info_tgz_fn = read_string("seurat_info.fofn")
  }
}

task isoseq_bcstats {
  input {
    File corrected_transcripts
    String cell_barcode_finding_method = "knee"
    Int cell_barcode_percentile_cutoff = 99
    Int number_of_sreads_m
    Int base_memory_mb

    String log_level = "INFO"
    Int nproc = 1
  }
  command {
    set -e
    isoseq3 bcstats \
      --log-level ${log_level} \
      --log-file isoseq_bcstats.log \
      -j ${nproc} \
      --alarms alarms.json \
      --method ${cell_barcode_finding_method} \
      --percentile ${cell_barcode_percentile_cutoff} \
      -o bcstats_report.tsv \
      --json bcstats.report.json \
      `readlink -f ${corrected_transcripts}`
    gzip bcstats_report.tsv
  }
  # FIXME this is too much guesswork
  Int SCALE_FACTOR = 8 * nproc
  Int total_mem_mb = base_memory_mb + 4096 + (number_of_sreads_m * SCALE_FACTOR)
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File report = "bcstats.report.json"
    File report_tsv = "bcstats_report.tsv.gz"
    File? alarms = "alarms.json"
  }
}

task pbreports_sc_isoseq {
  input {
    File ccsxml
    File report_isoseq_refine
    File report_isoseq_correct
    File report_isoseq_bcstats
    File report_mapping_stats
    File report_pigeon_classify
    File report_pigeon_filter
    File pigeon_classification
    File pigeon_classification_filtered
    File pigeon_saturation_txt
    File bcstats_report_tsv
    File dedup_bam
    Int base_memory_mb

    String log_level = "INFO"
    Int nproc = 1
  }
  command {
    set -e
    # most of the memory overhead is here (because of input CCS)
    python3 -m pbreports.report.sc_isoseq_read_statistics \
      --log-level ${log_level} \
      --log-file pbreports_read_statistics.log \
      -o read_statistics.report.json \
      `readlink -f ${ccsxml}` \
      ${report_isoseq_refine} \
      ${report_isoseq_correct} \
      `readlink -f ${dedup_bam}`

    python3 -m pbreports.report.sc_isoseq_cell_statistics \
      --log-level ${log_level} \
      --log-file pbreports_cell_statistics.log \
      -o cell_statistics.report.json \
      ${report_isoseq_bcstats} \
      ${bcstats_report_tsv}

    python3 -m pbreports.report.sc_isoseq_transcript_statistics \
      --log-level ${log_level} \
      --log-file pbreports_transcript_statistics.log \
      -o transcript_statistics.report.json \
      ${report_mapping_stats} \
      ${report_pigeon_classify} \
      ${report_pigeon_filter} \
      ${pigeon_classification} \
      ${pigeon_classification_filtered} \
      ${pigeon_saturation_txt}
  }
  runtime {
    cpu: 1
    # this is a massive overestimate, just to be safe
    memory: "${base_memory_mb}MB"
  }
  output {
    File report_read_statistics = "read_statistics.report.json"
    File report_cell_statistics = "cell_statistics.report.json"
    File report_transcript_statistics = "transcript_statistics.report.json"
    Array[File?] plot_pngs = glob("*.png")
  }
}

workflow pb_sc_isoseq {
  input {
    File eid_ccs
    File eid_barcode
    File eid_ref_dataset
    # The default value for this is specified in task_options.json.
    # The '${SMRT_DATA}' component will be replaced with the actual path to
    # bundled data in the SL tarball by the pipeline templates service.
    # CLI users will need to specify this manually.
    File tenx_barcodes
    # these will normally be bundled with the reference XML
    #File? pigeon_poly_a
    #File? pigeon_cage

    String output_prefix = "scisoseq"
    # FIXME what is the default value?  is it an enumeration?
    String isoseq_design = "T-12U-16B"
    String cell_barcode_finding_method = "knee"
    Int cell_barcode_percentile_cutoff = 99
    Int add_memory_mb = 0

    String log_level = "INFO"
    Int nproc = 1
    Int max_nchunks = 0
    String? tmp_dir
  }
  Int BASE_MEMORY_APP = 512
  Int base_memory_mb = BASE_MEMORY_APP + add_memory_mb

  call lima_sc_isoseq {
    input:
      ccsxml = eid_ccs,
      primers = eid_barcode,
      prefix = output_prefix,
      nproc = nproc,
      log_level = log_level
  }
  Int number_of_sreads = read_int(lima_sc_isoseq.num_records_txt)
  Int number_of_sreads_m = number_of_sreads / 1048576

  call isoseq_tag {
    input:
      lima_reads = lima_sc_isoseq.xml,
      prefix = output_prefix,
      design = isoseq_design,
      number_of_sreads_m = number_of_sreads_m,
      base_memory_mb = base_memory_mb,
      nproc = nproc,
      log_level = log_level
  }

  call isoseq_refine {
    input:
      tagged_reads = isoseq_tag.xml,
      primers = eid_barcode,
      prefix = output_prefix,
      base_memory_mb = base_memory_mb,
      nproc = nproc,
      log_level = log_level
  }

  call isoseq_correct {
    input:
      reads_in = isoseq_refine.xml,
      tenx_barcodes = tenx_barcodes,
      prefix = output_prefix,
      cell_barcode_finding_method = cell_barcode_finding_method,
      cell_barcode_percentile_cutoff = cell_barcode_percentile_cutoff,
      number_of_sreads_m = number_of_sreads_m,
      base_memory_mb = base_memory_mb,
      nproc = nproc,
      log_level = log_level
  }

  call isoseq_dedup {
    input:
      sorted_transcripts = isoseq_correct.bam,
      number_of_sreads_m = number_of_sreads_m,
      base_memory_mb = base_memory_mb,
      prefix = output_prefix,
      nproc = nproc,
      log_level = log_level
  }

  call pbmm2_align {
    input:
      transcripts = isoseq_dedup.unmapped_bam_fn,
      reference = eid_ref_dataset,
      base_memory_mb = base_memory_mb,
      prefix = output_prefix,
      nproc = nproc,
      log_level = log_level
  }

  call isoseq_collapse {
    input:
      mapped_transcripts = pbmm2_align.mapped,
      prefix = output_prefix,
      base_memory_mb = base_memory_mb,
      nproc = nproc,
      log_level = log_level
  }

  call isoseq_bcstats {
    input:
      corrected_transcripts = isoseq_correct.bam,
      number_of_sreads_m = number_of_sreads_m,
      cell_barcode_finding_method = cell_barcode_finding_method,
      cell_barcode_percentile_cutoff = cell_barcode_percentile_cutoff,
      base_memory_mb = base_memory_mb,
      nproc = nproc,
      log_level = log_level
  }

  call pigeon {
    input:
      gff_in = isoseq_collapse.gff,
      abundance = isoseq_collapse.abundance_txt,
      reference = eid_ref_dataset,
      dedup_fasta = isoseq_dedup.fasta_fn,
      groups = isoseq_collapse.groups_txt,
      prefix = output_prefix,
      number_of_sreads_m = number_of_sreads_m,
      base_memory_mb = base_memory_mb,
      nproc = nproc,
      log_level = log_level
  }

  call pbreports_sc_isoseq {
    input:
      log_level = log_level,
      ccsxml = eid_ccs,
      report_mapping_stats = pbmm2_align.report,
      report_isoseq_refine = isoseq_refine.report,
      report_isoseq_correct = isoseq_correct.report,
      report_isoseq_bcstats = isoseq_bcstats.report,
      report_pigeon_classify = pigeon.report_classify,
      report_pigeon_filter = pigeon.report_filter,
      pigeon_classification = pigeon.classification_txt,
      pigeon_classification_filtered = pigeon.classification_filtered_txt,
      pigeon_saturation_txt = pigeon.saturation_txt,
      bcstats_report_tsv = isoseq_bcstats.report_tsv,
      base_memory_mb = base_memory_mb,
      dedup_bam = isoseq_dedup.unmapped_bam_fn,
      log_level = log_level
  }

  output {
    File report_sc_isoseq_read_statistics = pbreports_sc_isoseq.report_read_statistics
    File report_sc_isoseq_cell_statistics = pbreports_sc_isoseq.report_cell_statistics
    File report_sc_isoseq_transcript_statistics = pbreports_sc_isoseq.report_transcript_statistics
    File bcstats_report_tsv = isoseq_bcstats.report_tsv
    File dedup_fasta = isoseq_dedup.fasta_fn
    File collapse_groups = isoseq_collapse.groups_txt
    # pigeon outputs
    File mapped_transcripts_gff = pigeon.sorted_transcripts_gff
    File mapped_transcripts_filtered_gff = pigeon.filtered_gff
    File classification_txt = pigeon.classification_txt
    File junctions_txt = pigeon.junctions_txt
    File classification_filtered_txt = pigeon.classification_filtered_txt
    File junctions_filtered_txt = pigeon.junctions_filtered_txt
    String seurat_info_tgz = pigeon.seurat_info_tgz_fn
    # BAM file names
    String mapped_bam_fn = pbmm2_align.mapped_bam_fn
    String mapped_bam_bai_fn = pbmm2_align.mapped_bam_bai_fn
    String unmapped_bam_fn = isoseq_dedup.unmapped_bam_fn
  }
}
