# Common workflow for pb_sv_ccs and pb_sv_clr

version 1.0

import "tasks/pbcoretools.wdl" as pbcoretools
import "tasks/chunking.wdl" as chunking

task pbsv_scatter_tandem_repeat_finder {
  input {
    File reference
    Int max_nchunks
    Int base_memory_mb = 0

    Int? nproc = 1
    String log_level = "INFO"
  }
  command {
    set -e
    touch .NOARCHIVE
    python3 \
      -m pbsvtools.tasks.split_ref_to_chrs \
      --log-level ${log_level} \
      --log-file split_ref_to_chrs.log \
      `readlink -f ${reference}` \
      pbsv.split_ref.csv
    python3 \
      -m pbsvtools.tasks.scatter_tandem_repeat_finder2 \
      --log-level ${log_level} \
      --log-file scatter_tandem_repeat_finder2.log \
      pbsv.split_ref.csv \
      . \
      ${max_nchunks}
  }
  Int total_mem_mb = 1024 + base_memory_mb
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File split_ref_csv = "pbsv.split_ref.csv"
    Array[File] chunks = glob("scatter_trf_chunk*.csv")
    Int nchunks = length(chunks)
  }
}

task pbsv_tandem_repeat_finder {
  input {
    File reference
    File split_ref_csv
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 \
      -m pbsvtools.tasks.tandem_repeat_finder \
      --log-level ${log_level} \
      --log-file tandem_repeat_finder.log \
      `readlink -f ${reference}` \
      ${split_ref_csv} \
      pbsv.trf.bed
  }
  Int total_mem_mb = 2048 + base_memory_mb
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File tandem_repeat_bed = "pbsv.trf.bed"
  }
}

task pbsv_scatter_alignments {
  input {
    File mapped
    Int max_nchunks
    Int base_memory_mb = 0

    Int? nproc = 1
    String? log_level = "INFO"
  }
  command {
    set -e
    python3 \
      -m pbcoretools.tasks.dataset_to_datastore \
      --log-level ${log_level} \
      --log-file dataset_to_datastore.log \
      `readlink -f ${mapped}` \
      mapped.datastore.json
    python3 \
      -m pbsvtools.tasks.scatter_align_datastore \
      --log-level ${log_level} \
      --log-file scatter_align_datastore.log \
      mapped.datastore.json \
      . \
      ${max_nchunks} \
  }
  Int total_mem_mb = 8192 + base_memory_mb
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    Array[File] chunks = glob("chunk*.datastore.json")
    File mapped_datastore = "mapped.datastore.json"
    Int nchunks = length(chunks)
  }
}

# calls 'pbsv discover'
task pbsv_align_json_to_svsig {
  input {
    File align_datastore_json
    File tandem_repeat_bed
    Boolean hifi_mode = true
    Int base_memory_mb = 0

    Int nproc = 1
    String? log_level = "INFO"
  }
  command {
    python3 \
      -m pbsvtools.tasks.align_json_to_svsig \
      --log-level ${log_level} \
      --log-file align_json_to_svsig.log \
      ${true="--hifi" false="" hifi_mode} \
      ${align_datastore_json} \
      ${tandem_repeat_bed} \
      pbsv.svsig.gz.fofn
  }
  Int total_mem_mb = 8192 + base_memory_mb
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File svsig_gz_fofn = "pbsv.svsig.gz.fofn"
  }
}

task pbsv_call {
  input {
    File reference
    File svsig_gz_fofn
    String chunk_length = "1M"
    Int min_sv_length = 20
    Int min_percent_reads = 20
    Int min_reads_one_sample = 2
    Int min_reads_all_samples = 2
    Int min_reads_per_strand_all_samples = 1
    Boolean hifi_mode = true
    String? other_args
    Int base_memory_mb = 0

    Int nproc
    String log_level = "INFO"
  }
  String variant_types = "DEL,INS,INV,DUP,BND"
  command {
    set -e
    pbsv call \
      --log-level ${log_level} \
      --log-file pbsv-call.log \
      -j ${nproc} \
      --alarms alarms.json \
      ${true="--hifi" false="" hifi_mode} \
      --types ${variant_types} \
      --chunk-length ${chunk_length} \
      --min-sv-length ${min_sv_length} \
      --call-min-read-perc-one-sample ${min_percent_reads} \
      --call-min-reads-all-samples ${min_reads_all_samples} \
      --call-min-reads-one-sample ${min_reads_one_sample} \
      --call-min-reads-per-strand-all-samples ${min_reads_per_strand_all_samples} \
      ${other_args} \
      `readlink -f ${reference}` \
      ${svsig_gz_fofn} \
      pbsv.vcf \
      structural_variants.lengths.json \
      structural_variants.summary.json
    bgzip pbsv.vcf
    tabix pbsv.vcf.gz
  }
  # best guess for human scale
  Int total_mem_mb = (1024 * 64) + base_memory_mb
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File variants = "pbsv.vcf.gz"
    File variants_index = "pbsv.vcf.gz.tbi"
    File lengths_report = "structural_variants.lengths.json"
    File summary_report = "structural_variants.summary.json"
  }
}


task pbsv_scatter_align_datastore_by_sample {
  input {
    File align_datastore_json
    Int max_nchunks
    Int base_memory_mb = 0

    Int? nproc
    String? log_level = "INFO"
  }
  command {
    touch .NOARCHIVE
    python3 \
      -m pbsvtools.tasks.scatter_align_datastore_by_sample \
      --log-level ${log_level} \
      --log-file scatter_align_datastore_by_sample.log \
      ${align_datastore_json} \
      . \
      ${max_nchunks}
  }
  Int total_mem_mb = 8192 + base_memory_mb
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    Array[File] chunks = glob("chunk_by_sample*.datastore.json")
    Int nchunks = length(chunks)
  }
}


task pbsv_merge_alignments_by_sample {
  input {
    File align_datastore_json
    Int i_chunk = 0
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    export PB_WORKFLOW_CHUNK_ID=${i_chunk}
    python3 \
      -m pbsvtools.tasks.merge_alignments_by_sample \
      --nproc ${nproc} \
      --log-level ${log_level} \
      --log-file merge_alignments_by_sample.log \
      ${align_datastore_json} \
      alignments_by_sample.datastore.json
  }
  Int total_mem_mb = 1024 * (nproc + 1) + base_memory_mb
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File alignments_by_sample_datastore = "alignments_by_sample.datastore.json"
  }
}

task report_structural_variants_2 {
  input {
    File lengths_json
    File summary_json
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 \
      -m pbreports.report.structural_variants_2 \
      --log-level ${log_level} \
      --log-file pbreports.log \
      ${summary_json} \
      ${lengths_json} \
      structural_variants.report.json
  }
  Int total_mem_mb = 1024 + base_memory_mb
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File report = "structural_variants.report.json"
    Array[File?] plot_pngs = glob("*.png")
  }
}

# Support scatter + gather for mapping, tandem repeat finder, pbsv discover, bam merge

workflow pbsv_intl {
  input {
    # note these may have relative paths and so always need to be resolved to
    # the true path with `readlink -f`
    File mapped
    File reference

    # pbsv call parameters
    String chunk_length = "1M"
    Int min_sv_length = 20
    Int min_percent_reads = 20
    Int min_reads_one_sample = 2
    Int min_reads_all_samples = 2
    Int min_reads_per_strand_all_samples = 1
    String? pbsv_override_args
    Boolean hifi_mode
    Int base_memory_mb = 0

    # Workflow resource configuration parameters
    Int nproc = 1
    String log_level = "INFO"
    Int max_nchunks = 24
  }

  call pbsv_scatter_tandem_repeat_finder {
    input:
      reference = reference,
      max_nchunks = max_nchunks,
      base_memory_mb = base_memory_mb,
      log_level = log_level,
      nproc = 1
  }

  scatter (chunk in pbsv_scatter_tandem_repeat_finder.chunks) {
    call pbsv_tandem_repeat_finder {
      input:
        reference = reference,
        split_ref_csv = chunk,
        base_memory_mb = base_memory_mb,
        log_level = log_level,
        nproc = nproc
    }
  }

  call chunking.gather_generic as gather_bed {
    input:
      chunks = pbsv_tandem_repeat_finder.tandem_repeat_bed,
      output_file_name = "gathered.pbsv.trf.bed",
      log_level = log_level
  }

  call pbsv_scatter_alignments {
    input:
      mapped = mapped,
      max_nchunks = max_nchunks,
      base_memory_mb = base_memory_mb,
      log_level = log_level
  }

  scatter (chunk in pbsv_scatter_alignments.chunks) {
    call pbsv_align_json_to_svsig {  # alignments -> svsig
      input:
        align_datastore_json = chunk,
        tandem_repeat_bed = gather_bed.gathered,
        hifi_mode = hifi_mode,
        base_memory_mb = base_memory_mb,
        log_level = log_level,
        nproc = nproc
    }
  }

  call chunking.gather_generic as pbsv_gather_svsig {
    input:
      chunks = pbsv_align_json_to_svsig.svsig_gz_fofn,
      output_file_name = "gathered.pbsv.svsig.fofn",
      log_level = log_level
  }

  call pbsv_call {  # svsig -> structural_variants.vcf, and json reports
    input:
      reference = reference,
      svsig_gz_fofn = pbsv_gather_svsig.gathered,
      chunk_length = chunk_length,
      min_sv_length = min_sv_length,
      min_percent_reads = min_percent_reads,
      min_reads_one_sample = min_reads_one_sample,
      min_reads_all_samples = min_reads_all_samples,
      min_reads_per_strand_all_samples = min_reads_per_strand_all_samples,
      other_args = pbsv_override_args,
      hifi_mode = hifi_mode,
      base_memory_mb = base_memory_mb,
      nproc = nproc,
      log_level = log_level
  }

  call pbsv_scatter_align_datastore_by_sample as scatter_datastore {
    # one chunk for each biosample
    input:
      align_datastore_json = pbsv_scatter_alignments.mapped_datastore,
      max_nchunks = max_nchunks,
      base_memory_mb = base_memory_mb,
      log_level = log_level
  }

  scatter (i_chunk in range(length(scatter_datastore.chunks))) {
    call pbsv_merge_alignments_by_sample {
      # output one alignments.bam for each sample
      input:
        align_datastore_json = scatter_datastore.chunks[i_chunk],
        i_chunk = i_chunk,
        base_memory_mb = base_memory_mb,
        nproc = nproc,
        log_level = log_level
    }
  }

  call chunking.gather_generic as gather_alignments {
    input:
      chunks = pbsv_merge_alignments_by_sample.alignments_by_sample_datastore,
      output_file_name = "alignments_by_sample.datastore.json",
      log_level = log_level
  }

  call report_structural_variants_2 {
    input:
      lengths_json = pbsv_call.lengths_report,
      summary_json = pbsv_call.summary_report,
      base_memory_mb = base_memory_mb,
      nproc = 1,
      log_level = log_level
  }

  output {
    File variants = pbsv_call.variants
    File variants_index = pbsv_call.variants_index
    File alignments_by_sample_datastore = gather_alignments.gathered
    File report = report_structural_variants_2.report
    Array[File?] report_plot_pngs = report_structural_variants_2.plot_pngs
  }
}
