version 1.0

# Alignment (should be chunked by ZMWs)
task pbmm2_align {
  input {
    File unmapped
    File? reference
    String aln_ext = ".alignmentset.xml"
    String? preset_mode
    Float min_concordance = 70
    Int min_length = 50
    Boolean hq_mode = false
    Boolean zmw_mode = false
    String? biosample_name
    String? pbmm2_overrides
    Boolean? median_filter = false
    Boolean? strip = false
    Boolean? split_by_sample = false

    Int nproc = 1
    String log_level = "DEBUG"
    String? tmp_dir

    #FIXME these are all bad approximations
    Int genome_length_mb = 3300
    Int mem_scale_factor = 4
    Int base_memory_mb = 0
    Int min_memory_mb = 4096
  }
  command {
    set -vex
    ${"TMPDIR=" + tmp_dir} \
    pbmm2 \
      align \
      `readlink -f ${reference}` \
      `readlink -f ${unmapped}` \
      mapped${aln_ext} \
      --sort \
      --min-gap-comp-id-perc ${min_concordance} \
      --min-length ${min_length} \
      --sample "${biosample_name}" \
      ${true="--median-filter" false="" median_filter} \
      ${true="--zmw" false="" zmw_mode} \
      ${true="--hqregion" false="" hq_mode} \
      ${"--preset " + preset_mode} \
      ${true="--strip" false="" strip} \
      ${pbmm2_overrides} \
      ${true="--split-by-sample" false="" split_by_sample} \
      -j ${nproc} \
      --log-level ${log_level} \
      --log-file pbmm2.log \
      --alarms alarms.json
    # special handling for split_by_sample=true
    (ls mapped${aln_ext}) || ( \
    cp mapped.json mapped.datastore.json && \
    python3 \
      -m pbcoretools.tasks.dataset_to_datastore \
      mapped.datastore.json \
      mapped${aln_ext}
    )
    if [ -f "mapped.bam" ]; then
      echo "`pwd`/mapped.bam" > mapped_bam.fofn
      echo "`pwd`/mapped.bam.bai" > mapped_bam_bai.fofn
    fi
  }
  #FIXME this is poorly understood or validated
  # we do not use precalculated indices in our references, so we need to index
  # on the fly which consumes more memory
  Int NPROC_SCALE_FACTOR = 256
  Int total_mem_mb = min_memory_mb + base_memory_mb + (genome_length_mb * mem_scale_factor) + (nproc * NPROC_SCALE_FACTOR)
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File mapped = "mapped${aln_ext}"
    File? bam_datastore = "mapped.datastore.json"
    File? mapped_bam_fofn = "mapped_bam.fofn"
    File? mapped_bam_bai_fofn = "mapped_bam_bai.fofn"
  }
}
