version 1.0

task dataset_filter {
  input {
    File dataset
    String dataset_ext = ".subreadset.xml"
    String filters
    Int? downsample_factor
    Int? min_read_length
    Float? min_rq
    Int? min_qv

    Int nproc = 2
    String log_level = "INFO"
    Int memory_gb = 8
  }

  command {
    touch .NOARCHIVE
    set -vex
    python3 \
      -m pbcoretools.tasks.dataset_filter \
      --log-level ${log_level} \
      `readlink -f ${dataset}` \
      filtered${dataset_ext} \
      '${filters}' \
      ${"--min-rq " + min_rq} \
      ${"--min-qv " + min_qv} \
      ${"--min-read-length " + min_read_length} \
      ${"--downsample " + downsample_factor}
  }
  runtime {
    cpu: 1
    memory: "${memory_gb}GB"
  }
  output {
    File filtered = "filtered${dataset_ext}"
  }
}

task auto_ccs_outputs {
  input {
    File ccsxml
    String mode
    Float? min_rq
    Int? min_qv
    Int add_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }

  command {
    python3 \
      -m pbcoretools.tasks.auto_ccs_outputs \
      --log-level ${log_level} \
      ${"--min-rq " + min_rq} \
      ${"--min-qv " + min_qv} \
      ${mode} \
      ${ccsxml} \
      ${mode}.datastore.json
  }
  Int total_mem_mb = add_memory_mb + 4096
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File datastore = "${mode}.datastore.json"
  }
}

# More flexible consolidation task that will run only if combined BAM files
# are less than 10GB
task auto_consolidate_alignments {
  input {
    File mapped
    Boolean force_consolidate = false
    String dataset_ext = ".alignmentset.xml"
    String? tmp_dir

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    # We specify TMPDIR despite --noTmp b/c the XML will be written to TMPDIR anyway.
    # --noTmp refers only to BAM.
    set -vex

    ${"TMPDIR=" + tmp_dir} \
    python3 -m pbcoretools.tasks.auto_consolidate \
      --log-level ${log_level} \
      --noTmp \
      ${true="--force" false="" force_consolidate} \
      ${mapped} \
      --datastore mapped.datastore.json \
      mapped.bam
  }
  runtime {
    cpu: 2
    memory: "8GB"
  }
  output {
    # this avoids call-caching bulky files in Cromwell
    File? datastore = "mapped.datastore.json"
    # if we have a consolidated dataset we will use that as the primary
    # mapping output and delete the chunked BAMs
    File? dataset = "mapped${dataset_ext}"
    File? was_consolidated_flag = "dataset_was_consolidated.txt"
  }
}

task extract_unmapped_reads {
  input {
    File mapped
    File unmapped
    String bam_ext = ".subreads.bam"

    Int nproc = 1
    String log_level = "INFO"
    Int memory_gb = 4
  }
  command {
    bamsieve \
      --log-level ${log_level} \
      --blacklist ${mapped} \
      --subreads \
      ${unmapped} \
      unaligned${bam_ext}
  }
  runtime {
    cpu: 1
    memory: "${memory_gb}GB"
  }
  output {
    File unmapped_bam = "unaligned${bam_ext}"
  }
}

task pbvalidate {
  input {
    File data
    Boolean quick_mode = false

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    pbvalidate \
      --log-level ${log_level} \
      --alarms alarms.json \
      ${true="--quick" false="" quick_mode} \
      ${data}
  }
  runtime {
    cpu: 1
    memory: "4GB"
  }
  # pbvalidate doesn't actually write anything if it succeeds
  output {
    File? alarms = "alarms.json"
  }
}
