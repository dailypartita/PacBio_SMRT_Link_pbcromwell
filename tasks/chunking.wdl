version 1.0

# Split any BAM-based dataset by holeNumber range (real for Subreads or CCS,
# synthetic for Transcripts)
task split_dataset {
  input {
    File ds_in
    String split_mode = "zmws"
    String? split_args

    Int? max_nchunks
    Int? target_size

    Int nproc = 1
    String log_level = "DEBUG"
    Int memory_gb = 4
    Int base_memory_mb = 0
  }
  command {
    touch .NOARCHIVE
    dataset \
      --log-level ${log_level} \
      split --${split_mode} \
      --simple-chunk-ids \
      ${split_args} \
      ${"--maxChunks " + max_nchunks} \
      ${"--targetSize " + target_size} \
      --outdir . \
      `readlink -f ${ds_in}`
  }
  # splitting the dataset adds overhead, but there is also a big spike in
  # memory when the index is loaded
  Int total_mem_mb = base_memory_mb + 1024 * (1 + memory_gb * 3)
  runtime {
    # single-processor but may have a large memory footprint for subread-based
    # applications
    # FIXME we should just use one processor and pass memory_gb to the JMS
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    # Note that unlike pbsmrtpipe tooling, these currently retain the basename
    # prefix of the input dataset.  We can override this with the new --prefix
    # option if desired.
    Array[File] chunks = glob("*set.xml")
    Int nchunks = length(chunks)
  }
}

task split_alignments {
  input {
    File mapped
    Int? max_nchunks
    Int? target_size

    Int nproc = 1
    String log_level = "DEBUG"
    Int memory_gb = 4
  }
  # splitting the dataset adds overhead
  Int memory_gb_final = memory_gb * 2
  command {
    touch .NOARCHIVE
    dataset \
      --log-level ${log_level} \
      split \
      ${"--maxChunks " + max_nchunks} \
      ${"--targetSize " + target_size} \
      --contigs --breakContigs \
      --outdir . \
      ${mapped}
  }
  runtime {
    cpu: 1
    memory: "${memory_gb_final}GB"
  }
  output {
    Array[File] chunks = glob("*set.xml")
    Int nchunks = length(chunks)
  }
}

# lightweight gather task that uses the counts in the input XMLs and skips
# reading index files, which is much faster and less memory-intensive
task gather_datasets_lite {
  input {
    Array[File] chunks
    String dataset_ext
    String? dataset_name
    String prefix = "gathered"

    Int nproc = 1
    String log_level = "DEBUG"
  }
  String dataset_out = "${prefix}${dataset_ext}"
  command {
    dataset \
      --log-level ${log_level} \
      --strict \
      create \
      --trustCounts \
      --no-sub-datasets \
      --unique-collections \
      ${"--name " + dataset_name} \
      ${dataset_out} \
      ${sep=" " chunks}
  }
  runtime {
    cpu: 1
    # even this is probably excessive but I am concerned about corner cases
    memory: "2GB"
  }
  output {
    File gathered = "${dataset_out}"
  }
}

task gather_generic {
  input {
    Array[File] chunks
    String output_file_name
    String? gather_args

    Int nproc = 1
    String log_level = "DEBUG"
  }
  command {
    python3 \
      -m pbcoretools.tasks.gather \
      --log-level ${log_level} \
      ${gather_args} \
      ${output_file_name} \
      ${sep=" " chunks}
  }
  runtime {
    cpu: 1
    memory: "1GB"
  }
  output {
    File gathered = "${output_file_name}"
  }
}

task gather_lima_datasets {
  input {
    Array[File] datastores

    Int nproc
    String log_level = "INFO"
  }
  command {
    python3 -m \
      pbcoretools.tasks.gather_lima_datasets \
      --log-level ${log_level} \
      --nproc ${nproc} \
      lima.datastore.json \
      ${sep=" " datastores} \
  }
  runtime {
    cpu: nproc
  }
  output {
    File gathered = "lima.datastore.json"
  }
}

task gather_bam {
  input {
    Array[File] bam_files
    String output_file_name = "merged.bam"

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    pbmerge -o ${output_file_name} ${sep=" " bam_files}
  }
  runtime {
    cpu: 1
  }
  output {
    File merged = "${output_file_name}"
  }
}

task gather_fofn {
  input {
    Array[File] fofn_files

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 -m \
      pbcoretools.tasks.gather_fofn \
      ${sep=" " fofn_files}
  }
  runtime {
    cpu: 1
    memory: "1GB"
  }
  output {
    File merged = "merged.bam"
  }
}
