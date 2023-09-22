# Utility for getting input sizes in order to determine approximate memory
# consumption.

version 1.0

task get_dataset_size {
  input {
    # not actually optional
    File? dataset
    Boolean skip_counts = true
    Boolean get_index_size = false

    Int nproc = 1
    String log_level = "INFO"
  }
  File dataset_xml = select_first([dataset])
  command {
    python3 -m pbcoretools.tasks.memory.get_dataset_size \
      `readlink -f ${dataset}` \
      ${true="--skip-counts" false="--no-skip-counts" skip_counts} \
      ${true="--get-index-size" false="" get_index_size}
  }
  runtime {
    cpu: 1
    # this is assuming skipCounts=True for BAM datasets
    memory: "1GB"
  }
  output {
    # workaround for SL-5656
    File num_records_txt = "numrecords.txt"
    File total_length_txt = "totallength.txt"
    File index_size_txt = "indexsize.txt"
    File num_files_txt = "numresources.txt"
    File num_filters_txt = "numfilters.txt"
  }
}

workflow get_input_sizes {
  input {
    File bam_dataset
    File? fasta_dataset
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }

  call get_dataset_size as get_bam_size {
    input:
      dataset = bam_dataset,
      skip_counts = true,
      get_index_size = true,
      log_level = log_level
  }

  Int bam_num_records = read_int(get_bam_size.num_records_txt)
  Int bam_total_length_mb = read_int(get_bam_size.total_length_txt)
  Int bam_index_size = read_int(get_bam_size.index_size_txt)
  Int bam_num_files = read_int(get_bam_size.num_files_txt)
  Int bam_num_filters = read_int(get_bam_size.num_filters_txt)

  if (defined(fasta_dataset)) {
    call get_dataset_size as get_ref_size {
      input:
        dataset = fasta_dataset,
        skip_counts = false,
        get_index_size = false,
        log_level = log_level
    }
    Int ref_num_records = read_int(get_ref_size.num_records_txt)
    Int ref_total_length = read_int(get_ref_size.total_length_txt)
  }

  output {
    Int bam_size = bam_num_records
    Int index_memory_gb = bam_index_size
    Int bam_length_mb = bam_total_length_mb
    Int bam_files = bam_num_files
    Int bam_filters = bam_num_filters
    Int? genome_length_mb = ref_total_length
    Int? num_contigs = ref_num_records
  }
}
