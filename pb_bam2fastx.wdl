version 1.0

import "tasks/pbcoretools.wdl" as pbcoretools
import "tasks/memory.wdl" as memory

task bam2fasta_archive {
  input {
    File reads
    String prefix = "reads"

    Int nproc = 1
    String log_level = "DEBUG"
    Int memory_gb = 8
  }
  command {
    python3 \
      -m pbcoretools.tasks.bam2fasta_archive \
      --log-level ${log_level} \
      ${reads} \
      ${prefix}.fasta.tar.gz
  }
  runtime {
    cpu: 1
    memory: "${memory_gb}GB"
  }
  output {
    File fasta_tgz = "${prefix}.fasta.tar.gz"
  }
}

task bam2fastq_archive {
  input {
    File reads
    String prefix = "reads"

    Int nproc = 1
    String log_level = "DEBUG"
    Int memory_gb = 8
  }
  command {
    python3 \
      -m pbcoretools.tasks.bam2fastq_archive \
      --log-level ${log_level} \
      ${reads} \
      ${prefix}.fastq.tar.gz
  }
  runtime {
    cpu: 1
    memory: "${memory_gb}GB"
  }
  output {
    File fastq_tgz = "${prefix}.fastq.tar.gz"
  }
}

workflow pb_bam2fastx {
  input {
    File eid_subread
    String dataset_filters = ""
    Int downsample_factor = 0

    # ignored
    Int max_nchunks = 0
    Int nproc = 1
    String log_level = "INFO"
    String? tmp_dir
  }

  call memory.get_input_sizes {
    input:
      bam_dataset = eid_subread,
      log_level = log_level
  }

  call pbcoretools.dataset_filter {
    input:
      dataset = eid_subread,
      filters = dataset_filters,
      downsample_factor = downsample_factor,
      memory_gb = get_input_sizes.index_memory_gb,
      log_level = log_level
  }

  call bam2fasta_archive {
    input:
      reads = dataset_filter.filtered,
      memory_gb = get_input_sizes.index_memory_gb,
      log_level = log_level
  }

  call bam2fastq_archive {
    input:
      reads = dataset_filter.filtered,
      memory_gb = get_input_sizes.index_memory_gb,
      log_level = log_level
  }

  output {
    File fasta_tgz = bam2fasta_archive.fasta_tgz
    File fastq_tgz = bam2fastq_archive.fastq_tgz
  }
}
