# Mock version of Export Reads workflow

version 1.0

task mock_export {
  input {
    File eid_ccs
    Boolean output_bam = false
    Boolean output_fastq = true
    Boolean output_fasta = true
  }
  command {
    if [ "${output_bam}" = "true" ]; then
      echo "1" > hifi_reads.bam
    fi
    if [ "${output_fastq}" = "true" ]; then
      echo "2" > hifi_reads.fastq
    fi
    if [ "${output_fasta}" = "true" ]; then
      echo "3" > hifi_reads.fasta
    fi
  }
  runtime {
    cpu: 1
    memory: "1GB"
  }
  output {
    File? bam = "hifi_reads.bam"
    File? fastq = "hifi_reads.fastq"
    File? fasta = "hifi_reads.fasta"
  }
}

workflow dev_mock_export_ccs {
  input {
    File eid_ccs
    Boolean output_bam = false
    Boolean output_fastq = true
    Boolean output_fasta = true

    Int filter_min_qv = 20

    # Included here for consistency but this is not a parallel workflow
    Int max_nchunks = 1
    Int nproc = 1
    String log_level = "INFO"
    String? tmp_dir
  }

  call mock_export {
    input:
      eid_ccs = eid_ccs,
      output_bam = output_bam,
      output_fastq = output_fastq,
      output_fasta = output_fasta
  }

  output {
    File? bam = mock_export.bam
    File? fastq = mock_export.fastq
    File? fasta = mock_export.fasta
  }
}
