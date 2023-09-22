# Variant calling using gcpp

version 1.0

import "tasks/chunking.wdl" as chunking
import "tasks/pbreports.wdl" as pbreports

task guess_optimal_max_nchunks {
  input {
    File reference
    Int max_nchunks

    String log_level = "INFO"
  }
  Int default_nchunks = 1
  command <<<
    python3 <<EOF
    from pbcore.io import ReferenceSet
    from pbcoretools.chunking.chunk_utils import guess_optimal_max_nchunks_for_consensus
    import os.path
    ds = ReferenceSet(os.path.realpath("~{reference}"))
    nchunks = guess_optimal_max_nchunks_for_consensus(ds.totalLength, ~{max_nchunks})
    with open("nchunks.txt", "w") as txt_out:
      txt_out.write(str(nchunks))
    EOF
  >>>
  runtime {
    cpu: 1
    backend: "Local"
    memory: "1GB"
  }
  output {
    File nchunks_txt = "nchunks.txt"
  }
}

task gcpp {
  input {
    File mapped
    File reference
    String consensus_algorithm = "arrow"
    Int min_coverage = 5
    Int max_coverage = 100
    Int min_confidence = 40
    Boolean report_effective_coverage = false
    Int genome_length_mb
    String? advanced_args

    Int nproc
    String log_level = "INFO"
  }
  command {
    set -e
    gcpp \
      -j ${nproc} \
      --algorithm ${consensus_algorithm} \
      ${advanced_args} \
      ${"--min-coverage " + min_coverage} \
      ${"--coverage " + max_coverage} \
      ${"--min-confidence " + min_confidence} \
      ${true="--report-effective-coverage" false="" report_effective_coverage} \
      --alarms alarms.json \
      -r `readlink -f ${reference}` \
      `readlink -f ${mapped}` \
      -o consensus.fasta,consensus.fastq,consensus.gff,consensus.vcf

    # GCCP does not generate all the output files when the input is empty.
    touch consensus.vcf
    touch consensus.gff
    touch consensus.fastq
    touch consensus.fasta
  }
  Int base_mem_mb = 2048
  Int total_mem_mb = base_mem_mb + (genome_length_mb * 3)
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File fasta = "consensus.fasta"
    File fastq = "consensus.fastq"
    File gff = "consensus.gff"
    File vcf = "consensus.vcf"
  }
}

workflow consensus {
  input {
    File alignments
    File reference
    String consensus_algorithm = "arrow"
    Int min_coverage = 5
    Int max_coverage = 100
    Int min_confidence = 40
    String? gcpp_advanced_args
    Boolean report_effective_coverage = false
    Int? genome_length_mb

    Int nproc = 8
    Int max_nchunks = 100
    Int? target_size
    String log_level = "INFO"
  }

  call guess_optimal_max_nchunks {
    input:
      reference = reference,
      max_nchunks = max_nchunks,
      log_level = log_level
  }

  Int nchunks = read_int(guess_optimal_max_nchunks.nchunks_txt)
  Int genome_length_default = 3300
  # FIXME this should probably take nchunks into account
  Int genome_length_mb_final = select_first([genome_length_mb,
                                             genome_length_default])

  call chunking.split_alignments {
    input:
      mapped = alignments,
      max_nchunks = nchunks,
      target_size = target_size,
      log_level = log_level
  }

  scatter (alignment_chunk in split_alignments.chunks) {
    call gcpp as genomic_consensus {
      input:
        mapped = alignment_chunk,
        reference = reference,
        consensus_algorithm = consensus_algorithm,
        min_coverage = min_coverage,
        max_coverage = max_coverage,
        min_confidence = min_confidence,
        advanced_args = gcpp_advanced_args,
        report_effective_coverage = report_effective_coverage,
        genome_length_mb = genome_length_mb_final,
        nproc = nproc,
        log_level = log_level
    }
  }

  call chunking.gather_generic as gather_gff {
    input:
      chunks = genomic_consensus.gff,
      output_file_name = "variants.gff",
      log_level = log_level
  }

  call chunking.gather_generic as gather_vcf {
    input:
      chunks = genomic_consensus.vcf,
      output_file_name = "variants.vcf",
      log_level = log_level
  }

  call chunking.gather_generic as gather_fasta {
    input:
      chunks = genomic_consensus.fasta,
      output_file_name = "consensus.fasta",
      gather_args = "--join-contigs",
      log_level = log_level
  }

  call chunking.gather_generic as gather_fastq {
    input:
      chunks = genomic_consensus.fastq,
      output_file_name = "consensus.fastq",
      gather_args = "--join-contigs",
      log_level = log_level
  }

  output {
    File consensus_fasta = gather_fasta.gathered
    File consensus_fastq = gather_fastq.gathered
    File variants_gff = gather_gff.gathered
    File variants_vcf = gather_vcf.gathered
  }
}
