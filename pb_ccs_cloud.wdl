# Standalone CCS workflow for running on cloud services without a shared
# filesystem.  Unlike our production workflows, this ditches the XML dataset
# format entirely, and uses the built-in chunking in the CCS CLI.

version 1.0

task ccs {
  input {
    File subreads_bam
    File subreads_bam_pbi
    Int min_passes = 3
    Float min_read_score = 0.65
    Float min_predicted_accuracy = 0.9
    Float min_snr = 2.5
    Int min_length = 10
    Int max_length = 50000
    Boolean by_strand = false

    Int nproc
    String log_level = "INFO"
    Int i_chunk = 1
    Int nchunks = 1
  }
  command {
    ln -s ${subreads_bam} subreads.bam
    ln -s ${subreads_bam_pbi} subreads.bam.pbi
    ccs \
      subreads.bam \
      ccs.bam \
      --chunk ${i_chunk + 1}/${nchunks} \
      --log-level ${log_level} \
      --minLength ${min_length} \
      --maxLength ${max_length} \
      --minPasses ${min_passes} \
      --minSnr ${min_snr} \
      --minPredictedAccuracy ${min_predicted_accuracy} \
      ${true="--by-strand" false="" by_strand} \
      --task-report task-report.json \
      -j ${nproc}
  }
  runtime {
    docker: "quay.io/biocontainers/pbccs:4.0.0--0"
    cpu: nproc
  }
  output {
    File bam = "ccs.bam"
  }
}

task gather_bam {
  input {
    Array[File] chunked_bams
  }
  command {
    pbmerge -o ccs.bam ${sep=" " chunked_bams}
  }
  runtime {
    docker: "quay.io/biocontainers/pbbam:1.0.6--hc16d5b3_1"
    cpu: 1
  }
  output {
    File merged_bam = "ccs.bam"
    File merged_bam_pbi = "ccs.bam.pbi"
  }
}

task export_fastq {
  input {
    File ccs_bam
    File ccs_bam_pbi
  }
  command {
    ln -s ${ccs_bam} ccs.bam
    ln -s ${ccs_bam_pbi} ccs.bam.pbi
    bam2fastq -o ccs ccs.bam
  }
  runtime {
    docker: "quay.io/biocontainers/bam2fastx:1.3.0--he1c1bb9_8"
    cpu: 1
  }
  output {
    # bam2fastx always compresses files
    File fastq = "ccs.fastq.gz"
  }
}

workflow pb_ccs_cloud {
  input {
    File subreads_bam
    File subreads_bam_pbi
    Float ccs_min_predicted_accuracy = 0.99
    Int ccs_min_passes = 3
    Float ccs_min_snr = 2.5
    Int ccs_min_length = 10
    Int ccs_max_length = 50000
    Boolean ccs_by_strand = false
    String dataset_filters = ""
    Int downsample_factor = 0

    Int nproc = 1
    Int nchunks = 1
  }

  scatter (i_chunk in range(nchunks)) {
    call ccs {
      input:
        subreads_bam = subreads_bam,
        subreads_bam_pbi = subreads_bam_pbi,
        min_passes = ccs_min_passes,
        min_predicted_accuracy = ccs_min_predicted_accuracy,
        min_snr = ccs_min_snr,
        min_length = ccs_min_length,
        max_length = ccs_max_length,
        by_strand = ccs_by_strand,
        nproc = nproc,
        i_chunk = i_chunk,
        nchunks = nchunks
    }
  }

  call gather_bam {
    input:
      chunked_bams = ccs.bam
  }

  call export_fastq {
    input:
      ccs_bam = gather_bam.merged_bam,
      ccs_bam_pbi = gather_bam.merged_bam_pbi
  }

  output {
    File ccs_bam = gather_bam.merged_bam
    File ccs_bam_pbi = gather_bam.merged_bam_pbi
    File ccs_fastq = export_fastq.fastq
  }
}
