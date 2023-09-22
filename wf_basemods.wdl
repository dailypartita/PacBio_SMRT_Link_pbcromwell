# Core base modification analysis, starting from mapped subreads

version 1.0

import "tasks/pbcoretools.wdl" as pbcoretools
import "tasks/chunking.wdl" as chunking

task ipdsummary {
  input {
    File mapped
    File reference
    Float p_value = 0.001
    String? identify_mods
    Int? max_length
    Boolean write_csv = false
    Boolean write_bigwig = true
    Int base_memory_mb = 0

    Int nproc
    String log_level = "DEBUG"
  }
  command {
    ipdSummary \
      --log-level ${log_level} \
      --log-file ipdSummary.log \
      --numWorkers ${nproc} \
      ${mapped} \
      --reference `readlink -f ${reference}` \
      --gff basemods.gff \
      ${true="--csv basemods.csv" false="" write_csv} \
      ${true="--bigwig ipds.bw" false="" write_bigwig} \
      --pvalue ${p_value} \
      --alignmentSetRefWindows \
      ${"--identify " + identify_mods} \
      ${"--maxLength " + max_length}
  }
  Int total_mem_mb = base_memory_mb + 4096 * nproc
  runtime {
    cpu: nproc
    # this is wildly excessive in most cases
    memory: "${total_mem_mb}MB"
  }
  output {
    File basemods_gff = "basemods.gff"
    File? basemods_csv = "basemods.csv"
    File? bigwig_file = "ipds.bw"
  }
}

task summarize_modifications {
  input {
    File basemods_gff
    File alignment_gff
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "DEBUG"
  }
  command {
    python3 \
      -m kineticsTools.summarizeModifications \
      --log-level ${log_level} \
      ${basemods_gff} \
      ${alignment_gff} \
      alignment_summary_with_basemods.gff
  }
  Int total_mem_mb = base_memory_mb + 2048
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File summary_basemods_gff = "alignment_summary_with_basemods.gff"
  }
}

task find_motifs {
  input {
    File basemods_gff
    File reference_set
    Int? min_qmod_score
    Int base_memory_mb = 0

    Int nproc = 4
    String log_level = "DEBUG"
  }
  # FIXME motifMaker won't accept ReferenceSet directly
  command <<<
    python3 <<EOF
    from pbcore.io import ReferenceSet
    import os.path
    fasta_file = ReferenceSet(os.path.realpath("~{reference_set}"),
                              skipCounts=True).toExternalFiles()[0]
    os.symlink(fasta_file, "reference.fasta")
    EOF
    motifMaker find \
      -j ~{nproc} \
      --gff ~{basemods_gff} \
      --fasta reference.fasta \
      --output motifs.csv \
      ~{"--minScore " + min_qmod_score}
  >>>
  Int total_mem_mb = base_memory_mb + 8192
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File motifs_csv = "motifs.csv"
    File reference_fasta = "reference.fasta"
  }
}

task reprocess_motifs {
  input {
    File basemods_gff
    File motifs_csv
    File reference_fasta
    Float? min_fraction
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "DEBUG"
  }
  command {
    motifMaker reprocess \
      --gff ${basemods_gff} \
      --motifs ${motifs_csv} \
      --fasta ${reference_fasta} \
      --output motifs.gff \
      ${"--minFraction " + min_fraction}
  }
  Int total_mem_mb = base_memory_mb + 4096
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File motifs_gff = "motifs.gff"
  }
}

task gather_bigwig {
  input {
    Array[File] bigwig_files
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "DEBUG"
  }
  command {
    python3 \
      -m kineticsTools.tasks.gather_bigwig \
      --log-level ${log_level} \
      ipds.bw \
      ${sep=" " bigwig_files}
  }
  Int total_mem_mb = base_memory_mb + 8192
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File gathered = "ipds.bw"
  }
}

task modifications_report {
  input {
    File basemods_csv
    Int? modqv_cutoff
    Int genome_length_mb
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 \
      -m pbreports.report.modifications \
      --log-level ${log_level} \
      ${"--modqv-cutoff " + modqv_cutoff} \
      ${basemods_csv} \
      modifications.report.json
  }
  # the high overhead appears to be a numpy deficiency
  Int total_mem_mb = base_memory_mb + 2048 + genome_length_mb * 25
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File report = "modifications.report.json"
    Array[File?] plot_pngs = glob("*.png")
  }
}

task motifs_report {
  input {
    File motifs_gff
    File motifs_csv
    File? kinetics_csv
    File? reference
    Int max_motifs = 10
    Int? modqv_cutoff
    Int genome_length_mb
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 \
      -m pbreports.report.motifs \
      --log-level ${log_level} \
      ${motifs_gff} ${motifs_csv} \
      motifs.report.json \
      --maxMotifs ${max_motifs} \
      ${"--modqv-cutoff " + modqv_cutoff} \
      ${"--kinetics-csv " + kinetics_csv} \
      ${"--reference `readlink -f " + reference + "`"}
  }
  Int total_mem_mb = base_memory_mb + 2048 + genome_length_mb * 640
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File report = "motifs.report.json"
    Array[File?] plot_pngs = glob("*.png")
  }
}

workflow wf_basemods {
  input {
    File alignments
    File reference
    Boolean run_find_motifs = false

    String kineticstools_identify_mods = "m4C,m6A"
    Float kineticstools_p_value = 0.001
    Int motif_min_score = 100
    Float motif_min_fraction = 0.30

    Int? genome_length_mb
    Int base_memory_mb = 0

    Int nproc = 1
    Int max_nchunks = 100
    Int target_size = 20000
    String log_level = "INFO"
    String? tmp_dir
  }

  Boolean write_bigwig = true
  Int max_genome_length_mb_for_csv = 20

  Int dataset_size = select_first([genome_length_mb,
                                   max_genome_length_mb_for_csv])
  Boolean write_csv = dataset_size <= max_genome_length_mb_for_csv

  call chunking.split_alignments {
    input:
      mapped = alignments,
      max_nchunks = max_nchunks,
      target_size = target_size
  }

  scatter (aligned_chunk in split_alignments.chunks) {
    call ipdsummary {
      input:
        mapped = aligned_chunk,
        reference = reference,
        identify_mods = kineticstools_identify_mods,
        p_value = kineticstools_p_value,
        write_bigwig = write_bigwig,
        write_csv = write_csv,
        base_memory_mb = base_memory_mb,
        nproc = nproc,
        log_level = log_level
    }
  }

  call chunking.gather_generic as gather_basemods_gff {
    input:
      chunks = ipdsummary.basemods_gff,
      output_file_name = "basemods.gff",
      log_level = log_level
  }

  if (write_bigwig) {
    call gather_bigwig {
      input:
        bigwig_files = select_all(ipdsummary.bigwig_file),
        base_memory_mb = base_memory_mb,
        log_level = log_level
    }
  }

  if (write_csv) {
    call chunking.gather_generic as gather_basemods_csv {
      input:
        chunks = select_all(ipdsummary.basemods_csv),
        output_file_name = "basemods.csv",
        log_level = log_level
    }

    call modifications_report {
      input:
        basemods_csv = gather_basemods_csv.gathered,
        modqv_cutoff = motif_min_score,
        genome_length_mb = dataset_size,
        base_memory_mb = base_memory_mb,
        nproc = nproc,
        log_level = log_level
    }
  }

  if (run_find_motifs) {
    Int MOTIF_NPROC_MAX = 4
    Int motif_nproc = if (nproc > MOTIF_NPROC_MAX) then MOTIF_NPROC_MAX else nproc
    call find_motifs {
      input:
        basemods_gff = gather_basemods_gff.gathered,
        reference_set = reference,
        min_qmod_score = motif_min_score,
        base_memory_mb = base_memory_mb,
        nproc = motif_nproc,
        log_level = log_level
    }

    call reprocess_motifs {
      input:
        basemods_gff = gather_basemods_gff.gathered,
        motifs_csv = find_motifs.motifs_csv,
        reference_fasta = find_motifs.reference_fasta,
        min_fraction = motif_min_fraction,
        base_memory_mb = base_memory_mb
    }

    call motifs_report {
      input:
        motifs_gff = reprocess_motifs.motifs_gff,
        motifs_csv = find_motifs.motifs_csv,
        modqv_cutoff = motif_min_score,
        kinetics_csv = gather_basemods_csv.gathered,
        reference = reference,
        genome_length_mb = dataset_size,
        base_memory_mb = base_memory_mb,
        log_level = log_level
    }
  }

  output {
    File basemods_gff = gather_basemods_gff.gathered
    File? basemods_csv = gather_basemods_csv.gathered
    File? report_modifications = modifications_report.report
    File? bigwig_file = gather_bigwig.gathered
    File? motifs_csv = find_motifs.motifs_csv
    File? motifs_gff = reprocess_motifs.motifs_gff
    File? report_motifs = motifs_report.report
  }
}
