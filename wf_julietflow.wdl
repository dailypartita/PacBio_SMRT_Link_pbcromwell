# standalone single-sample minor variants workflow

version 1.0

import "tasks/pbmm2.wdl" as mapping_tasks
import "tasks/pbreports.wdl" as pbreports

task fuse {
  input {
    File alignments

    Int nproc = 1
  }
  command {
    set -e
    fuse \
      --alarms alarms.json \
      `readlink -f ${alignments}` \
      fuse.fasta
    samtools faidx fuse.fasta
    dataset --strict create --type ReferenceSet --name "Alignment-based reference" fuse.referenceset.xml fuse.fasta
  }
  runtime {
    cpu: nproc
    memory: "8GB"
  }
  output {
    File reference_new = "fuse.referenceset.xml"
  }
}

task cleric {
  input {
    File alignments
    File reference
    File reference_fuse

    Int nproc = 1
  }
  command {
    cleric \
      --alarms alarms.json \
      `readlink -f ${alignments}` \
      `readlink -f ${reference}` \
      `readlink -f ${reference_fuse}` \
      cleric.consensusalignmentset.xml
  }
  runtime {
    cpu: 1
    memory: "8GB"
  }
  output {
    File mapped = "cleric.consensusalignmentset.xml"
  }
}

task juliet {
  input {
    File alignments
    String target_config = "none"
    Boolean mode_phasing = false
    String genomic_region = ""
    Float minimal_percentage = 0.1
    Float maximal_percentage = 100
    Boolean only_known_drms = false
    Float substitution_rate = 0.0
    Float deletion_rate = 0.0
    Boolean debug = false

    Int nproc = 1
  }
  # FIXME
  #   ${"--region " + genomic_region} \
  command {
    juliet \
      --alarms alarms.json \
      ${"--config " + target_config} \
      ${true="--mode-phasing" false="" mode_phasing} \
      ${true="--drm-only" false="" only_known_drms} \
      ${"--min-perc " + minimal_percentage} \
      ${"--max-perc " + maximal_percentage} \
      ${"--sub " + substitution_rate} \
      ${"--del " + deletion_rate} \
      ${true="--debug" false="" debug} \
      `readlink -f ${alignments}`
  }
  runtime {
    cpu: 1
    memory: "8GB"
  }
  output {
    File html = "cleric.html"
    File json = "cleric.json"
  }
}

workflow julietflow {
  input {
    File ccs
    File reference
    Boolean juliet_debug = false
    Float deletion_rate = 0.0
    Float minimal_percentage = 0.1
    Float maximal_percentage = 100.0
    Boolean mode_phasing = true
    Boolean only_known_drms = false
    String genomic_region = ""
    Float substitution_rate = 0.0
    String target_config = "none"
    String? pbmm2_overrides

    Int nproc
    String log_level = "DEBUG"
    String? tmp_dir
  }

  call mapping_tasks.pbmm2_align as align1 {
    input:
      unmapped = ccs,
      reference = reference,
      aln_ext = ".consensusalignmentset.xml",
      preset_mode = "CCS",
      pbmm2_overrides = pbmm2_overrides,
      nproc = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir
  }

  call fuse {
    input:
      alignments = align1.mapped,
      nproc = nproc
  }

  call mapping_tasks.pbmm2_align as align2 {
    input:
      unmapped = ccs,
      reference = fuse.reference_new,
      aln_ext = ".consensusalignmentset.xml",
      preset_mode = "CCS",
      nproc = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir
  }

  call cleric {
    input:
      alignments = align2.mapped,
      reference = reference,
      reference_fuse = fuse.reference_new
  }

  call juliet {
    input:
      alignments = cleric.mapped,
      target_config = target_config,
      mode_phasing = mode_phasing,
      genomic_region = genomic_region,
      minimal_percentage = minimal_percentage,
      maximal_percentage = maximal_percentage,
      only_known_drms = only_known_drms,
      substitution_rate = substitution_rate,
      deletion_rate = deletion_rate,
      debug = juliet_debug,
      nproc = nproc
  }

  output {
    File juliet_html = juliet.html
    File juliet_json = juliet.json
    File mapped = align2.mapped
  }
}
