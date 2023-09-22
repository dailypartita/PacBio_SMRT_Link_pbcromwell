# HGAP4
version 1.0

import "wf_falcon.wdl"
import "wf_consensus.wdl"
import "wf_mapping.wdl"
import "wf_coverage_reports.wdl"
import "tasks/pbcoretools.wdl" as pbcoretools
import "tasks/pbreports.wdl" as pbreports

task task_get_fastas {
  input {
    File fofn # for now, all abspaths
  }
  Array[File] fns = read_lines(fofn)

  command <<<
  >>>
  output {
    Array[File] array_of_fasta = fns
  }
}

task task_get_dextas {
  input {
    File subreadset_xml # can have multiple BAMs
    Int nproc = 1
  }
  command <<<
    set -vex
    env | sort

    # dataset absolutize --outdir . ... # would fail because it will not traverse the symlink
    rfn=$(readlink -f ~{subreadset_xml})
    python3 -m falcon_kit.mains.bam2dexta --nproc=~{nproc}  split-apply-combine \
        --bam=${rfn} --dexta=./dexta.fofn

    cat ./dexta.fofn
    ls -l uow-*/*.dexta
  >>>
  runtime {
    cpu: nproc
  }
  output {
    Array[File] array_of_dexta = glob("uow-*/*.dexta")
  }
}
task task_gen_config {
  input {
    Array[File] array_of_fasta # can be .dexta or .fasta files
    # These are effectively workflow-options.
    Boolean HGAP_AggressiveAsm_bool
    String HGAP_GenomeLength_str
    String HGAP_SeedCoverage_str
    String HGAP_SeedLengthCutoff_str
    String HGAP_FalconAdvanced_str
  }
  Object object_opts = object {
    HGAP_AggressiveAsm_bool: HGAP_AggressiveAsm_bool,
    HGAP_GenomeLength_str: HGAP_GenomeLength_str,
    HGAP_SeedCoverage_str: HGAP_SeedCoverage_str,
    HGAP_SeedLengthCutoff_str: HGAP_SeedLengthCutoff_str,
    HGAP_FalconAdvanced_str: HGAP_FalconAdvanced_str,
  }
  File hgap4_options_json = write_json(object_opts)
  command <<<
    set -vex
    ln -f ~{hgap4_options_json} ./hgap4_options.json
    python3 -m pbfalcon.gen_config --i-options-fn ./hgap4_options.json --o-cfg-fn ./fc_run.cfg --o-json-fn ./General_config.json
  >>>
  output {
    File cfg = "fc_run.cfg"
    File config = "General_config.json"
  }
}

task fasta_to_reference {
  input {
    File fasta

    String log_level = "INFO"
  }
  command {
    set -vex
    ln -s ${fasta} falcon_draft_assembly.fasta
    dataset \
      --log-level ${log_level} \
      create \
      --generateIndices \
      --type ReferenceSet \
      --name falcon_draft_assembly \
      falcon_draft_assembly.referenceset.xml \
      falcon_draft_assembly.fasta
  }
  runtime {
    cpu: 1
  }
  output {
    File reference_xml = "falcon_draft_assembly.referenceset.xml"
  }
}

workflow pb_hgap4 {
  input {
    File eid_subread
    #File ifile_input_fofn = "/scratch/cdunn/pbcromwell/abs_input.fofn"
    ##File ifile_config = "/scratch/cdunn/pbcromwell/General_config.json"
    ##File ifile_config = "General_config.json"
    Boolean hgap4_aggressive_asm = false
    Int hgap4_genome_length     = 5000000
    Int hgap4_seed_coverage     = 30
    Int hgap4_seed_length_cutoff = -1
    String hgap4_falcon_advanced = ""
    String consensus_algorithm = "arrow"
    String dataset_filters = ""
    Int downsample_factor = 0
    Float mapping_min_concordance = 70
    Int mapping_min_length = 50
    String? mapping_biosample_name
    String? mapping_pbmm2_overrides

    Int nproc = 8
    Int max_nchunks = 100
    String? tmp_dir
    String log_level = "INFO"
  }

  ##Map[String, String] cfg = read_json(ifile_config)
  ##Array[File] array_of_fasta = read_lines(ifile_input_fofn)

  call pbcoretools.dataset_filter {
    input:
      dataset = eid_subread,
      filters = dataset_filters,
      downsample_factor = downsample_factor,
      #nproc = 1 # currently ignored anyway
  }
  #call task_get_fastas {
  #  input:
  #    fofn = ifile_input_fofn,
  #}
  call task_get_dextas {
    input:
      #subreadset_xml = pbcoretools.dataset_filter.filtered, # does not work
      subreadset_xml = dataset_filter.filtered,
      nproc          = nproc,
  }
  call task_gen_config {
    input:
      #array_of_fasta = task_get_fastas.array_of_fasta,
      array_of_fasta = task_get_dextas.array_of_dexta,  # can be dexta or fasta
      HGAP_AggressiveAsm_bool = hgap4_aggressive_asm,
      HGAP_GenomeLength_str = hgap4_genome_length,
      HGAP_SeedCoverage_str = hgap4_seed_coverage,
      HGAP_SeedLengthCutoff_str = hgap4_seed_length_cutoff,
      HGAP_FalconAdvanced_str = hgap4_falcon_advanced
  }
  call wf_falcon.falcon {
    input:
      ifile_config   = task_gen_config.config,
      #array_of_fasta = task_get_fastas.array_of_fasta,
      array_of_fasta = task_get_dextas.array_of_dexta,
      nproc          = nproc,
      max_nchunks    = max_nchunks
  }

  call fasta_to_reference {
    input:
      fasta = falcon.ofile_p_ctg_fasta,
      log_level = log_level
  }

  call wf_mapping.mapping {
    input:
      reads = dataset_filter.filtered,
      reference = fasta_to_reference.reference_xml,
      mapping_stats_module = "mapping_stats_hgap",
      min_concordance = mapping_min_concordance,
      min_length = mapping_min_length,
      biosample_name = mapping_biosample_name,
      pbmm2_overrides = mapping_pbmm2_overrides,
      run_coverage = false,
      nproc = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir
  }

  call wf_coverage_reports.coverage_reports {
    input:
      mapped = mapping.mapped,
      reference = fasta_to_reference.reference_xml,
      nproc = nproc,
      log_level = log_level
  }

  call wf_consensus.consensus {
    input:
      alignments = mapping.mapped,
      reference = fasta_to_reference.reference_xml,
      consensus_algorithm = "arrow",
      log_level = log_level,
      nproc = nproc,
      max_nchunks = max_nchunks
  }

  call pbreports.polished_assembly {
    input:
      coverage_gff = coverage_reports.coverage_gff,
      consensus_fastq = consensus.consensus_fastq,
      nproc = nproc,
      log_level = log_level
  }

  output {
    File ofile_p_ctg_fasta = falcon.ofile_p_ctg_fasta
    File ofile_a_ctg_fasta = falcon.ofile_a_ctg_fasta
    File mapped = mapping.mapped
    File? report_mapping_stats = mapping.report_mapping_stats
    File coverage_gff = coverage_reports.coverage_gff
    File report_coverage = coverage_reports.report_coverage
    File consensus_fasta = consensus.consensus_fasta
    File consensus_fastq = consensus.consensus_fastq
    File report_preassembly = falcon.report_preassembly
    File report_polished_assembly = polished_assembly.report
  }
}
