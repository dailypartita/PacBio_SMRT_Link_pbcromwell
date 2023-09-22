# Combined HiFi microbial assembly and basemods analysis

version 1.0

import "pb_assembly_hifi_microbial.wdl"
import "pb_basemods_hifi.wdl"

task make_reference {
  input {
    File assembly_fasta
    File? assembly_fasta_fai

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    set -e
    if [ -z "${assembly_fasta_fai}" ]; then
      cp ${assembly_fasta} assembly_reference.fasta
      samtools faidx assembly_reference.fasta
    else
      ln -s ${assembly_fasta} assembly_reference.fasta
      ln -s ${assembly_fasta_fai} assembly_reference.fasta.fai
    fi
    dataset \
      --log-level DEBUG \
      create \
      --type ReferenceSet \
      --name "Final Assembly" \
      assembly.referenceset.xml \
      assembly_reference.fasta
  }
  runtime {
    cpu: 1
    memory: "1GB"
  }
  output {
    File reference = "assembly.referenceset.xml"
  }
}

task evaluate_kinetics {
  input {
    File assembly_fasta
    File eid_ccs

    Int nproc = 1
    String log_level = "DEBUG"
  }
  command {
    python3 -m pbcoretools.tasks.evaluate_kinetics \
      --log-level ${log_level} \
      --log-file evaluate_kinetics.log \
      `readlink -f ${eid_ccs}` \
      `readlink -f ${assembly_fasta}`
  }
  runtime {
    cpu: 1
    memory: "1GB"
  }
  output {
    File? success_fn = "success.txt"
    File? alarms_json = "alarms.json"
  }
}

workflow pb_microbial_analysis {
  input {
    File eid_ccs

    Int filter_min_qv = 20
    String dataset_filters = ""

    String ipa2_genome_size = "10M"
    Int ipa2_downsampled_coverage = 100
    # the real defaults are set below
    String ipa2_advanced_options_chrom = ""
    String ipa2_advanced_options_plasmid = ""
    Boolean ipa2_cleanup_intermediate_files = true
    Int microasm_plasmid_contig_len_max = 300000 # 300kB
    Boolean microasm_run_secondary_polish = true

    Boolean run_basemods = true
    String kineticstools_identify_mods = "m4C,m6A"
    Float kineticstools_p_value = 0.001
    # aka Qmod, this is 100 in the subreads workflow, but empirically 35 is a
    # good cutoff for HiFi (tuned for sensitivity).  plots in jira:SL-7479
    Int motif_min_score = 35
    Float motif_min_fraction = 0.30
    Boolean run_find_motifs = true
    Int add_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
    # Use max_nchunks = 1 for a local run to avoid resource competition. On a cluster, the larger the better.
    Int max_nchunks = 40
    String tmp_dir = "/tmp"
  }

  # Internal defaults for advanced options - if the workflow options are
  # null or empty, these will be applied automatically
  String DEFAULT_OPTIONS_CHROM = "config_block_size = 100; config_seeddb_opt = -k 28 -w 20 --space 0 --use-hpc-seeds-only; config_ovl_opt = --one-hit-per-target --min-idt 98 --traceback --mask-hp --mask-repeats --trim --trim-window-size 30 --trim-match-frac 0.75;"
  String DEFAULT_OPTIONS_PLASMID = "config_block_size = 100; config_ovl_filter_opt = --max-diff 80 --max-cov 100 --min-cov 2 --bestn 10 --min-len 500 --gapFilt --minDepth 4 --idt-stage2 98; config_ovl_min_len = 500; config_seeddb_opt = -k 28 -w 20 --space 0 --use-hpc-seeds-only; config_ovl_opt = --one-hit-per-target --min-idt 98 --min-map-len 500 --min-anchor-span 500 --traceback --mask-hp --mask-repeats --trim --trim-window-size 30 --trim-match-frac 0.75 --smart-hit-per-target --secondary-min-ovl-frac 0.05; config_layout_opt = --allow-circular;"

  if (ipa2_advanced_options_chrom == "") {
    String internal_advanced_options_chrom = DEFAULT_OPTIONS_CHROM
  }
  if (ipa2_advanced_options_plasmid == "") {
    String internal_advanced_options_plasmid = DEFAULT_OPTIONS_PLASMID
  }

  String ipa2_advanced_options_chrom_final = select_first(
    [internal_advanced_options_chrom, ipa2_advanced_options_chrom])
  String ipa2_advanced_options_plasmid_final = select_first(
    [internal_advanced_options_plasmid, ipa2_advanced_options_plasmid])

  call pb_assembly_hifi_microbial.pb_assembly_hifi_microbial {
    input:
      eid_ccs = eid_ccs,
      ipa2_genome_size = ipa2_genome_size,
      ipa2_downsampled_coverage = ipa2_downsampled_coverage,
      ipa2_advanced_options_chrom = ipa2_advanced_options_chrom_final,
      ipa2_advanced_options_plasmid = ipa2_advanced_options_plasmid_final,
      ipa2_cleanup_intermediate_files = ipa2_cleanup_intermediate_files,
      microasm_plasmid_contig_len_max = microasm_plasmid_contig_len_max,
      microasm_run_secondary_polish = microasm_run_secondary_polish,
      add_memory_mb = add_memory_mb,
      nproc = nproc,
      max_nchunks = max_nchunks,
      log_level = log_level,
      tmp_dir = tmp_dir
  }

  if (run_basemods) {
    call evaluate_kinetics {
      input:
        assembly_fasta = pb_assembly_hifi_microbial.assembly_fasta,
        eid_ccs = eid_ccs
    }
  }

  # only run basemods on non-empty reference (with index)
  if (defined(evaluate_kinetics.success_fn)) {
    call make_reference {
      input:
        assembly_fasta = pb_assembly_hifi_microbial.assembly_fasta,
        assembly_fasta_fai = pb_assembly_hifi_microbial.assembly_fasta_fai
    }

    call pb_basemods_hifi.pb_basemods_hifi {
      input:
        eid_ccs = eid_ccs,
        eid_ref_dataset = make_reference.reference,
        kineticstools_identify_mods = kineticstools_identify_mods,
        kineticstools_p_value = kineticstools_p_value,
        run_find_motifs = run_find_motifs,
        motif_min_score = motif_min_score,
        motif_min_fraction = motif_min_fraction,
        max_nchunks = max_nchunks,
        add_memory_mb = add_memory_mb,
        nproc = nproc,
        log_level = log_level,
        tmp_dir = tmp_dir
    }
  }

  output {
    File assembly_fasta = pb_assembly_hifi_microbial.assembly_fasta
    File? assembly_fasta_fai = pb_assembly_hifi_microbial.assembly_fasta_fai
    File ncbi_fasta = pb_assembly_hifi_microbial.ncbi_fasta
    File circular_list = pb_assembly_hifi_microbial.circular_list
    File mapped = pb_assembly_hifi_microbial.mapped
    File mapped_target_fasta = pb_assembly_hifi_microbial.mapped_target_fasta
    File? mapped_target_fasta_fai = pb_assembly_hifi_microbial.mapped_target_fasta_fai
    File report_polished_assembly = pb_assembly_hifi_microbial.report_polished_assembly
    File coverage_gff = pb_assembly_hifi_microbial.coverage_gff
    File report_mapping_stats = pb_assembly_hifi_microbial.report_mapping_stats
    File report_coverage = pb_assembly_hifi_microbial.report_coverage
    File? mapped_bam_datastore = pb_assembly_hifi_microbial.mapped_bam_datastore
    File? basemods_gff = pb_basemods_hifi.basemods_gff
    File? basemods_csv = pb_basemods_hifi.basemods_csv
    File? report_modifications = pb_basemods_hifi.report_modifications
    File? bigwig_file = pb_basemods_hifi.bigwig_file
    File? motifs_csv = pb_basemods_hifi.motifs_csv
    File? motifs_gff = pb_basemods_hifi.motifs_gff
    File? report_motifs = pb_basemods_hifi.report_motifs
  }
}
