# Combined CLR microbial assembly and basemods, for testing against HiFi
# implementation

version 1.0

import "pb_basemods.wdl"
import "pb_assembly_microbial.wdl"
import "pb_microbial_analysis.wdl"

workflow dev_clr_assembly_basemods {
  input {
    File eid_subread

    Int microasm_length_cutoff = -1
    Int microasm_downsampled_coverage = 100
    String microasm_genome_size = "5M"
    Int microasm_coverage = 30
    String microasm_advanced_options = ""
    Int microasm_plasmid_contig_len_max = 300000
    String? mapping_biosample_name

    Boolean run_find_motifs = true
    String kineticstools_identify_mods = "m4C,m6A"
    Float kineticstools_p_value = 0.001
    # TODO this cutoff needs to be properly validated
    Int motif_min_score = 20
    Float motif_min_fraction = 0.30

    Int nproc = 1
    String log_level = "INFO"
    Int max_nchunks = 16
    String tmp_dir = "/tmp"
  }

  call pb_assembly_microbial.pb_assembly_microbial {
    input:
      eid_subread = eid_subread,
      microasm_length_cutoff = microasm_length_cutoff,
      microasm_downsampled_coverage = microasm_downsampled_coverage,
      microasm_genome_size = microasm_genome_size,
      microasm_coverage = microasm_coverage,
      microasm_advanced_options = microasm_advanced_options,
      microasm_plasmid_contig_len_max = microasm_plasmid_contig_len_max,
      mapping_biosample_name = mapping_biosample_name,
      nproc = nproc,
      max_nchunks = max_nchunks,
      tmp_dir = tmp_dir,
      log_level = log_level
  }

  call pb_microbial_analysis.make_reference {
    input:
      assembly_fasta = pb_assembly_microbial.assembled_fasta
  }

  call pb_basemods.pb_basemods {
    input:
      eid_subread = eid_subread,
      eid_ref_dataset = make_reference.reference,
      run_find_motifs = run_find_motifs,
      kineticstools_identify_mods = kineticstools_identify_mods,
      kineticstools_p_value = kineticstools_p_value,
      motif_min_score = motif_min_score,
      motif_min_fraction = motif_min_fraction,
      nproc = nproc,
      max_nchunks = max_nchunks,
      tmp_dir = tmp_dir,
      log_level = log_level
  }

  output {
    File assembled_fasta = pb_assembly_microbial.assembled_fasta
    File mapped = pb_basemods.mapped
    File report_polished_assembly = pb_assembly_microbial.report_polished_assembly
    File? report_mapping_stats = pb_basemods.report_mapping_stats
    File coverage_gff = pb_basemods.coverage_gff
    File report_coverage = pb_basemods.report_coverage
    File basemods_gff = pb_basemods.basemods_gff
    File? basemods_csv = pb_basemods.basemods_csv
    File? report_modifications = pb_basemods.report_modifications
    File? bigwig_file = pb_basemods.bigwig_file
    File? motifs_csv = pb_basemods.motifs_csv
    File? motifs_gff = pb_basemods.motifs_gff
    File? report_motifs = pb_basemods.report_motifs
    File? mapped_bam_datastore = pb_basemods.mapped_bam_datastore
  }
}
