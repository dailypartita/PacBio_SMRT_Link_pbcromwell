# IPA HiFi Assembler.
# Authors: Ivan Sovic, Zev Kronenberg, Christopher Dunn, Derek Barnett

version 1.0

import "tasks/pbcoretools.wdl" as pbcoretools
import "tasks/pbreports.wdl" as pbreports
import "wf_prepare_input.wdl" as wf_prepare_input

task generate_config {
  input {
    String advanced_opt_str
    Int downsampled_coverage
    String genome_size
    Boolean run_polishing
    Boolean run_phasing
    Boolean run_purge_dups = true
    Float m4filt_high_copy_sample_rate = 1.0
    String purge_dups_calcuts = ""
    String purge_dups_get_seqs = ""
    String log_level
    String tmp_dir
  }
  command {
    params_advanced_opt="${advanced_opt_str}" \
    params_coverage="${downsampled_coverage}" \
    params_genome_size="${genome_size}" \
    params_polish_run=${true="1" false="0" run_polishing} \
    params_phase_run=${true="1" false="0" run_phasing} \
    params_purge_dups_run=${true="1" false="0" run_purge_dups} \
    params_purge_dups_calcuts="${purge_dups_calcuts}" \
    params_purge_dups_get_seqs="${purge_dups_get_seqs}" \
    params_m4filt_high_copy_sample_rate="${m4filt_high_copy_sample_rate}" \
    params_log_level="${log_level}" \
    params_tmp_dir="${tmp_dir}" \
    output_fn="generated.config.sh" \
    sentinel_fn="generated.config" \
      ipa2-task generate_config_from_workflow
  }
  runtime {
    cpu: 1
    memory: "1024MB"
  }
  output {
    File config_sh = "generated.config"
  }
}

task build_db {
  input {
    File reads_fn
    File config_sh_fn
    String db_prefix
    Int num_threads
    String log_level
    String tmp_dir
    Int mem_scale_factor = 8
    Int base_memory_mb = 0
  }
  command <<<
    set -e
    input_reads_fn="`readlink -f ~{reads_fn}`" \
    params_db_prefix="~{db_prefix}" \
    params_config_sh_fn="~{config_sh_fn}" \
    params_num_threads=~{num_threads} \
    params_log_level="~{log_level}" \
    params_tmp_dir="~{tmp_dir}" \
      ipa2-task build_db
    echo `pwd`/~{db_prefix}.seqdb > seqdb_fn
    echo `pwd`/~{db_prefix}.seeddb > seeddb_fn
    echo `pwd`/~{db_prefix}.seqdb.0.seq > seqdb_seqs_fn
    echo `pwd`/~{db_prefix}.seeddb.0.seeds > seeddb_seeds_fn
  >>>
  Int mem_total_mb = base_memory_mb + 256 + (num_threads * 256 * mem_scale_factor)
  runtime {
    cpu: num_threads
    memory: "${mem_total_mb}MB"
  }
  output {
    String seqdb_fn = read_string("seqdb_fn")
    String seeddb_fn = read_string("seeddb_fn")
    String seqdb_seqs_fn = read_string("seqdb_seqs_fn")
    String seeddb_seeds_fn = read_string("seeddb_seeds_fn")
    File input_fofn = "input.fofn"
  }
}

task ovl_asym_prepare {
  input {
    String seqdb_fn
    Int max_nchunks
    String log_level
    String tmp_dir
  }
  command {
    input_db="${seqdb_fn}" \
    output_shard_ids="./shard_ids" \
    output_pwd="./pwd.txt" \
    params_max_nchunks="${max_nchunks}" \
    params_log_level="${log_level}" \
    params_tmp_dir="${tmp_dir}" \
      ipa2-task ovl_asym_prepare
  }
  runtime {
    cpu: 1
    memory: "1024MB"
  }
  output {
    File shard_ids = "shard_ids"
    #String shard_ids_dn = read_string("pwd.txt")
    File shard_ids_pwd = "pwd.txt"
  }
}

task ovl_asym_run {
  input {
    String seqdb_fn
    String seeddb_fn
    String seqdb_seqs_fn
    String seeddb_seeds_fn
    File config_sh_fn
    File shard_ids_fn # not needed, but shows the dependency-graph
    String shard_ids_dn # for the directory-name
    String shard_id
    String db_prefix
    Int num_threads
    String log_level
    String tmp_dir
    Int mem_scale_factor = 8
    Int base_memory_mb = 0
  }
  command <<<
    set -e
    input_seqdb="~{seqdb_fn}" \
    params_ovl_asym_prepare_dn="~{shard_ids_dn}" \
    params_shard_id=~{shard_id} \
    params_db_prefix=~{db_prefix} \
    params_num_threads=~{num_threads} \
    params_config_sh_fn="~{config_sh_fn}" \
    params_log_level="~{log_level}" \
    params_tmp_dir="~{tmp_dir}" \
      ipa2-task ovl_asym_run
    echo `pwd`/ovl.m4 > out_m4_fn
    echo `pwd`/ovl.sorted.m4 > out_sorted_m4_fn
  >>>
  Int total_mem_mb = base_memory_mb + 8192 * mem_scale_factor
  runtime {
    cpu: num_threads
    memory: "${total_mem_mb}MB"
  }
  output {
    String out_m4_fn = read_string("out_m4_fn")
    String out_sorted_m4_fn = read_string("out_sorted_m4_fn")
  }
}

task ovl_asym_merge {
  input {
    Array[String] in_fns
    File config_sh_fn
    Int base_memory_mb = 0
    Int num_threads
    String db_prefix
    String log_level
    String tmp_dir
  }
  command <<<
    set -vex

    echo ~{sep=' ' in_fns} | xargs -n 1 > merged.fofn

    input_fofn=./merged.fofn \
    params_num_threads="~{num_threads}" \
    params_db_prefix="~{db_prefix}" \
    params_config_sh_fn="~{config_sh_fn}" \
    params_log_level="~{log_level}" \
    params_tmp_dir="~{tmp_dir}" \
      ipa2-task ovl_asym_merge
    echo `pwd`/ovl.merged.m4 > m4_merged_raw_fn
    echo `pwd`/ovl.nonlocal.m4 > m4_filtered_nonlocal_fn
  >>>
  Int total_mem_mb = base_memory_mb + 1024
  runtime {
    cpu: num_threads
    memory: "${total_mem_mb}MB"
  }
  output {
    String m4_merged_raw_fn = read_string("m4_merged_raw_fn")
    String m4_filtered_nonlocal_fn = read_string("m4_filtered_nonlocal_fn")
  }
}

#################
### Phasing #####
#################
task phasing_prepare {
  input {
    String seqdb_fn
    String m4_fn
    File config_sh_fn
    Int max_nchunks
    String log_level
    String tmp_dir
    Int mem_scale_factor = 8
    Int base_memory_mb = 0
  }
  command {
    output_shard_ids="./shard_ids" \
    output_pwd="./pwd.txt" \
    input_m4="${m4_fn}" \
    params_config_sh_fn="${config_sh_fn}" \
    params_max_nchunks="${max_nchunks}" \
    params_log_level="${log_level}" \
    params_tmp_dir="${tmp_dir}" \
      ipa2-task phasing_prepare
  }
  # XXX this is overkill - it usually needs <2GB but it can bloat sometimes
  Int total_mem_mb = base_memory_mb + 2048 + (512 * mem_scale_factor)
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File shard_ids = "shard_ids"
    #String shard_ids_dn = read_string("pwd.txt")
    File shard_ids_pwd = "pwd.txt"
  }
}

task phasing_run {
  input {
    String seqdb_fn
    String seqdb_seqs_fn
    File config_sh_fn
    File shard_ids_fn # not needed, but shows the dependency-graph
    String shard_ids_dn # for the directory-name
    String shard_id
    Int num_threads
    String log_level
    String tmp_dir
    Int mem_scale_factor = 8
    Int base_memory_mb = 0
  }
  command <<<
    set -e
    input_seqdb="~{seqdb_fn}" \
    output_keep_m4="ovl.phased.m4" \
    output_scraps_m4="ovl.phased.m4.scraps" \
    output_outdir_fn="outdir_fn.txt" \
    params_phasing_prepare_dn="~{shard_ids_dn}" \
    params_shard_id=~{shard_id} \
    params_num_threads=~{num_threads} \
    params_config_sh_fn="~{config_sh_fn}" \
    params_log_level="~{log_level}" \
    params_tmp_dir="~{tmp_dir}" \
      ipa2-task phasing_run
    echo `pwd`/ovl.phased.m4 > out_keep_m4_fn
    echo `pwd`/ovl.phased.m4.scraps > out_scraps_m4_fn
  >>>
  Int total_mem_mb = base_memory_mb + 256 + (num_threads * 256 * mem_scale_factor)
  runtime {
    cpu: num_threads
    memory: "${total_mem_mb}MB"
  }
  output {
    String out_keep_m4_fn = read_string("out_keep_m4_fn")
    String out_scraps_m4_fn = read_string("out_scraps_m4_fn")
    # Workaround to allow multiple files to be gathered in the next task.
    File outdir_fn = "outdir_fn.txt"
  }
}

task phasing_merge {
  input {
    Array[String] in_fns
    String original_m4_fn
    File config_sh_fn
    Int num_threads
    String log_level
    String tmp_dir
    Int mem_scale_factor = 8
    Int base_memory_mb = 0
  }
  command <<<
    set -vex

    # Slight difference from the Snakemake workflow in order to get multiple input files.
    # Instead of in_fns pointing to a particular of either of the "ovl.phased.m4" or "ovl.phased.m4.scraps",
    # we are using a sentinel file "outdir_fn.txt" which holds the pwd of the phasing_run folder which
    # generated these files. This is required because Cromwell changes the folder for the input files,
    # and we need to get to the original folder in order to be able to collect the FOFNs.
    echo ~{sep=' ' in_fns} | xargs -n 1 cat | awk '{ print $1"/ovl.phased.m4" }' > merged.keep.fofn
    echo ~{sep=' ' in_fns} | xargs -n 1 cat | awk '{ print $1"/ovl.phased.m4.scraps" }' > merged.scraps.fofn
    # Note: Empty .m4 is legal and can occur when a run was a no-op.

    input_keep_fofn="./merged.keep.fofn" \
    input_scraps_fofn="./merged.scraps.fofn" \
    input_original_m4="~{original_m4_fn}" \
    output_m4="ovl.phased.m4" \
    params_num_threads=~{num_threads} \
    params_config_sh_fn="~{config_sh_fn}" \
    params_log_level="~{log_level}" \
    params_tmp_dir="~{tmp_dir}" \
      ipa2-task phasing_merge
    echo `pwd`/ovl.phased.m4 > gathered_m4_fn
    echo `pwd`/all.keep.m4 > all_keeps_m4_fn
    echo `pwd`/all.scraps.m4 > all_scraps_m4_fn
  >>>
  # FIXME should this depend on num_threads?
  Int total_mem_mb = base_memory_mb + 4096 + 2048 * mem_scale_factor
  runtime {
    cpu: num_threads
    memory: "${total_mem_mb}MB"
  }
  output {
    String gathered_m4_fn = read_string("gathered_m4_fn")
    String all_keeps_m4_fn = read_string("all_keeps_m4_fn")
    String all_scraps_m4_fn = read_string("all_scraps_m4_fn")
  }
}

task ovl_filter {
  input {
    String m4_fn
    File config_sh_fn
    Int num_threads
    String log_level
    String tmp_dir
    Int mem_scale_factor = 8
    Int base_memory_mb = 0
  }
  command {
    set -vex

    input_m4="${m4_fn}" \
    output_m4_final="ovl.final.m4" \
    output_m4_chimerfilt="ovl.chimerfilt.m4" \
    params_num_threads=${num_threads} \
    params_config_sh_fn="${config_sh_fn}" \
    params_log_level="${log_level}" \
    params_tmp_dir="${tmp_dir}" \
      ipa2-task ovl_filter
    echo `pwd`/ovl.final.m4 > m4_final_fn
    echo `pwd`/ovl.chimerfilt.m4 > m4_chimerfilt_fn
  }
  # FIXME this is probably an extreme upper bounds; num_threads dependent?
  Int total_mem_mb = base_memory_mb + 7 * 1024 * mem_scale_factor
  runtime {
    cpu: num_threads
    memory: "${total_mem_mb}MB"
  }
  output {
    String m4_final_fn = read_string("m4_final_fn")
    String m4_chimerfilt_fn = read_string("m4_chimerfilt_fn")
  }
}

task assemble {
  input {
    File reads
    String seqdb_fn
    String seqdb_seqs_fn
    String m4_fn
    String m4_phasing_merge_fn
    File config_sh_fn
    String ctg_prefix
    Int num_threads
    String log_level
    String tmp_dir
    Int mem_scale_factor = 8
    Int base_memory_mb = 0
  }

  command {
    output_read_to_contig="./read_to_contig.paf" \
    input_seqdb="${seqdb_fn}" \
    input_m4="${m4_fn}" \
    input_m4_phasing_merge="${m4_phasing_merge_fn}" \
    input_reads="`readlink -f ${reads}`" \
    params_ctg_prefix="${ctg_prefix}" \
    params_num_threads=${num_threads} \
    params_config_sh_fn="${config_sh_fn}" \
    params_log_level="${log_level}" \
    params_tmp_dir="${tmp_dir}" \
      ipa2-task assemble
    echo `pwd`/p_ctg.fasta > p_ctg_fasta_fn
    echo `pwd`/a_ctg.fasta > a_ctg_fasta_fn
  }
  # FIXME this is probably an extreme upper bounds; num_threads dependent?
  Int total_mem_mb = base_memory_mb + 8192 * mem_scale_factor
  runtime {
    cpu: num_threads
    memory: "${total_mem_mb}MB"
  }
  output {
    String p_ctg_fasta_fn = read_string("p_ctg_fasta_fn")
    String a_ctg_fasta_fn = read_string("a_ctg_fasta_fn")
    File p_ctg_tiling_path = "p_ctg_tiling_path"
    File a_ctg_tiling_path = "a_ctg_tiling_path"
    File p_ctg_fasta_fai = "p_ctg.fasta.fai"
    File a_ctg_fasta_fai = "a_ctg.fasta.fai"
    File read_to_contig = "read_to_contig.paf"
    # File asm_gfa = "asm.gfa"
    # File sg_gfa = "sg.gfa"
    # File contig_gfa2 = "contig.gfa2"
    File circular_contigs = "circular_contigs.csv"
  }
}

#################
### Polishing ###
#################
task polish_prepare {
  input {
    File read_to_contig
    File p_ctg_fasta
    File a_ctg_fasta
    File p_ctg_fasta_fai
    File a_ctg_fasta_fai
    File config_sh_fn
    Int max_nchunks
    String log_level
    String tmp_dir
    Int mem_scale_factor = 8
    Int base_memory_mb = 0
  }
  command {
    output_shard_ids="./shard_ids" \
    output_pwd="./pwd.txt" \
    input_read_to_contig="${read_to_contig}" \
    input_p_ctg_fasta="${p_ctg_fasta}" \
    input_a_ctg_fasta="${a_ctg_fasta}" \
    input_p_ctg_fasta_fai="${p_ctg_fasta_fai}" \
    input_a_ctg_fasta_fai="${a_ctg_fasta_fai}" \
    params_config_sh_fn="${config_sh_fn}" \
    params_max_nchunks="${max_nchunks}" \
    params_log_level="${log_level}" \
    params_tmp_dir="${tmp_dir}" \
      ipa2-task polish_prepare
  }
  # FIXME should this depend on num_threads?
  Int total_mem_mb = base_memory_mb + 2048 + (1280 * mem_scale_factor)
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File shard_ids = "shard_ids"
    #String shard_ids_dn = read_string("pwd.txt")
    File shard_ids_pwd = "pwd.txt"
  }
}

task polish_run {
  input {
    File input_fofn
    String seqdb_fn
    String seqdb_seqs_fn
    File p_ctg_fasta
    File a_ctg_fasta
    File p_ctg_fasta_fai
    File a_ctg_fasta_fai
    File config_sh_fn
    File shard_ids_fn # not needed, but shows the dependency-graph
    String shard_ids_dn # for the directory-name
    String shard_id
    Int num_threads
    String log_level
    String tmp_dir
    Int mem_scale_factor = 8
    Int base_memory_mb = 0
  }
  command <<<
    set -e
    input_fofn="~{input_fofn}" \
    input_seqdb="~{seqdb_fn}" \
    input_p_ctg_fasta="~{p_ctg_fasta}" \
    input_a_ctg_fasta="~{a_ctg_fasta}" \
    output_consensus="consensus.fasta" \
    params_polish_prepare_dn="~{shard_ids_dn}" \
    params_shard_id=~{shard_id} \
    params_num_threads=~{num_threads} \
    params_config_sh_fn="~{config_sh_fn}" \
    params_log_level="~{log_level}" \
    params_tmp_dir="~{tmp_dir}" \
      ipa2-task polish_run
    echo `pwd`/consensus.fasta > consensus_fn
    echo `pwd` > dirname.txt
  >>>
  # FIXME should this depend on num_threads?
  Int total_mem_mb = base_memory_mb + 2048 * mem_scale_factor
  runtime {
    cpu: num_threads
    memory: "${total_mem_mb}MB"
  }
  output {
    String consensus_fn = read_string("consensus_fn")
    String dn = read_string("dirname.txt")
  }
}

task polish_merge {
  input {
    Array[String] in_fns
    Array[String] in_dns
    File p_ctg_fasta
    File a_ctg_fasta
    File config_sh_fn
    Int base_memory_mb = 0
    Int num_threads = 1
    String log_level
    String tmp_dir
  }
  command {
    set -vex

    echo ${sep=' ' in_fns} | xargs -n 1 > merged.fofn
    echo ${sep=' ' in_dns} | xargs -n 1 > merged.fodn

    input_fofn=./merged.fofn \
    input_fodn=./merged.fodn \
    input_p_ctg_fasta="${p_ctg_fasta}" \
    input_a_ctg_fasta="${a_ctg_fasta}" \
    output_fasta="assembly.merged.fasta" \
    output_p_paf="p.paf" \
    output_a_paf="a.paf" \
    params_num_threads=${num_threads} \
    params_config_sh_fn="${config_sh_fn}" \
    params_log_level="${log_level}" \
    params_tmp_dir="${tmp_dir}" \
      ipa2-task polish_merge
    echo `pwd`/assembly.merged.fasta > consensus_merged_fn
  }
  runtime {
    cpu: 1
    memory: "1024MB"
  }
  output {
    String consensus_merged_fn = read_string("consensus_merged_fn")
    File consensus_merged_fai = "assembly.merged.fasta.fai"
    File p_paf = "p.paf"
    File a_paf = "a.paf"
  }
}

task separate_p_from_a {
  input {
    String assembly_merged_fasta_fn
    File assembly_merged_fai
    File assembly_p_ctg_fasta
    File assembly_a_ctg_fasta
    Int num_threads  # unused
    String log_level
    String tmp_dir
  }
  command {
    set -vex

    input_assembly_merged_fasta="${assembly_merged_fasta_fn}" \
    input_p_ctg_fasta="${assembly_p_ctg_fasta}" \
    input_a_ctg_fasta="${assembly_a_ctg_fasta}" \
    params_log_level="${log_level}" \
    params_tmp_dir="${tmp_dir}" \
    params_polish_run=1 \
    output_p_ctg_fasta="p_ctg.fasta" \
    output_a_ctg_fasta="a_ctg.fasta" \
        ipa2-task separate_p_from_a
    echo `pwd`/p_ctg.fasta > p_ctg_fasta_fn
    echo `pwd`/a_ctg.fasta > a_ctg_fasta_fn
  }
  runtime {
    cpu: num_threads
    memory: "4096MB"
  }
  output {
    String p_ctg_fasta_fn = read_string("p_ctg_fasta_fn")
    String a_ctg_fasta_fn = read_string("a_ctg_fasta_fn")
  }
}

##################
### Purge dups ###
##################
task purge_dups_map_prepare {
  input {
    String seqdb_fn
    String p_ctg_fasta
    File config_sh_fn
    Int num_threads
    Int max_nchunks
    String log_level
    String tmp_dir
    Int mem_scale_factor = 8
    Int base_memory_mb = 0
  }
  command <<<
    output_shard_ids="./shard_ids" \
    output_pwd="./pwd.txt" \
    input_reads_db="~{seqdb_fn}" \
    input_primary_fasta="~{p_ctg_fasta}" \
    params_config_sh_fn="~{config_sh_fn}" \
    params_num_threads="~{num_threads}" \
    params_max_nchunks="~{max_nchunks}" \
    params_log_level="~{log_level}" \
        ipa2-task purge_dups_map_prepare
  >>>
  Int total_mem_mb = base_memory_mb + 512 + (num_threads * 256 * mem_scale_factor)
  runtime {
    cpu: num_threads
    memory: "${total_mem_mb}MB"
  }
  output {
    File shard_ids = "shard_ids"
    #String shard_ids_dn = read_string("pwd.txt")
    File shard_ids_pwd = "pwd.txt"
  }
}

task purge_dups_map_run {
  input {
    String seqdb_fn
    File config_sh_fn
    File shard_ids_fn # not needed, but shows the dependency-graph
    String shard_ids_dn # for the directory-name
    String shard_id
    Int num_threads
    String log_level
    String tmp_dir = "/tmp"  # ignored for now
    Int mem_scale_factor = 8
    Int base_memory_mb = 0
  }
  command <<<
    set -vex

    input_seqdb="~{seqdb_fn}" \
    params_purge_dups_map_prepare_dn="~{shard_ids_dn}" \
    params_shard_id=~{shard_id} \
    params_num_threads=~{num_threads} \
    params_config_sh_fn="~{config_sh_fn}" \
    params_log_level="~{log_level}" \
    output_paf="out.paf" \
        ipa2-task purge_dups_map_run
    echo `pwd`/out.paf > paf_fn.txt
    echo `pwd` > dirname.txt
  >>>
  # FIXME should this depend on num_threads? or input size?
  Int total_mem_mb = base_memory_mb + 4096 * mem_scale_factor
  runtime {
    cpu: num_threads
    memory: "${total_mem_mb}MB"
  }
  output {
    String paf_fn = read_string("paf_fn.txt")
    String dn = read_string("dirname.txt")
  }
}

task purge_dups_map_merge {
  input {
    Array[String] in_fns
    Array[String] in_dns #unused
    File config_sh_fn
    Int num_threads = 1 #unused
    String log_level #unused
    String tmp_dir #unused
  }
  command {
    set -vex

    echo ${sep=' ' in_fns} | xargs -n 1 > merged.fofn
    echo ${sep=' ' in_dns} | xargs -n 1 > merged.fodn

    input_fofn=./merged.fofn \
    input_fodn=./merged.fodn \
    output_paf="merged.paf" \
    params_config_sh_fn="${config_sh_fn}" \
      ipa2-task purge_dups_map_merge
    echo `pwd`/merged.paf > merged_paf_fn.txt
  }
  runtime {
    cpu: 1
    memory: "1024MB"
  }
  output {
    String merged_paf_fn = read_string("merged_paf_fn.txt")
  }
}

task purge_dups {
  input {
    String unpurged_p_ctg_fasta
    String unpurged_a_ctg_fasta
    String merged_paf_fn
    File config_sh_fn
    Int num_threads = 1
    String log_level
    Int mem_scale_factor = 8
    Int base_memory_mb = 0
  }
  command <<<
    set -vex

    input_paf="~{merged_paf_fn}" \
    input_primary_fasta="~{unpurged_p_ctg_fasta}" \
    input_haplotigs_fasta="~{unpurged_a_ctg_fasta}" \
    params_config_sh_fn="~{config_sh_fn}" \
    params_num_threads="~{num_threads}" \
    params_log_level="~{log_level}" \
    output_primary_fasta="final_purged_primary.fasta" \
    output_haplotigs_fasta="final_purged_haplotigs.fasta" \
        ipa2-task purge_dups_paf
    echo `pwd` > dirname.txt
  >>>
  # FIXME should this depend on num_threads? or input size?
  Int total_mem_mb = base_memory_mb + 4096 * mem_scale_factor
  runtime {
    cpu: num_threads
    memory: "${total_mem_mb}MB"
  }
  output {
    File p_ctg_fasta = "final_purged_primary.fasta"
    File a_ctg_fasta = "final_purged_haplotigs.fasta"
    String dn = read_string("dirname.txt")
  }
}

task cleanup_files {
  input {
    Array[String] singular_fns
    Array[String] ovl_run_m4_fns
    Array[String] ovl_run_sorted_m4_fns
    Array[String] phasing_run_keep_fns
    Array[String] phasing_run_scraps_fns
    Array[String] polish_run_consensus_fns
    Array[String] purge_dups_map_run_fns
    String phasing_prepare_dn
    File purge_dups_fn # strictly for task-ordering
    Int num_threads  # unused
  }
  command <<<
    set -vex

    echo ~{sep=' ' singular_fns} | xargs -n 1 > removed_files.fofn
    echo ~{sep=' ' ovl_run_m4_fns} | xargs -n 1 >> removed_files.fofn
    echo ~{sep=' ' ovl_run_sorted_m4_fns} | xargs -n 1 >> removed_files.fofn
    echo ~{sep=' ' phasing_run_keep_fns} | xargs -n 1 >> removed_files.fofn
    echo ~{sep=' ' phasing_run_scraps_fns} | xargs -n 1 >> removed_files.fofn
    echo ~{sep=' ' polish_run_consensus_fns} | xargs -n 1 >> removed_files.fofn
    echo ~{sep=' ' purge_dups_map_run_fns} | xargs -n 1 >> removed_files.fofn

    for fn in $(find ~{phasing_prepare_dn} -name "chunk.*.m4")
    do
        echo ${fn} >> removed_files.fofn
    done

    input_fofn=./removed_files.fofn \
        ipa2-task cleanup_files
  >>>
  runtime {
    cpu: 1
    memory: "1024MB"
  }
  output {
    File removed_files_fofn = "removed_files.fofn"
  }
}

task report_assembly2 {
  input {
    File contigs_fasta
    File haplotigs_fasta
    File? circular_contigs_txt
    Int mem_scale_factor = 8
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 \
      -m pbreports.report.assembly2 \
      --log-level ${log_level} \
      ${"--circular-contigs " + circular_contigs_txt} \
      ${contigs_fasta} ${haplotigs_fasta} \
      polished_assembly.report.json
  }
  Int total_mem_mb = base_memory_mb + 1024 * mem_scale_factor
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File report = "polished_assembly.report.json"
  }
}

workflow pb_assembly_hifi {
  input {
    File? eid_ccs
    File? reads
    String ipa2_genome_size = "0k"
    Int ipa2_downsampled_coverage = 0
    String ipa2_advanced_options = ""
    Boolean ipa2_run_polishing = true
    Boolean ipa2_run_phasing = true
    Boolean ipa2_run_purge_dups = true
    String ipa2_ctg_prefix = "ctg"
    String ipa2_reads_db_prefix = "reads"  # this will soon be implicit
    Boolean ipa2_cleanup_intermediate_files = true

    String dataset_filters = ""
    Int filter_min_qv = 20
    Int downsample_factor = 0
    Int mem_scale_factor = 8
    Int add_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
    Int max_nchunks = 40
    String tmp_dir = "/tmp"
  }

  call wf_prepare_input.prepare_input {
    input:
      dataset_xml = eid_ccs,
      reads = reads,
      dataset_filters = dataset_filters,
      downsample_factor = downsample_factor,
      filter_min_qv = filter_min_qv,
      nproc = 1,
      log_level = log_level
  }

  call generate_config {
    input:
      advanced_opt_str = ipa2_advanced_options,
      downsampled_coverage = ipa2_downsampled_coverage,
      genome_size = ipa2_genome_size,
      run_polishing = ipa2_run_polishing,
      run_phasing = ipa2_run_phasing,
      run_purge_dups = ipa2_run_purge_dups,
      log_level = log_level,
      tmp_dir = tmp_dir,
  }

  call build_db {
    input:
      reads_fn = prepare_input.reads_file,
      db_prefix = ipa2_reads_db_prefix,
      config_sh_fn = generate_config.config_sh,
      num_threads = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir,
      mem_scale_factor = mem_scale_factor,
      base_memory_mb = add_memory_mb
  }

  call ovl_asym_prepare {
    input:
      seqdb_fn = build_db.seqdb_fn,
      max_nchunks = max_nchunks,
      log_level = log_level,
      tmp_dir = tmp_dir,
  }

  Array[String] shard_ids_ovl = read_lines(ovl_asym_prepare.shard_ids)
  String shard_ids_ovl_dn = read_string(ovl_asym_prepare.shard_ids_pwd)

  scatter (shard_id in shard_ids_ovl) {
    call ovl_asym_run {
      input:
        seqdb_fn = build_db.seqdb_fn,
        seeddb_fn = build_db.seeddb_fn,
        seqdb_seqs_fn = build_db.seqdb_seqs_fn,
        seeddb_seeds_fn = build_db.seeddb_seeds_fn,
        config_sh_fn = generate_config.config_sh,
        shard_ids_fn = ovl_asym_prepare.shard_ids,
        shard_ids_dn = shard_ids_ovl_dn,
        shard_id = shard_id,
        db_prefix = ipa2_reads_db_prefix,
        num_threads = nproc,
        log_level = log_level,
        tmp_dir = tmp_dir,
        mem_scale_factor = mem_scale_factor,
        base_memory_mb = add_memory_mb
    }
  }

  call ovl_asym_merge {
    input:
      in_fns = ovl_asym_run.out_sorted_m4_fn,
      config_sh_fn = generate_config.config_sh,
      num_threads = nproc,
      db_prefix = ipa2_reads_db_prefix,
      base_memory_mb = add_memory_mb,
      log_level = log_level,
      tmp_dir = tmp_dir,
  }

  call phasing_prepare {
    input:
      seqdb_fn = build_db.seqdb_fn,
      m4_fn = ovl_asym_merge.m4_filtered_nonlocal_fn,
      config_sh_fn = generate_config.config_sh,
      max_nchunks = max_nchunks,
      log_level = log_level,
      tmp_dir = tmp_dir,
      mem_scale_factor = mem_scale_factor,
      base_memory_mb = add_memory_mb,
  }

  Array[String] shard_ids_phasing = read_lines(phasing_prepare.shard_ids)
  String shard_ids_phasing_dn = read_string(phasing_prepare.shard_ids_pwd)

  scatter (shard_id in shard_ids_phasing) {
    call phasing_run {
      input:
        seqdb_fn = build_db.seqdb_fn,
        seqdb_seqs_fn = build_db.seqdb_seqs_fn,
        config_sh_fn = generate_config.config_sh,
        shard_ids_fn = phasing_prepare.shard_ids,
        shard_ids_dn = shard_ids_phasing_dn,
        shard_id = shard_id,
        num_threads = nproc,
        log_level = log_level,
        tmp_dir = tmp_dir,
        mem_scale_factor = mem_scale_factor,
        base_memory_mb = add_memory_mb,
    }
  }

  call phasing_merge {
    input:
      in_fns = phasing_run.outdir_fn,
      original_m4_fn = ovl_asym_merge.m4_filtered_nonlocal_fn,
      config_sh_fn = generate_config.config_sh,
      num_threads = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir,
      mem_scale_factor = mem_scale_factor,
      base_memory_mb = add_memory_mb,
  }

  call ovl_filter{
    input:
      m4_fn = phasing_merge.gathered_m4_fn,
      config_sh_fn = generate_config.config_sh,
      num_threads = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir,
      mem_scale_factor = mem_scale_factor,
      base_memory_mb = add_memory_mb,
  }

  call assemble {
    input:
      reads = prepare_input.reads_file,
      seqdb_fn = build_db.seqdb_fn,
      seqdb_seqs_fn = build_db.seqdb_seqs_fn,
      m4_fn = ovl_filter.m4_final_fn,
      m4_phasing_merge_fn = phasing_merge.gathered_m4_fn,
      config_sh_fn = generate_config.config_sh,
      ctg_prefix = ipa2_ctg_prefix,
      num_threads = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir,
      mem_scale_factor = mem_scale_factor,
      base_memory_mb = add_memory_mb,
  }

  call polish_prepare {
    input:
      read_to_contig = assemble.read_to_contig,
      p_ctg_fasta = assemble.p_ctg_fasta_fn,
      a_ctg_fasta = assemble.a_ctg_fasta_fn,
      p_ctg_fasta_fai = assemble.p_ctg_fasta_fai,
      a_ctg_fasta_fai = assemble.a_ctg_fasta_fai,
      config_sh_fn = generate_config.config_sh,
      max_nchunks = max_nchunks,
      log_level = log_level,
      tmp_dir = tmp_dir,
      mem_scale_factor = mem_scale_factor,
      base_memory_mb = add_memory_mb,
  }

  Array[String] shard_ids_polish = read_lines(polish_prepare.shard_ids)
  String shard_ids_polish_dn = read_string(polish_prepare.shard_ids_pwd)

  scatter (shard_id in shard_ids_polish) {
    call polish_run {
      input:
        input_fofn = build_db.input_fofn,
        seqdb_fn = build_db.seqdb_fn,
        seqdb_seqs_fn = build_db.seqdb_seqs_fn,
        p_ctg_fasta = assemble.p_ctg_fasta_fn,
        a_ctg_fasta = assemble.a_ctg_fasta_fn,
        p_ctg_fasta_fai = assemble.p_ctg_fasta_fai,
        a_ctg_fasta_fai = assemble.a_ctg_fasta_fai,
        config_sh_fn = generate_config.config_sh,
        shard_ids_fn = polish_prepare.shard_ids,
        shard_ids_dn = shard_ids_polish_dn,
        shard_id = shard_id,
        num_threads = nproc,
        log_level = log_level,
        tmp_dir = tmp_dir,
        mem_scale_factor = mem_scale_factor,
      base_memory_mb = add_memory_mb,
    }
  }

  call polish_merge {
    input:
      in_fns = polish_run.consensus_fn,
      in_dns = polish_run.dn,
      p_ctg_fasta = assemble.p_ctg_fasta_fn,
      a_ctg_fasta = assemble.a_ctg_fasta_fn,
      config_sh_fn = generate_config.config_sh,
      num_threads = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir,
  }

  call separate_p_from_a {
    input:
      assembly_merged_fasta_fn = polish_merge.consensus_merged_fn,
      assembly_merged_fai = polish_merge.consensus_merged_fai,
      assembly_p_ctg_fasta = assemble.p_ctg_fasta_fn,
      assembly_a_ctg_fasta = assemble.a_ctg_fasta_fn,
      num_threads = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir,
  }

  call purge_dups_map_prepare {
    input:
      seqdb_fn = build_db.seqdb_fn,
      p_ctg_fasta = separate_p_from_a.p_ctg_fasta_fn,
      config_sh_fn = generate_config.config_sh,
      num_threads = nproc,
      max_nchunks = max_nchunks,
      log_level = log_level,
      tmp_dir = tmp_dir,
      mem_scale_factor = mem_scale_factor,
      base_memory_mb = add_memory_mb,
  }

  Array[String] shard_ids_purge_dups_map = read_lines(purge_dups_map_prepare.shard_ids)
  String shard_ids_purge_dups_map_dn = read_string(purge_dups_map_prepare.shard_ids_pwd)

  scatter (shard_id in shard_ids_purge_dups_map) {
    call purge_dups_map_run {
      input:
        seqdb_fn = build_db.seqdb_fn,
        shard_ids_fn = purge_dups_map_prepare.shard_ids,
        shard_ids_dn = shard_ids_purge_dups_map_dn,
        shard_id = shard_id,
        config_sh_fn = generate_config.config_sh,
        num_threads = nproc,
        log_level = log_level,
        tmp_dir = tmp_dir,
        mem_scale_factor = mem_scale_factor,
        base_memory_mb = add_memory_mb,
    }
  }

  call purge_dups_map_merge {
    input:
      in_fns = purge_dups_map_run.paf_fn,
      in_dns = purge_dups_map_run.dn,
      config_sh_fn = generate_config.config_sh,
      num_threads = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir,
  }

  call purge_dups {
    input:
      unpurged_p_ctg_fasta = separate_p_from_a.p_ctg_fasta_fn,
      unpurged_a_ctg_fasta = separate_p_from_a.a_ctg_fasta_fn,
      merged_paf_fn = purge_dups_map_merge.merged_paf_fn,
      config_sh_fn = generate_config.config_sh,
      num_threads = nproc,  # for minimap2
      log_level = log_level,
      mem_scale_factor = mem_scale_factor,
      base_memory_mb = add_memory_mb,
  }

  call report_assembly2 {
    input:
      contigs_fasta = purge_dups.p_ctg_fasta,
      haplotigs_fasta = purge_dups.a_ctg_fasta,
      circular_contigs_txt = assemble.circular_contigs,
      log_level = log_level,
      mem_scale_factor = mem_scale_factor,
  }

  if (ipa2_cleanup_intermediate_files) {
    Array[String] files_to_remove = [
        #build_db.seqdb_seqs_fn,
        #build_db.seeddb_seeds_fn,
        ovl_asym_merge.m4_merged_raw_fn,
        ovl_asym_merge.m4_filtered_nonlocal_fn,
        phasing_merge.gathered_m4_fn,
        phasing_merge.all_keeps_m4_fn,
        phasing_merge.all_scraps_m4_fn,
        ovl_filter.m4_final_fn,
        ovl_filter.m4_chimerfilt_fn,
        polish_merge.consensus_merged_fn,
        #purge_dups.p_ctg_fasta_fn,
        #purge_dups.a_ctg_fasta_fn,
    ]

    call cleanup_files {
      input:
          singular_fns = files_to_remove,
          ovl_run_m4_fns = ovl_asym_run.out_m4_fn,
          ovl_run_sorted_m4_fns = ovl_asym_run.out_sorted_m4_fn,
          phasing_run_keep_fns = phasing_run.out_keep_m4_fn,
          phasing_run_scraps_fns = phasing_run.out_scraps_m4_fn,
          polish_run_consensus_fns = polish_run.consensus_fn,
          phasing_prepare_dn = shard_ids_phasing_dn,
          purge_dups_map_run_fns = purge_dups_map_run.paf_fn,
          purge_dups_fn = purge_dups.p_ctg_fasta,
          num_threads = 1,
    }
  }

  output {
    File final_primary_contigs_fasta = purge_dups.p_ctg_fasta
    File final_haplotigs_fasta = purge_dups.a_ctg_fasta
    # Maybe add .fai files someday.
    File circular_contigs = assemble.circular_contigs
    File report_polished_assembly = report_assembly2.report
    String full_seqdb_fn = build_db.seqdb_fn
  }
}
