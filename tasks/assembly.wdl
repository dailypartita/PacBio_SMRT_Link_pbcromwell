version 1.0

task build_raw_rdb {
    input {
        File reads_fn
        String db_prefix
        File? blacklist_fofn
        File config_sh_fn
    }
    command {
        set -e
        input_reads_fn="`readlink -f ${reads_fn}`" \
        input_blacklist_fofn="${blacklist_fofn}" \
        params_config_sh_fn=${config_sh_fn} \
            ipa-task build_raw_rdb

        ln -sf reads.rdb ${db_prefix}.rdb
        ln -sf reads.full.rdb ${db_prefix}.full.rdb
    }
    runtime {
        cpu: 1
        memory: "500MB"
    }
    output {
        File out_rdb = "${db_prefix}.rdb"
        File out_rdb_full = "${db_prefix}.full.rdb"
        File out_cutoff = "length_cutoff"
    }
}

task ovlcns_raw_prepare {
    input {
        String rdb
    }
    command {
        output_blocks=./blocks \
        input_rdb="${rdb}" \
            ipa-task ovlcns_raw_prepare
    }
    runtime {
        cpu: 1
        memory: "500MB"
    }
    output {
        # defer reading this until main workflow resumes (AWS workaround)
        File blocks_fofn = "blocks"
    }
}

task ovlcns_raw_run {
    input {
        File rdb
        File length_cutoff_fn
        String db_prefix
        String block_id
        Int nproc
        File config_sh_fn

        Int base_mem_gb = 8
    }
    # This is an approximate upper limit
    Int mem_gb_per_core = 2

    command {
        input_rdb="${rdb}" \
        input_length_cutoff_fn="$rel/${length_cutoff_fn}" \
        params_block_id=${block_id} \
        params_num_threads=${nproc} \
        params_config_sh_fn=${config_sh_fn} \
            ipa-task ovlcns_raw_run

        ln -f ovl.m4 ovl.${db_prefix}.${block_id}.m4
        ln -f ovl.flipped.sorted.filtered.m4 ovl.${db_prefix}.flipped.sorted.filtered.${block_id}.m4
        ln -f ercreads.fasta ercreads.${block_id}.fasta
    }
    Int total_mem_gb = base_mem_gb + (nproc * mem_gb_per_core)
    runtime {
        cpu: nproc
        memory: "${total_mem_gb}GB"
    }

    output {
        File out_m4 = "ovl.rawreads.${block_id}.m4"
        File out_filtered_m4 = "ovl.${db_prefix}.flipped.sorted.filtered.${block_id}.m4"
        File out_ercreads_fasta = "ercreads.${block_id}.fasta"
    }
}

task ovlcns_raw_merge {
    input {
        Array[File] in_fns
    }
    command {
        set -vex
        echo ${sep=' ' in_fns} | xargs -n 1 > ercreads.fofn
        cat ${sep=' ' in_fns} > ercreads.fasta
    }
    runtime {
        cpu: 1
        memory: "1GB"
    }
    output {
        File gathered_fofn = "ercreads.fofn"
        File gathered_fasta = "ercreads.fasta"
    }
}

task build_erc_rdb {
    input {
        File reads_fn
        String db_prefix
        File config_sh_fn
    }

    command {
        input_reads_fn="`readlink -f ${reads_fn}`" \
        params_db_prefix="${db_prefix}" \
        params_config_sh_fn=${config_sh_fn} \
            ipa-task build_erc_rdb
    }
    runtime {
        cpu: 1
        memory: "500MB"
    }
    output {
        File out_rdb = "${db_prefix}.rdb"
    }
}

task ovl_erc_prepare {
    input {
        String rdb
    }
    command {
        output_blocks=./blocks \
        input_rdb="${rdb}" \
            ipa-task ovl_erc_prepare
    }
    runtime {
        cpu: 1
        memory: "200MB"
    }
    output {
        File blocks_file = "blocks"
    }
}

task ovl_erc_run {
    input {
        String rdb
        String block_id
        Int nproc
        File config_sh_fn
    }

    command {
        input_rdb="${rdb}" \
        params_block_id=${block_id} \
        params_num_threads=${nproc} \
        params_config_sh_fn=${config_sh_fn} \
            ipa-task ovl_erc_run

        ln -f ovl.m4 ovl.${block_id}.m4
        ln -f ovl.flipped.sorted.dovetail.m4 ovl.flipped.sorted.dovetail.${block_id}.m4
    }
    Int total_mem_gb = 4 + nproc
    runtime {
        cpu: nproc
        memory: "${total_mem_gb}GB"
    }
    output {
        File out_m4 = "ovl.${block_id}.m4"
        File out_filtered_m4 = "ovl.flipped.sorted.dovetail.${block_id}.m4"
    }
}

task ovl_erc_merge {
    input {
        Array[File] in_fns
        String db_prefix
        Int nproc
        File config_sh_fn
    }
    command {
        echo "${sep='\n' in_fns}" > ./merged.fofn

        input_fofn=./merged.fofn \
        output_ovl="~{db_prefix}.ovl" \
        params_num_threads="~{nproc}" \
        params_config_sh_fn=${config_sh_fn} \
            ipa-task ovl_erc_merge
    }
    runtime {
        cpu: nproc
        memory: "1GB"
    }
    output {
        File gathered_m4 = "${db_prefix}.ovl"
    }
}

task assemble {
    input {
        File in_m4
        File in_ercreads
        String ctg_prefix
        File config_sh_fn
    }

    command {
        input_m4="${in_m4}" \
        input_ercreads="${in_ercreads}" \
        params_ctg_prefix=${ctg_prefix} \
        params_config_sh_fn=${config_sh_fn} \
            ipa-task assemble
    }
    runtime {
        cpu: 1
        memory: "1GB"
    }
    output {
        File p_ctg_fa = "p_ctg.fasta"
        File a_ctg_fa = "a_ctg.fasta"
        File p_ctg_tiling_path = "p_ctg_tiling_path"
        File a_ctg_tiling_path = "a_ctg_tiling_path"
        File asm_gfa = "asm.gfa"
        File sg_gfa = "sg.gfa"
        File contig_gfa2 = "contig.gfa2"
        File circular_contigs = "circular_contigs.csv"
    }
}

task create_datasets {
    input {
        File p_ctg_fa
        File reads_fofn
    }
    command {
        input_fasta="${p_ctg_fa}" \
        input_reads="`readlink -f ${reads_fofn}`" \
            ipa-task create_datasets
    }
    runtime {
        cpu: 1
        memory: "1GB"
    }
    output {
        File subreadset_xml = "raw.subreadset.xml"
        File contigset_xml = "contigs.referenceset.xml"
    }
}

task oric {
  input{
    File whitelist
    File p_ctg_fa
    File? p_ctg_fq
  }
  command <<<
    input_fa="~{p_ctg_fa}" \
    input_whitelist="~{whitelist}" \
    output_fa="p_ctg_oric.fasta" \
    output_json="contigs_report.json" \
        ipa-task circ_oric

    # Run again for FASTQ, if found.
    # Note: We will skip contigs-report for FASTQ.
    if [[ -s "~{p_ctg_fq}" ]]; then
        input_fa="~{p_ctg_fq}" \
        input_whitelist="~{whitelist}" \
        output_fa="p_ctg_oric.fastq" \
            ipa-task circ_oric
    else
        # Just a harmless placeholder.
        touch "p_ctg_oric.fastq"
    fi
  >>>
  runtime {
    cpu: 1
    memory: "1GB"
  }
  output{
    File p_ctg_oric_fa = "p_ctg_oric.fasta"
    File p_ctg_oric_fq = "p_ctg_oric.fastq"
  }
}

task pbreports_contig_table {
    input {
        File contig_fa_fn
        File circ_list_fn
        File gff_fn
    }
    command {
        input_contig_fa_fn="${contig_fa_fn}" \
        input_circ_fn="${circ_list_fn}" \
        input_gff_fn="${gff_fn}" \
        output_json="contigs_report.json" \
            ipa-task pbreports_contig_table
    }
    runtime {
        cpu: 1
        memory: "200MB"
    }
    output {
        File report = "contigs_report.json"
    }
}

task make_faidx {
  input {
    File fasta

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    ln -s ${fasta}
    if [[ -s ${fasta} ]]; then
        samtools faidx `basename ${fasta}`
    else
        touch `basename ${fasta}`.fai
    fi
  }
  runtime {
    cpu: 1
    backend: "Local"
    memory: "100MB"
  }
  output {
    File fasta_fai = glob("*.fai")[0]
  }
}
