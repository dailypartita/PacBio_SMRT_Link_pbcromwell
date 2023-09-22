version 1.0

import "pb_assembly_clr.wdl" as wf_clr
import "tasks/assembly.wdl" as ipa
import "wf_graph_mapping.wdl" as wf_raptor
import "wf_mapping.wdl" as wf_mapping
import "wf_consensus.wdl" as wf_consensus
import "tasks/pbcoretools.wdl" as pbcoretools
import "tasks/pbreports.wdl" as pbreports

task generate_config {
    input {
        String advanced_opt_str
        Int length_cutoff
        Int downsampled_coverage
        String genome_size
        Int coverage
        Int plasmid_contig_len_max
    }
    command {
        params_advanced_opt="${advanced_opt_str}" \
        params_length_cutoff="${length_cutoff}" \
        params_subsample_coverage="${downsampled_coverage}" \
        params_genome_size="${genome_size}" \
        params_coverage="${coverage}" \
        params_max_plasmid_len="${plasmid_contig_len_max}" \
            ipa-task generate_config ipa-task
    }
    runtime {
      cpu: 1
      memory: "100MB"
    }
    output {
        File config_stage1_sh_fn = "generated.config.stage1.bash"
        File config_stage1_json_fn = "generated.config.stage1.json"
        File config_stage2_sh_fn = "generated.config.stage2.bash"
        File config_stage2_json_fn = "generated.config.stage2.json"
    }
}

task filter_draft_contigs {
    input {
        File p_ctg_fa
        File config_sh_fn
    }
    command {
        input_fasta="${p_ctg_fa}" \
        output_fasta="p_ctg_filtered.fasta" \
        params_config_sh_fn=${config_sh_fn} \
            ipa-task filter_draft_contigs
    }
    runtime {
      cpu: 1
      memory: "200MB"
    }
    output {
        File p_ctg_filtered_fa = "p_ctg_filtered.fasta"
    }
}

task preassembly_stats {
    input {
        File rdb_rawreads_fn
        File rdb_ercreads_fn
        File length_cutoff_fn
        String genome_size
    }
    command {
        input_rawreads_rdb="${rdb_rawreads_fn}" \
        input_ercreads_rdb="${rdb_ercreads_fn}" \
        input_length_cutoff_fn=${length_cutoff_fn} \
        output_json="preassembly.report.json" \
        params_genome_size="${genome_size}" \
            ipa-task preassembly_stats
    }
    runtime {
      cpu: 1
      memory: "200MB"
    }
    output {
        File report_preassembly = "preassembly.report.json"
    }
}

task get_unaligned_reads_from_dataset {
    input {
        File dataset
        File config_sh_fn
    }
    command {
        input_dataset="${dataset}" \
        params_config_sh_fn=${config_sh_fn} \
            ipa-task get_unaligned_reads_from_dataset
    }
    runtime {
      cpu: 1
      memory: "1GB"
    }
    output {
        #gathered_fofn = "bams.fofn"
        File exclude_subreads = "exclude-subreads.txt"
    }
}

task collect_contigs {
    input {
        File chr_contigs
        File plasmid_contigs
    }
    command {
        output_collected_contigs_fa=collected_ctg.fasta \
        input_chr_contigs_fa="${chr_contigs}" \
        input_plasmid_contigs_fa="${plasmid_contigs}" \
            ipa-task collect_contigs
    }
    runtime {
      cpu: 1
      memory: "100MB"
    }
    output {
        File collected_ctgs = "collected_ctg.fasta"
    }
}

task collect_circ_list {
    input {
        File chr_circ_fn
        File plasmid_circ_fn
    }
    command {
        output_collected_circ=collected_circ.txt \
        input_chr_circ="${chr_circ_fn}" \
        input_plasmid_circ="${plasmid_circ_fn}" \
            ipa-task collect_circ_list
    }
    runtime {
      cpu: 1
      memory: "100MB"
    }
    output {
        File collected = "collected_circ.txt"
    }
}

task collect_gfa {
    input {
        File chr_gfa
        File plasmid_gfa
    }
    command {
        input_chr_gfa="${chr_gfa}" \
        input_plasmid_gfa="${plasmid_gfa}" \
            ipa-task collect_gfa
    }
    runtime {
      cpu: 1
      memory: "100MB"
    }
    output {
        File collected = "all_contig.gfa2"
    }
}

task dedup_plasmids {
    input {
        File chr_contigs
        File plasmid_contigs
        File config_sh_fn
    }
    command {
        output_deduped_contigs_fa=deduped_plasmid_ctgs.fasta \
        input_chr_contigs_fa="${chr_contigs}" \
        input_plasmid_contigs_fa="${plasmid_contigs}" \
        params_config_sh_fn=${config_sh_fn} \
            ipa-task dedup_plasmids
    }
    runtime {
      cpu: 1
      memory: "1GB"
    }
    output {
        File out_ctgs = "deduped_plasmid_ctgs.fasta"
    }
}

task rename_for_ncbi {
    input {
        File circ_list_fn
        File ctg_fa
    }
    command {
        input_circular_fn="${circ_list_fn}" \
        input_fa="${ctg_fa}" \
            ipa-task rename_for_ncbi
    }
    runtime {
      cpu: 1
      memory: "100MB"
    }
    output {
        File assembly = "assembly.rotated.polished.renamed.fsa"
    }
}

workflow pb_assembly_microbial {
    input {
        File eid_subread
        Int microasm_length_cutoff = -1
        Int microasm_downsampled_coverage = 100
        String microasm_genome_size = "5M"
        Int microasm_coverage = 30
        String microasm_advanced_options = ""

        Int microasm_plasmid_contig_len_max = 300000
        String? mapping_biosample_name

        Int nproc = 1
        String log_level = "INFO"
        # Use max_nchunks = 1 for a local run to avoid resource competition. On a cluster, the larger the better.
        Int max_nchunks = 40
        String tmp_dir = "/tmp"

    }
    # XXX for memory allocation in consensus workflow - we can just hardcode
    # an arbitrary max since the memory savings from a more exact number is
    # negligible
    Int genome_length_mb_max = 20

    call generate_config {
        input:
            advanced_opt_str = microasm_advanced_options,
            length_cutoff = microasm_length_cutoff,
            downsampled_coverage = microasm_downsampled_coverage,
            genome_size = microasm_genome_size,
            coverage = microasm_coverage,
            plasmid_contig_len_max = microasm_plasmid_contig_len_max,
    }

    ### Assemble long contigs.
    call wf_clr.pb_assembly_clr as asm_chrom {
        input:
            eid_subread=eid_subread,
            nproc=nproc,
            log_level=log_level,
            ctg_prefix="ctg.s1.",
            #blacklist_fofn=blacklist_fofn,

            blacklist_fofn=eid_subread, # phony input, as a sentinel

            config_sh_fn=generate_config.config_stage1_sh_fn,

            # wdl2json currently does not handle File, except eid_*
            # But we need an optional File for the first iteration of the sub-workflow.
            # So we will use eid_subread as a sentinel, for now.
    }

    call preassembly_stats {
        input:
            rdb_rawreads_fn = asm_chrom.raw_rdb,
            rdb_ercreads_fn = asm_chrom.erc_rdb,
            length_cutoff_fn = asm_chrom.raw_length_cutoff_fn,
            genome_size = microasm_genome_size,
    }

    call filter_draft_contigs {
        input:
            p_ctg_fa=asm_chrom.p_ctg_fa,
            config_sh_fn=generate_config.config_stage1_sh_fn,
    }

    call ipa.create_datasets {
        input:
            p_ctg_fa = filter_draft_contigs.p_ctg_filtered_fa,
            reads_fofn=eid_subread,
    }

    call wf_mapping.mapping as mapping_chr {
        input:
            reads = create_datasets.subreadset_xml,
            reference = create_datasets.contigset_xml,
            alignment_ext = ".alignmentset.xml",
            mapping_stats_module = "mapping_stats_hgap",
            coverage_report_module = "coverage_hgap",
            min_concordance = 80, # vs. 70 ???
            min_length = 50, # default=50?
            biosample_name = mapping_biosample_name,
            nproc = nproc,
            run_coverage = false,
            run_mapping_stats = false
    }

    call get_unaligned_reads_from_dataset {
        input:
            #dataset = rules.create_datasets.output.subreadset_xml,
            dataset = mapping_chr.mapped,
            config_sh_fn=generate_config.config_stage1_sh_fn,
    }

    ### Assemble plasmids.
    call wf_clr.pb_assembly_clr as asm_plasmid {
        input:
            eid_subread = eid_subread,
            nproc = nproc,
            log_level=log_level,
            blacklist_fofn=get_unaligned_reads_from_dataset.exclude_subreads,
            ctg_prefix="ctg.s2.",

            config_sh_fn=generate_config.config_stage2_sh_fn,
    }

    call dedup_plasmids {
        input:
            chr_contigs = filter_draft_contigs.p_ctg_filtered_fa,
            plasmid_contigs = asm_plasmid.p_ctg_fa,
            config_sh_fn=generate_config.config_stage1_sh_fn,
    }

    call collect_contigs {
        input:
            chr_contigs = filter_draft_contigs.p_ctg_filtered_fa,
            plasmid_contigs = dedup_plasmids.out_ctgs,
    }

    call collect_gfa {
        input:
            chr_gfa = asm_chrom.contig_gfa2,
            plasmid_gfa = asm_plasmid.contig_gfa2,
    }

    call collect_circ_list {
        input:
            chr_circ_fn = asm_chrom.circular_contigs,
            plasmid_circ_fn = asm_plasmid.circular_contigs,
    }

    call ipa.create_datasets as create_datasets_gathered {
        input:
            p_ctg_fa = collect_contigs.collected_ctgs,
            reads_fofn=eid_subread,
    }

    ### Polish.
    ### Alternative - Pbmm2.
#   TODO: We need to include both draft assemblies in the dataset for polishing.
#     call wf_mapping.mapping as mapping_all {
#         input:
#             reads = create_datasets_gathered.subreadset_xml,
#             reference = create_datasets_gathered.contigset_xml,
#             mapping_stats_module = "mapping_stats_hgap",
#             coverage_report_module = "coverage_hgap",
#             min_concordance = 0.80,
#             min_length = 500,
#             nproc = nproc,
#             max_nchunks = max_nchunks,
#             log_level = "INFO",
#     }

    # we need to use select_first here because the underlying implementation
    # has this as optional, but in this context it's always defined
    Int subreads_index_size_gb = select_first([mapping_chr.index_memory_gb])
    call wf_raptor.RaptorGraphMapping as mapping_all {
        input:
            reads_xml = create_datasets_gathered.subreadset_xml,
            reference = collect_contigs.collected_ctgs,
            referenceset_xml = create_datasets_gathered.contigset_xml,
            gfa = collect_gfa.collected,
            param_overrides = "--align",
            mapping_stats_module = "mapping_stats_hgap",
            coverage_report_module = "coverage_hgap",
            min_concordance = 80,
            min_length = 50,
            index_size_gb = subreads_index_size_gb,
            nproc = nproc,
            max_nchunks = 1
    }

    call pbcoretools.auto_consolidate_alignments {
      input:
        mapped = mapping_all.mapped,
        force_consolidate = true,
        tmp_dir = tmp_dir,
        log_level = log_level
    }

    call wf_consensus.consensus {
        input:
            alignments = mapping_all.mapped,
            reference = create_datasets_gathered.contigset_xml,
            consensus_algorithm = "arrow",
            genome_length_mb = genome_length_mb_max,
            log_level = "INFO",
            nproc = nproc,
            max_nchunks = max_nchunks,
    }

    call ipa.oric {
        input:
            whitelist=collect_circ_list.collected,
            p_ctg_fa=consensus.consensus_fasta,
            p_ctg_fq=consensus.consensus_fastq,
    }

    call ipa.pbreports_contig_table {
        input:
            contig_fa_fn=consensus.consensus_fasta,
            circ_list_fn=collect_circ_list.collected,
            gff_fn=mapping_all.coverage_gff,
    }

    call rename_for_ncbi {
        input:
            circ_list_fn = collect_circ_list.collected,
            ctg_fa = oric.p_ctg_oric_fa,
    }

    call pbreports.polished_assembly {
        input:
            coverage_gff = mapping_all.coverage_gff,
            consensus_fastq = consensus.consensus_fastq,
            report_contigs = pbreports_contig_table.report,
            nproc = nproc,
            log_level = "INFO"
    }

    call ipa.make_faidx {
      input:
        fasta = collect_contigs.collected_ctgs
    }

    output {
        File assembled_fasta = oric.p_ctg_oric_fa
        File assembled_fastq = oric.p_ctg_oric_fq
        File ncbi_fasta = rename_for_ncbi.assembly
        File circular_list = collect_circ_list.collected
        File mapped = mapping_all.mapped
        File coverage_gff = mapping_all.coverage_gff
        File report_mapping_stats = mapping_all.report_mapping_stats
        File report_coverage = mapping_all.report_coverage
        File consensus_fasta = consensus.consensus_fasta
        File consensus_fastq = consensus.consensus_fastq
        File report_polished_assembly = polished_assembly.report
        File? mapped_bam_datastore = auto_consolidate_alignments.datastore
        File draft_assembly = collect_contigs.collected_ctgs
        File draft_assembly_fai = make_faidx.fasta_fai
        File report_preassembly = preassembly_stats.report_preassembly
    }
}
