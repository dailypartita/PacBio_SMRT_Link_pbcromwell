version 1.0

import "pb_assembly_hifi.wdl" as wf_hifi
import "tasks/assembly.wdl" as ipa
import "wf_graph_mapping.wdl" as wf_raptor
import "wf_coverage_reports.wdl"
import "tasks/pbreports.wdl"
import "tasks/pbmm2.wdl"
import "tasks/pbcoretools.wdl" as pbcoretools
#import "tasks/pbreports.wdl" as pbreports
import "wf_prepare_input.wdl" as wf_prepare_input

task generate_config {
    input {
        String advanced_opt_str
        Int downsampled_coverage
        String genome_size
        Int plasmid_contig_len_max
    }
    command {
        params_advanced_opt="${advanced_opt_str}" \
        params_length_cutoff="-1" \
        params_subsample_coverage="${downsampled_coverage}" \
        params_genome_size="${genome_size}" \
        params_coverage="30" \
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

task get_unaligned_reads_from_dataset_and_seqdb {
    input {
        File dataset
        File config_sh_fn  # unused
        String full_seqdb_fn
    }
    command {
        set -vex
        input_dataset_fn="${dataset}" \
        input_full_seqdb_fn=${full_seqdb_fn} \
            ipa-task get_unaligned_reads_from_dataset_and_seqdb
        # Implicitly wrote into "./unaligned.fasta"
    }
    runtime {
      cpu: 1
      memory: "4GB"
    }
    output {
        #gathered_fofn = "bams.fofn"
        File unaligned_reads_fasta = "unaligned.fasta"
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
            ipa2-task collect_contigs
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
            ipa2-task collect_circ_list
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
            ipa2-task collect_gfa
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

task assembly_cleanup {
    input {
        File circ_list_fn
        File ctg_fa
    }
    command {
        set -vex
        ln -s "${ctg_fa}" final_assembly.fasta
        input_circular_fn="${circ_list_fn}" \
        input_fa="${ctg_fa}" \
            ipa-task rename_for_ncbi
        if [ -s "assembly.rotated.polished.renamed.fsa" ]; then
          samtools faidx assembly.rotated.polished.renamed.fsa
          samtools faidx final_assembly.fasta
        fi
        # FIXME pbmm2 won't recognize .fsa so we have to use ctg_fa instead
        dataset create \
          --type ReferenceSet \
          --name "Final Assembly" \
          assembly.referenceset.xml \
          final_assembly.fasta
    }
    runtime {
      cpu: 1
      memory: "100MB"
    }
    output {
        File assembly = "final_assembly.fasta"
        File? assembly_fai = "final_assembly.fasta.fai"
        File assembly_xml = "assembly.referenceset.xml"
        File ncbi_fasta = "assembly.rotated.polished.renamed.fsa"
    }
}

# copied from tasks/assembly.wdl
# added "?"
task pbreports_polished_assembly {
    input {
        File contig_fa_fn
        File circ_list_fn
        File? gff_fn # actually required
    }
    command {
        input_contig_fa_fn="${contig_fa_fn}" \
        input_circ_fn="${circ_list_fn}" \
        input_gff_fn="${gff_fn}" \
        output_json="polished_assembly.report.json" \
            ipa2-task pbreports_polished_assembly
    }
    runtime {
        cpu: 1
        memory: "200MB"
    }
    output {
        File report = "polished_assembly.report.json"
    }
}

task polish_contigs {
    input {
        File contig_fa_fn
        File circ_list_fn
        String full_seqdb_fn
        File mapped_bam_fn
        Boolean run_polishing
        Int nproc
        Int mem_scale_factor = 8
    }
    command <<<
        params_polish_run=~{true="1" false="0" run_polishing}
        echo "params_polish_run = ${params_polish_run}"

        if [[ ${params_polish_run} -eq 1 ]]; then
            # The input mapped file is in the XML format. This converts it to SAM.
            dataset summarize ~{mapped_bam_fn} | grep -E "*.bam$" | grep -v ":" > bam.fofn
            samtools merge -b bam.fofn --threads ~{nproc} mapped.merged.bam
            samtools view -h mapped.merged.bam > aln.mapped.merged.sam

            # Extract the FASTA read sequences.
            samtools fasta mapped.merged.bam > reads.fasta

            # Run consensus.
            racon -t ~{nproc} reads.fasta aln.mapped.merged.sam ~{contig_fa_fn} > polished_assembly.fasta

            # Cleanup.
            rm -f mapped.merged.bam
            rm -f aln.mapped.merged.sam
            rm -f reads.fasta
        else
            cp ~{contig_fa_fn} polished_assembly.fasta
        fi

        # Index, if not empty.
        touch polished_assembly.fasta.fai
        if [[ -s "polished_assembly.fasta" ]]; then
            samtools faidx polished_assembly.fasta
        fi

        # Rename the polished contigs for NCBI.
        touch assembly.rotated.polished.renamed.fsa
        touch assembly.rotated.polished.renamed.fsa.fai
        input_circular_fn="~{circ_list_fn}" \
        input_fa="~{contig_fa_fn}" \
            ipa-task rename_for_ncbi
        if [ -s "assembly.rotated.polished.renamed.fsa" ]; then
          samtools faidx assembly.rotated.polished.renamed.fsa
        fi
    >>>
    # Note: memory consumption of this task directly correlates with the size of the input data.
    Int mem_total_mb = 2048 * (mem_scale_factor + 1)
    runtime {
        cpu: nproc
        memory: "${mem_total_mb}MB"
    }
    output {
        File out_contigs_fa = "polished_assembly.fasta"
        File out_contigs_fa_fai = "polished_assembly.fasta.fai"
        File out_contigs_ncbi_fa = "assembly.rotated.polished.renamed.fsa"
        File out_contigs_ncbi_fa_fai = "assembly.rotated.polished.renamed.fsa.fai"
    }
}

workflow pb_assembly_hifi_microbial {
    input {
        File? eid_ccs
        File? reads

        String ipa2_genome_size = "10M"
        Int ipa2_downsampled_coverage = 100
        String ipa2_advanced_options_chrom = "config_block_size = 100; config_seeddb_opt = -k 28 -w 20 --space 0 --use-hpc-seeds-only; config_ovl_opt = --smart-hit-per-target --min-idt 98 --traceback --mask-hp --mask-repeats --trim --trim-window-size 30 --trim-match-frac 0.75;"
        String ipa2_advanced_options_plasmid = "config_block_size = 100; config_ovl_filter_opt = --max-diff 80 --max-cov 100 --min-cov 2 --bestn 10 --min-len 500 --gapFilt --minDepth 4 --idt-stage2 98; config_ovl_min_len = 500; config_seeddb_opt = -k 28 -w 20 --space 0 --use-hpc-seeds-only; config_ovl_opt = --smart-hit-per-target --min-idt 98 --min-map-len 500 --min-anchor-span 500 --traceback --mask-hp --mask-repeats --trim --trim-window-size 30 --trim-match-frac 0.75 --smart-hit-per-target --secondary-min-ovl-frac 0.05; config_layout_opt = --allow-circular;"
        Boolean ipa2_cleanup_intermediate_files = true
        Int microasm_plasmid_contig_len_max = 300000 # 300kB
        Boolean microasm_run_secondary_polish = true

        String? mapping_biosample_name

        String dataset_filters = ""
        Int filter_min_qv = 20
        Int downsample_factor = 0
        Int add_memory_mb = 0

        Int nproc = 1
        String log_level = "INFO"
        # Use max_nchunks = 1 for a local run to avoid resource competition. On a cluster, the larger the better.
        Int max_nchunks = 40
        String tmp_dir = "/tmp"

    }
    # XXX for memory allocation in the pbmm2 task - we can just hardcode
    # an arbitrary max since the memory savings from a more exact number is
    # negligible
    Int genome_length_mb_max = 20
    # default is 70
    Float pbmm2_min_concordance = 98
    Int pbmm2_min_length = 500

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
            advanced_opt_str = "",
            downsampled_coverage = ipa2_downsampled_coverage,
            genome_size = ipa2_genome_size,
            plasmid_contig_len_max = microasm_plasmid_contig_len_max,
    }

    ### Assemble long contigs.
    call wf_hifi.pb_assembly_hifi as asm_chrom {
        input:
            eid_ccs = prepare_input.reads_file,
            nproc = nproc,
            log_level = log_level,
            ipa2_ctg_prefix = "ctg.s1",
            ipa2_run_phasing = true,
            ipa2_run_polishing = true,
            ipa2_run_purge_dups = false,
            ipa2_genome_size = ipa2_genome_size,
            ipa2_downsampled_coverage = ipa2_downsampled_coverage,
            ipa2_advanced_options = ipa2_advanced_options_chrom,
            ipa2_cleanup_intermediate_files = ipa2_cleanup_intermediate_files,
            mem_scale_factor = 1,
            add_memory_mb = add_memory_mb
            #config_sh_fn = generate_config.config_stage1_sh_fn,
    }

    call filter_draft_contigs {
        input:
            p_ctg_fa=asm_chrom.final_primary_contigs_fasta,
            config_sh_fn=generate_config.config_stage1_sh_fn,
    }

    call ipa.create_datasets {
        input:
            p_ctg_fa = filter_draft_contigs.p_ctg_filtered_fa,
            reads_fofn=prepare_input.reads_file,
    }

    call pbmm2.pbmm2_align as mapping_chr {
        input:
            unmapped = prepare_input.reads_file,
            reference = create_datasets.contigset_xml,
            aln_ext = ".consensusalignmentset.xml",
            min_concordance = pbmm2_min_concordance,
            # XXX why is this different from mapping_all?
            min_length = 50, # default=50?
            biosample_name = mapping_biosample_name,
            genome_length_mb = genome_length_mb_max,
            base_memory_mb = add_memory_mb,
            preset_mode = "HIFI",
            nproc = nproc,
            log_level = log_level
    }

    call get_unaligned_reads_from_dataset_and_seqdb {
        input:
            dataset = mapping_chr.mapped,
            full_seqdb_fn = asm_chrom.full_seqdb_fn,
            config_sh_fn=generate_config.config_stage1_sh_fn,
    }

    ### Assemble plasmids.
    call wf_hifi.pb_assembly_hifi as asm_plasmid {
        input:
            #eid_ccs = prepare_input.reads_file,
            reads = get_unaligned_reads_from_dataset_and_seqdb.unaligned_reads_fasta,
            nproc = nproc,
            log_level = log_level,
            ipa2_ctg_prefix = "ctg.s2",
            ipa2_run_phasing = true,
            ipa2_run_polishing = true,
            ipa2_run_purge_dups = false,
            ipa2_genome_size = ipa2_genome_size,
            ipa2_downsampled_coverage = ipa2_downsampled_coverage,
            ipa2_advanced_options = ipa2_advanced_options_plasmid,
            ipa2_cleanup_intermediate_files = ipa2_cleanup_intermediate_files,
            mem_scale_factor = 1,
            add_memory_mb = add_memory_mb,

            #config_sh_fn = generate_config.config_stage2_sh_fn,
    }

    call dedup_plasmids {
        input:
            chr_contigs = filter_draft_contigs.p_ctg_filtered_fa,
            plasmid_contigs = asm_plasmid.final_primary_contigs_fasta,
            config_sh_fn=generate_config.config_stage1_sh_fn,
    }

    call collect_contigs {
        input:
            chr_contigs = filter_draft_contigs.p_ctg_filtered_fa,
            plasmid_contigs = dedup_plasmids.out_ctgs,
    }

    #call collect_gfa {
    #    input:
    #        chr_gfa = asm_chrom.contig_gfa2,
    #        plasmid_gfa = asm_plasmid.contig_gfa2,
    #}

    call collect_circ_list {
        input:
            chr_circ_fn = asm_chrom.circular_contigs,
            plasmid_circ_fn = asm_plasmid.circular_contigs,
    }

    call ipa.oric {
        input:
            whitelist=collect_circ_list.collected,
            p_ctg_fa=collect_contigs.collected_ctgs,
    }

    call assembly_cleanup {
        input:
            circ_list_fn = collect_circ_list.collected,
            ctg_fa = oric.p_ctg_oric_fa,
    }

    call pbreports_polished_assembly {
        input:
            contig_fa_fn = assembly_cleanup.assembly,
            circ_list_fn = collect_circ_list.collected,
            gff_fn = summarize_coverage.coverage_gff
    }

    # Unlike the old CLR workflow, this mapping is not used for polishing, it
    # just provides the user a BAM file to look at in IGV, plus the mapping
    # and coverage reports.  We can also use the mapping results for
    # subsequent basemods analysis (FUTURE).
    call pbmm2.pbmm2_align as mapping_all {
        input:
            unmapped = prepare_input.reads_file,
            reference = assembly_cleanup.assembly_xml,
            aln_ext = ".consensusalignmentset.xml",
            min_concordance = pbmm2_min_concordance,
            min_length = pbmm2_min_length,
            genome_length_mb = genome_length_mb_max,
            base_memory_mb = add_memory_mb,
            preset_mode = "HIFI",
            nproc = nproc,
            log_level = log_level
    }

    call polish_contigs {
        input:
            contig_fa_fn = assembly_cleanup.assembly,
            circ_list_fn = collect_circ_list.collected,
            full_seqdb_fn = asm_chrom.full_seqdb_fn,
            mapped_bam_fn = mapping_all.mapped,
            run_polishing = microasm_run_secondary_polish,
            nproc = nproc
    }

    call pbreports.mapping_stats {
        input:
            mapped = mapping_all.mapped,
            all_reads = prepare_input.reads_file,
            report_module = "mapping_stats",
            show_calibration_plot = false,
            min_concordance = pbmm2_min_concordance,
            index_memory_gb = 4,
            base_memory_mb = add_memory_mb,
            nproc = nproc,
            log_level = log_level
    }

    call wf_coverage_reports.gc_coverage_plot {
        input:
            alignments = mapping_all.mapped,
            reference = assembly_cleanup.assembly_xml,
            base_memory_mb = add_memory_mb,
            nproc = nproc,
            log_level = log_level
    }

    call wf_coverage_reports.summarize_coverage {
        input:
            mapped = mapping_all.mapped,
            reference = assembly_cleanup.assembly_xml,
            base_memory_mb = add_memory_mb,
            nproc = nproc,
            log_level = log_level
    }

    call wf_coverage_reports.pbreports_coverage {
        input:
            reference = assembly_cleanup.assembly_xml,
            coverage_gff = summarize_coverage.coverage_gff,
            gc_coverage_plot_png = gc_coverage_plot.plot_png,
            base_memory_mb = add_memory_mb,
            nproc = nproc,
            log_level = log_level
    }

    output {
        File assembly_fasta = polish_contigs.out_contigs_fa
        File? assembly_fasta_fai = polish_contigs.out_contigs_fa_fai
        File ncbi_fasta = polish_contigs.out_contigs_ncbi_fa
        File circular_list = collect_circ_list.collected
        File mapped = mapping_all.mapped
        File mapped_target_fasta = assembly_cleanup.assembly
        File? mapped_target_fasta_fai = assembly_cleanup.assembly_fai
        File report_polished_assembly = pbreports_polished_assembly.report
        File coverage_gff = summarize_coverage.coverage_gff
        File report_mapping_stats = mapping_stats.report
        File report_coverage = pbreports_coverage.report
        File? mapped_bam_datastore = mapping_all.bam_datastore
    }
}
