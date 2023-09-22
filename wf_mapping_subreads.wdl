# XXX deprecated chunked version of mapping workflow

version 1.0

import "tasks/chunking.wdl" as chunking
import "tasks/pbreports.wdl" as pbreports
import "tasks/pbcoretools.wdl" as pbcoretools
import "tasks/pbmm2.wdl" as mapping_tasks
import "tasks/memory.wdl" as memory
import "wf_coverage_reports.wdl"

workflow mapping_subreads {
  input {
    File reads
    File reference
    File? target_regions_bed
    String dataset_split_mode = "zmws"
    String preset_mode = "SUBREAD"
    Float min_concordance = 70
    Int min_length = 50
    Boolean hq_mode = false
    Boolean zmw_mode = false
    Boolean? median_filter = false
    Boolean? strip = false
    Boolean? split_by_sample = false
    String? biosample_name
    String? pbmm2_overrides
    Boolean run_coverage = true
    Boolean run_mapping_stats = true
    Boolean run_gc_coverage_plot = true
    Boolean run_bam_auto_consolidation = false # this just enables auto behavior
    Boolean consolidate_aligned_bam = false # this is the user-facing option
    # RQ vs. concordance, CCS mapping only
    Boolean? show_calibration_plot

    Int nproc
    Int max_nchunks
    Int target_size
    String log_level = "INFO"
    String? tmp_dir
    Int mem_scale_factor = 8
    Int add_memory_mb = 0
  }
  # this is set higher for internal purposes only
  Int base_memory_mb = add_memory_mb + 16384
  String alignment_ext = ".alignmentset.xml"
  String mapping_stats_module = "mapping_stats_subreads"
  String coverage_report_module = "coverage"

  call memory.get_input_sizes {
    input:
      bam_dataset = reads,
      fasta_dataset = reference,
      log_level = log_level
  }

  if (max_nchunks > 1) {
    call chunking.split_dataset as split_reads {
      input:
        ds_in = reads,
        max_nchunks = max_nchunks,
        target_size = target_size,
        split_mode = dataset_split_mode,
        memory_gb = get_input_sizes.index_memory_gb,
        log_level = log_level
    }

    scatter (chunk in split_reads.chunks) {
      call mapping_tasks.pbmm2_align as pbmm2_chunk {
        input:
          unmapped = chunk,
          reference = reference,
          aln_ext = alignment_ext,
          preset_mode = preset_mode,
          min_concordance = min_concordance,
          min_length = min_length,
          hq_mode = hq_mode,
          zmw_mode = zmw_mode,
          median_filter = median_filter,
          strip = strip,
          split_by_sample = split_by_sample,
          biosample_name = biosample_name,
          pbmm2_overrides = pbmm2_overrides,
          genome_length_mb = get_input_sizes.genome_length_mb,
          mem_scale_factor = mem_scale_factor,
          base_memory_mb = base_memory_mb,
          min_memory_mb = 8192,
          nproc = nproc,
          log_level = log_level,
          tmp_dir = tmp_dir
      }
    }

    call chunking.gather_datasets_lite as gather_alignments {
      input:
        chunks = pbmm2_chunk.mapped,
        dataset_name = "'Mapped Reads'",
        dataset_ext = alignment_ext,
        prefix = "mapped",
        log_level = log_level
    }
  }
  # unchunked version
  if (max_nchunks <= 1) {
    call mapping_tasks.pbmm2_align as pbmm2_all {
      input:
        unmapped = reads,
        reference = reference,
        aln_ext = alignment_ext,
        preset_mode = preset_mode,
        min_concordance = min_concordance,
        min_length = min_length,
        hq_mode = hq_mode,
        zmw_mode = zmw_mode,
        median_filter = median_filter,
        strip = strip,
        split_by_sample = split_by_sample,
        biosample_name = biosample_name,
        pbmm2_overrides = pbmm2_overrides,
        genome_length_mb = get_input_sizes.genome_length_mb,
        mem_scale_factor = mem_scale_factor,
        base_memory_mb = base_memory_mb,
        min_memory_mb = 8192,
        nproc = nproc,
        log_level = log_level,
        tmp_dir = tmp_dir
    }
  }
  File mapped_out = select_first([gather_alignments.gathered,
                                  pbmm2_all.mapped])

  # for workflows like pb_assembly_microbial that run two separate mapping
  # steps, we only want to generate one report
  if (run_mapping_stats) {
    call pbreports.mapping_stats {
      input:
        mapped = mapped_out,
        all_reads = reads,
        report_module = mapping_stats_module,
        show_calibration_plot = show_calibration_plot,
        min_concordance = min_concordance,
        index_memory_gb = get_input_sizes.index_memory_gb,
        base_memory_mb = base_memory_mb,
        nproc = nproc,
        log_level = log_level
    }
  }

  # warning: this is slow, so some workflows run these tasks separately
  if (run_coverage) {
    call wf_coverage_reports.coverage_reports {
      input:
        mapped = mapped_out,
        reference = reference,
        target_regions_bed = target_regions_bed,
        index_memory_gb = get_input_sizes.index_memory_gb,
        base_memory_mb = base_memory_mb,
        run_gc_coverage_plot = run_gc_coverage_plot,
        nproc = nproc,
        log_level = log_level
    }
  }

  # this is even slower, so it is also usually called elsewhere
  if (run_bam_auto_consolidation) {
    call pbcoretools.auto_consolidate_alignments {
      input:
        mapped = mapped_out,
        force_consolidate = consolidate_aligned_bam,
        dataset_ext = alignment_ext,
        log_level = log_level,
        tmp_dir = tmp_dir
    }
  }

  File mapped_final = select_first([auto_consolidate_alignments.dataset,
                                    mapped_out])

  output {
    File mapped = mapped_final
    File? report_mapping_stats = mapping_stats.report
    File? coverage_gff = coverage_reports.coverage_gff
    File? report_coverage = coverage_reports.report_coverage
    File? report_target_coverage = coverage_reports.report_target_coverage
    File? mapped_bam_datastore = auto_consolidate_alignments.datastore
    # for repeat analysis this is used in a different report
    File? plot_target_coverage_png = coverage_reports.plot_target_coverage_png
    Int? genome_length_mb = get_input_sizes.genome_length_mb
    Int? index_memory_gb = get_input_sizes.index_memory_gb
  }
}
