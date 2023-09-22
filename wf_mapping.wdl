# simple mapping pipeline, also used in resequencing

version 1.0

import "tasks/pbreports.wdl" as pbreports
import "tasks/pbmm2.wdl" as mapping_tasks
import "tasks/memory.wdl" as memory
import "wf_coverage_reports.wdl"

workflow mapping {
  Int TARGET_SIZE_CCS = 300000
  input {
    File reads
    File reference
    File? target_regions_bed
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
    # RQ vs. concordance, CCS mapping only
    Boolean? show_calibration_plot
    # the default is CCS now
    String preset_mode = "HiFi"
    # FIXME this is still used by both the CCS and CLR implementations of SV
    String alignment_ext = ".consensusalignmentset.xml"
    String mapping_stats_module = "mapping_stats"
    # this can also be coverage_hgap
    String coverage_report_module = "coverage"

    Int nproc
    String log_level = "INFO"
    String? tmp_dir
    Int mem_scale_factor = 6
    Int base_memory_mb = 0
  }

  call memory.get_input_sizes {
    input:
      bam_dataset = reads,
      fasta_dataset = reference,
      base_memory_mb = base_memory_mb,
      log_level = log_level
  }

  call mapping_tasks.pbmm2_align as pbmm2 {
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
      nproc = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir
  }

  # for workflows like pb_assembly_microbial that run two separate mapping
  # steps, we only want to generate one report
  if (run_mapping_stats) {
    call pbreports.mapping_stats {
      input:
        mapped = pbmm2.mapped,
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
        mapped = pbmm2.mapped,
        reference = reference,
        target_regions_bed = target_regions_bed,
        index_memory_gb = get_input_sizes.index_memory_gb,
        base_memory_mb = base_memory_mb,
        run_gc_coverage_plot = run_gc_coverage_plot,
        nproc = nproc,
        log_level = log_level
    }
  }

  if (defined(pbmm2.mapped_bam_fofn)) {
    String mapped_bam_fn = read_string(select_first([pbmm2.mapped_bam_fofn]))
    String mapped_bam_bai_fn = read_string(select_first([pbmm2.mapped_bam_bai_fofn]))
  }

  output {
    File mapped = pbmm2.mapped
    File? report_mapping_stats = mapping_stats.report
    File? coverage_gff = coverage_reports.coverage_gff
    File? report_coverage = coverage_reports.report_coverage
    File? report_target_coverage = coverage_reports.report_target_coverage
    File? mapped_bam_datastore = pbmm2.bam_datastore
    # for repeat analysis this is used in a different report
    File? plot_target_coverage_png = coverage_reports.plot_target_coverage_png
    Int? genome_length_mb = get_input_sizes.genome_length_mb
    Int? index_memory_gb = get_input_sizes.index_memory_gb
    String? mapped_bam = mapped_bam_fn
    String? mapped_bam_bai = mapped_bam_bai_fn
  }
}
