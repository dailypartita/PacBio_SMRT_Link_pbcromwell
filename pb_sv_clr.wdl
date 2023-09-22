version 1.0

import "tasks/pbcoretools.wdl" as pbcoretools
import "wf_pbsv_intl.wdl"
import "wf_mapping_subreads.wdl"

workflow pb_sv_clr {
  input {
    File eid_subread
    File eid_ref_dataset

    # Filter option
    String dataset_filters = ""
    Int downsample_factor = 0

    Int mapping_min_length = 50
    Float mapping_min_concordance = 70
    String? mapping_biosample_name
    String? mapping_pbmm2_overrides

    # pbsv call parameters
    String chunk_length = "1M"
    Int min_sv_length = 20
    String? pbsv_override_args
    # pbsv v2.2.2, CLR optimized parameters -A 2 -O 2 -S 1 -P 20
    Int min_percent_reads = 20                # -P 20
    Int min_reads_one_sample = 2              # -O 2
    Int min_reads_all_samples = 2             # -A 2
    Int min_reads_per_strand_all_samples = 1  # -S 1

    # Workflow resource configuration parameters
    Int nproc = 1
    String log_level = "INFO"
    Int max_nchunks = 100
    String? tmp_dir
  }

  # Mapping parameters
  Boolean mapping_median_filter = true
  Boolean mapping_strip = true
  Boolean mapping_split_by_sample = true

  call pbcoretools.dataset_filter {
    input:
      dataset = eid_subread,
      filters = dataset_filters,
      downsample_factor = downsample_factor,
      nproc = 1,
      log_level = log_level
  }

  call wf_mapping_subreads.mapping_subreads as mapping {
    input:
      reads = dataset_filter.filtered,
      reference = eid_ref_dataset,
      min_concordance = mapping_min_concordance,
      median_filter = mapping_median_filter,
      strip = mapping_strip,
      split_by_sample = mapping_split_by_sample,
      min_length = mapping_min_length,
      biosample_name = mapping_biosample_name,
      pbmm2_overrides = mapping_pbmm2_overrides,
      mem_scale_factor = 4,
      nproc = nproc,
      max_nchunks = max_nchunks,
      target_size = 100000,
      log_level = log_level
  }

  call wf_pbsv_intl.pbsv_intl {
    input:
      mapped = mapping.mapped,
      reference = eid_ref_dataset,
      chunk_length = chunk_length,
      min_sv_length = min_sv_length,
      min_percent_reads = min_percent_reads,
      min_reads_one_sample = min_reads_one_sample,
      min_reads_all_samples = min_reads_all_samples,
      min_reads_per_strand_all_samples = min_reads_per_strand_all_samples,
      pbsv_override_args = pbsv_override_args,
      hifi_mode = false,
      nproc = nproc,
      log_level = log_level,
      max_nchunks = max_nchunks
  }

  output {
    File variants = pbsv_intl.variants
    File alignments_by_sample_datastore = pbsv_intl.alignments_by_sample_datastore
    File report_structural_variants_2 = pbsv_intl.report
  }
}
