version 1.0

import "tasks/pbcoretools.wdl" as pbcoretools
import "wf_pbsv_intl.wdl"
import "wf_prepare_input.wdl"
import "wf_mapping.wdl"

workflow pb_sv_ccs {
  input {
    File eid_ccs
    File eid_ref_dataset

    # Filter option
    Int filter_min_qv = 20
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
    # pbsv v2.8.0, equivalent to --hifi
    Int min_percent_reads = 10                # -P 10
    Int min_reads_one_sample = 3              # -O 3
    Int min_reads_all_samples = 3             # -A 3
    Int min_reads_per_strand_all_samples = 0  # -S 0

    Int add_memory_mb = 0

    # Workflow resource configuration parameters
    Int nproc = 1
    String log_level = "INFO"
    Int max_nchunks = 100
    String? tmp_dir
  }

  # Mapping parameters
  Boolean mapping_median_filter = false
  Boolean mapping_strip = true
  Boolean mapping_split_by_sample = true

  call wf_prepare_input.prepare_input {
    input:
      dataset_xml = eid_ccs,
      dataset_ext = ".consensusreadset.xml",
      dataset_filters = dataset_filters,
      downsample_factor = downsample_factor,
      filter_min_qv = filter_min_qv,
      nproc = 1,
      log_level = log_level
  }

  call wf_mapping.mapping {
    input:
      reads = prepare_input.reads_file,
      reference = eid_ref_dataset,
      min_concordance = mapping_min_concordance,
      median_filter = mapping_median_filter,
      strip = mapping_strip,
      split_by_sample = mapping_split_by_sample,
      min_length = mapping_min_length,
      biosample_name = mapping_biosample_name,
      pbmm2_overrides = mapping_pbmm2_overrides,
      mem_scale_factor = 6,
      base_memory_mb = add_memory_mb,
      nproc = nproc,
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
      hifi_mode = true,
      base_memory_mb = add_memory_mb,
      nproc = nproc,
      log_level = log_level,
      max_nchunks = max_nchunks,
  }

  output {
    File variants = pbsv_intl.variants
    File alignments_by_sample_datastore = pbsv_intl.alignments_by_sample_datastore
    File report_structural_variants_2 = pbsv_intl.report
  }
}
