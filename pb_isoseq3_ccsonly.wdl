# Legacy version of isoseq3 HiFi workflow

version 1.0

import "pb_isoseq3.wdl" as pb_isoseq3

workflow pb_isoseq3_ccsonly {
  input {
    File eid_ccs
    File eid_barcode
    File? eid_ref_dataset

    File? biosamples_csv
    String? biosamples_csv_str
    String dataset_filters = ""
    Int filter_min_qv = 20
    Int add_memory_mb = 0

    # main window options
    Boolean run_clustering = true
    Boolean cluster_separately = false

    # Options exposed in smrtlink:
    Boolean refine_require_polya = true

    Int mapping_min_length = 50
    Float mapping_min_concordance = 95
    Float mapping_min_coverage = 99
    String? mapping_pbmm2_overrides

    Int isocollapse_max_fuzzy_junction = 5

    Int nproc = 1
    String log_level = "INFO"
    # this is completely ignored now
    Int max_nchunks = 0
    String? tmp_dir
  }

  call pb_isoseq3.pb_isoseq3 {
    input: 
      eid_ccs = eid_ccs,
      eid_barcode = eid_barcode,
      eid_ref_dataset = eid_ref_dataset,
      run_clustering = run_clustering,
      cluster_separately = cluster_separately,
      refine_require_polya = refine_require_polya,
      mapping_min_length = mapping_min_length,
      mapping_min_concordance = mapping_min_concordance,
      mapping_min_coverage = mapping_min_coverage,
      isocollapse_max_fuzzy_junction = isocollapse_max_fuzzy_junction,
      filter_min_qv = filter_min_qv,
      mapping_pbmm2_overrides = mapping_pbmm2_overrides,
      add_memory_mb = add_memory_mb,
      nproc = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir
  }

  output {
    File report_isoseq_primers = pb_isoseq3.report_isoseq_primers
    File? report_isoseq = pb_isoseq3.report_isoseq
    File? report_isoseq_mapping = pb_isoseq3.report_isoseq_mapping
    Array[File] datastore_refine = pb_isoseq3.datastore_refine
    Array[File]? datastore_cluster = pb_isoseq3.datastore_cluster
    File? collapse_abundance_joint = pb_isoseq3.collapse_abundance_joint
    #File? lima_summary = pb_isoseq3.lima_summary
    #File? lima_infer_log = pb_isoseq3.lima_infer_log
  }
}
