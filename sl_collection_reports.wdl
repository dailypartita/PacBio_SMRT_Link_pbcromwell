# Standalone report generation for running on cluster only.

version 1.0

import "sl_dataset_reports.wdl"

# TODO this should eventually migrate to pb_demux_ccs.wdl
task pbreports_barcode {
  input {
    File barcoded

    Int nproc = 1
    String log_level = "INFO"
  }
  String reports_datastore_fn = "per_barcode_reports.datastore.json"
  command {
    python3 \
      -m pbreports.report.barcode \
      --log-level ${log_level} \
      --log-file barcode_report.log \
      --nproc ${nproc} \
      --per-barcode-reports ${reports_datastore_fn} \
      "`readlink -f '${barcoded}'`" \
      barcode.report.json
  }
  Int base_memory_mb = 16384
  Int total_mem_mb = base_memory_mb + (nproc * 512)
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File report = "barcode.report.json"
    #File summary_csv = "barcode_ccs_summary.csv"
    File? reports_datastore = "${reports_datastore_fn}"
    Array[File?] plot_pngs = glob("*.png")
    Array[File?] plot_jsons = glob("*.json.gz")
  }
}

workflow sl_collection_reports {
  input {
    File? eid_ccs
    File? eid_subread
    Boolean has_barcodes = false
    # FIXME all of this is still run on the head node right now
    Boolean make_runqc_reports = false

    Int nproc = 1
    String log_level = "INFO"
    String? tmp_dir
    Int max_nchunks = 1
  }
  Array[File?] datasets_ = [eid_subread, eid_ccs]
  File dataset = select_first(datasets_)

  if (has_barcodes && defined(eid_ccs)) {
    Int NPROC_MAX = 8
    Int nproc_final = if (nproc > NPROC_MAX) then NPROC_MAX else nproc
    call pbreports_barcode {
      input:
        barcoded = dataset,
        nproc = nproc_final,
        log_level = log_level
    }
  }

  if (make_runqc_reports) {
    call sl_dataset_reports.import_dataset_reports as collection_reports {
      input:
        dataset_xml = dataset,
        nproc = nproc,
        log_level = log_level
    }
  }

  output {
    File? report_raw_data = collection_reports.report_raw_data
    File? report_adapters = collection_reports.report_adapter
    File? report_loading = collection_reports.report_loading
    File? report_control = collection_reports.report_control
    File? report_subread_stats = collection_reports.report_subread_stats
    File? report_ccs2 = collection_reports.report_ccs2
    File? report_barcode = pbreports_barcode.report
    File? barcoded_reports_datastore = pbreports_barcode.reports_datastore
    File? report_detect_cpg_methyl = collection_reports.report_detect_cpg_methyl
#    File? barcode_summary_csv = pbreports_barcode.summary_csv
    File dataset_in = dataset
  }
}
