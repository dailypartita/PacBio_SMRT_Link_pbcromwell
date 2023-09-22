# Combined demux + Iso-Seq analysis for internal testing

version 1.0

import "pb_isoseq3.wdl"
import "pb_demux_ccs.wdl"

workflow dev_demux_isoseq {
  input {
    File eid_ccs
    File eid_barcode
    File eid_barcode_2
    File? biosamples_csv
    Boolean lima_symmetric_barcodes = true
    Boolean run_clustering = true
    Boolean cluster_separately = false

    Int nproc = 1
    String log_level = "INFO"
    String? tmp_dir
    Int max_nchunks = 1
  }

  call pb_demux_ccs.pb_demux_ccs {
    input:
      eid_ccs = eid_ccs,
      eid_barcode = eid_barcode,
      biosamples_csv = biosamples_csv,
      lima_symmetric_barcodes = lima_symmetric_barcodes,
      nproc = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir,
      max_nchunks = max_nchunks
  }

  call pb_isoseq3.pb_isoseq3 {
    input:
      eid_ccs = pb_demux_ccs.barcoded_reads,
      eid_barcode = eid_barcode_2,
      run_clustering = run_clustering,
      cluster_separately = cluster_separately,
      nproc = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir,
      max_nchunks = max_nchunks
  }

  output {
    File barcoded_reads = pb_demux_ccs.barcoded_reads
    File demuxed_files_datastore = pb_demux_ccs.demuxed_files_datastore
    File report_barcodes = pb_demux_ccs.report_barcodes
    File report_isoseq_primers = pb_isoseq3.report_isoseq_primers
    File? report_isoseq = pb_isoseq3.report_isoseq
    File? report_isoseq_mapping = pb_isoseq3.report_isoseq_mapping
  }
}
