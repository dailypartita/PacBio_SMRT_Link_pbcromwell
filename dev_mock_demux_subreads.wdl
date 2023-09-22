# Mock demultiplexing workflow for testing scaling and end-to-end integration

version 1.0

import "dev/tasks.wdl" as dev_tasks

workflow dev_mock_demux_subreads {
  input {
    File eid_subread
    File eid_barcode
    Boolean lima_symmetric_barcodes = true
    Boolean lima_peek_guess = true
    Int lima_min_score = 26
    Boolean use_barcode_uuids = false

    Int nproc = 1
    Int max_nchunks = 100
    String log_level = "INFO"
    String? tmp_dir
  }

  call dev_tasks.mock_demux {
    input:
      dataset = eid_subread,
      dataset_ext = ".subreadset.xml",
      log_level = log_level,
      nproc = 1
  }

  call dev_tasks.mock_update_barcoded_sample_metadata {
    input:
      lima_datastore = mock_demux.datastore,
      input_reads = eid_subread,
      barcodes = eid_barcode,
      use_barcode_uuids = use_barcode_uuids,
      nproc = 1,
      log_level = log_level
  }

  output {
    File datastore = mock_update_barcoded_sample_metadata.datastore
  }
}
