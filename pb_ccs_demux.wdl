version 1.0

import "pb_ccs.wdl" as ccs
import "pb_demux_ccs.wdl"

workflow pb_ccs_demux {
  input {
    File eid_subread
    File eid_barcode
    File? biosamples_csv
    String? biosamples_csv_str
    String? new_dataset_name
    Float ccs_min_predicted_accuracy = 0.99
    Int ccs_min_passes = 3
    Float ccs_min_snr = 2.5
    Int ccs_min_length = 10
    Int ccs_max_length = 50000
    Int filter_min_qv = 20
    Boolean lima_symmetric_barcodes = true
    Boolean lima_peek_guess = true
    Int lima_min_score = 80
    String? ccs_model_args
    Boolean pb_test_mode = false
    Boolean ccs_process_all = false
    Boolean ccs_include_kinetics = false
    # these will be set to true by Run Design
    # FIXME this should be a single option
    Boolean ccs_use_run_design_uuid = false
    Boolean use_barcode_uuids = false

    Int nproc
    Int max_nchunks = 300
    String log_level = "INFO"
    String? tmp_dir
  }

  call ccs.pb_ccs {
    input:
      eid_subread = eid_subread,
      ccs_min_predicted_accuracy = ccs_min_predicted_accuracy,
      ccs_min_passes = ccs_min_passes,
      ccs_min_snr = ccs_min_snr,
      ccs_min_length = ccs_min_length,
      ccs_max_length = ccs_max_length,
      ccs_use_run_design_uuid = ccs_use_run_design_uuid,
      ccs_model_args = ccs_model_args,
      ccs_process_all = ccs_process_all,
      ccs_include_kinetics = ccs_include_kinetics,
      nproc = nproc,
      max_nchunks = max_nchunks,
      log_level = log_level,
      tmp_dir = tmp_dir
  }

  call pb_demux_ccs.pb_demux_ccs {
    input:
      eid_ccs = pb_ccs.ccsxml,
      eid_barcode = eid_barcode,
      biosamples_csv = biosamples_csv,
      biosamples_csv_str = biosamples_csv_str,
      lima_symmetric_barcodes = lima_symmetric_barcodes,
      lima_peek_guess = lima_peek_guess,
      lima_min_score = lima_min_score,
      new_dataset_name = new_dataset_name,
      use_barcode_uuids = use_barcode_uuids,
      filter_min_qv = filter_min_qv,
      log_level = log_level
  }

  output {
    File ccsxml = pb_ccs.ccsxml
    File report_ccs2 = pb_ccs.report_ccs2
    File report_ccs_processing = pb_ccs.report_ccs_processing
    File ccs_csv = pb_ccs.ccs_csv
    File ccs_zmws = pb_ccs.ccs_zmws
    File bam_datastore = pb_ccs.bam_datastore
    File fasta_datastore = pb_ccs.fasta_datastore
    File fastq_datastore = pb_ccs.fastq_datastore
    File barcoded_reads = pb_demux_ccs.barcoded_reads
    File report_barcodes = pb_demux_ccs.report_barcodes
    File summary_csv = pb_demux_ccs.summary_csv
    File lima_summary = pb_demux_ccs.lima_summary
    File? lima_infer_log = pb_demux_ccs.lima_infer_log
    File? barcoded_reports_datastore = pb_demux_ccs.barcoded_reports_datastore
    File demuxed_files_datastore = pb_demux_ccs.demuxed_files_datastore
    File? unbarcoded = pb_demux_ccs.unbarcoded
    String? reads_bam = pb_ccs.reads_bam
    File demultiplexing_files = pb_demux_ccs.demultiplexing_files
    File fastx_files = pb_demux_ccs.fastx_files
#    File? user_barcodes_csv = create_input_dataset.user_barcodes_csv
  }
}
