# workflow version of SMRT Link import-fasta-barcodes job

version 1.0

import "tasks/pbcoretools.wdl" as pbcoretools
import "tasks/sl_tasks.wdl" as sl_tasks

workflow sl_import_fasta_barcodes {
  input {
    File eid_fasta
    String dataset_name

    Int nproc = 1
    String log_level = "INFO"
    String? tmp_dir
  }

  call pbcoretools.pbvalidate {
    input:
      data = eid_fasta,
      quick_mode = false,
      log_level = log_level
  }

  if (!defined(pbvalidate.alarms)) {
    call sl_tasks.fasta_to_barcodes {
      input:
        barcodes_fasta = eid_fasta,
        dataset_name = dataset_name,
        log_level = log_level
    }

    # This is here mostly just for testing purposes
    call sl_tasks.simple_dataset_report {
      input:
        dataset = fasta_to_barcodes.barcodeset,
        skip_counts = true,
        nproc = 1,
        log_level = log_level
    }
  }

  output {
    File? barcodeset = fasta_to_barcodes.barcodeset
    File? report_dataset = simple_dataset_report.report
  }
}
