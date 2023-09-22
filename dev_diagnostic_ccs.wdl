version 1.0

import "tasks/pbreports.wdl" as pbreports
import "pb_ccs.wdl"
import "dev/tasks.wdl" as dev_tasks
import "tasks/sl_tasks.wdl" as sl_tasks

workflow dev_diagnostic_ccs {
  input {
    File eid_ccs
    File? eid_ref_dataset
    Int? exit_code
    Boolean emit_warn_alarm = false
    Boolean emit_error_alarm = false

    Int nproc = 1
    Int max_nchunks = 100
    String log_level = "INFO"
    String? tmp_dir
  }

  call dev_tasks.dump_environment {
    input:
      nproc = 1,
      log_level = "INFO"
  }

  call pb_ccs.pbreports_ccs2 as ccs2_report {
    input:
      ccsxml = eid_ccs,
      nproc = nproc,
      log_level = log_level
  }

  call dev_tasks.samples_report {
    input:
      dataset = eid_ccs,
      nproc = nproc,
      log_level = log_level
  }

  call sl_tasks.simple_dataset_report {
    input:
      dataset = eid_ccs,
      nproc = nproc,
      log_level = log_level
  }

  if (defined(exit_code)) {
    call dev_tasks.dev_failing_task {
      input:
        exit_code = exit_code,
        nproc = nproc,
        log_level = log_level
    }
  }

  if (defined(eid_ref_dataset)) {

    call sl_tasks.simple_dataset_report as reference_report {
      input:
        dataset = eid_ref_dataset,
        nproc = nproc,
        log_level = log_level
    }
  }

  if (emit_warn_alarm) {
    call dev_tasks.emit_alarm {
      input:
        dataset = eid_ccs
    }
  }

  if (emit_error_alarm) {
    call dev_tasks.raise_alarm {
      input:
        dataset = eid_ccs
    }
  }

  output {
    File report_ccs2 = ccs2_report.report
    File? report_dataset = simple_dataset_report.report
    File report_samples = samples_report.report
    File? report_reference = reference_report.report
  }
}
