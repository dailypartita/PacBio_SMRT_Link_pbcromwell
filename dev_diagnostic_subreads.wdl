version 1.0

import "tasks/pbreports.wdl" as pbreports
import "dev/tasks.wdl" as dev_tasks
import "tasks/sl_tasks.wdl" as sl_tasks

task echo_txt {
  input {
    File? optional_txt_input
  }
  command {
    cat ${optional_txt_input} > "hello.txt"
  }
  runtime {
    cpu: 1
  }
  output {
    File txt_out = "hello.txt"
  }
}

task echo_genome_size {
  input {
    Int? genome_size
  }
  command {
    echo "${genome_size}" > genome_size.txt
  }
  runtime {
    cpu: 1
  }
  output {
    File txt_out = "genome_size.txt"
  }
}

workflow dev_diagnostic_subreads {
  input {
    File eid_subread
    File? eid_ref_dataset
    File? optional_txt_input
    Int? exit_code
    # this is just for testing int overflow
    Int? genome_size
    Int sleep_time = 0
    Boolean emit_warn_alarm = false
    Boolean emit_error_alarm = false

    Int nproc = 1
    Int max_nchunks = 100
    String log_level = "INFO"
    String? tmp_dir
  }

  call dev_tasks.dump_environment {
    input:
      sleep_time = sleep_time,
      nproc = 1,
      log_level = "INFO"
  }

  call pbreports.subread_stats {
    input:
      subreads = eid_subread,
      nproc = nproc,
      log_level = log_level
  }

  call dev_tasks.samples_report {
    input:
      dataset = eid_subread,
      nproc = nproc,
      log_level = log_level
  }

  call sl_tasks.simple_dataset_report {
    input:
      dataset = eid_subread,
      nproc = nproc,
      log_level = log_level
  }

  if (defined(exit_code)) {
    call dev_tasks.dev_failing_task {
      input:
        exit_code = exit_code,
        sleep_time = sleep_time,
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
        dataset = eid_subread
    }
  }

  if (emit_error_alarm) {
    call dev_tasks.raise_alarm {
      input:
        dataset = eid_subread
    }
  }

  if (defined(optional_txt_input)) {
    scatter (x in [1, 2]) {
      call echo_txt {
        input:
          optional_txt_input = optional_txt_input
      }
    }
  }

  if (defined(genome_size)) {
    call echo_genome_size {
      input:
        genome_size = genome_size
    }
  }

  output {
    File report_subreads = subread_stats.report
    File? report_dataset = simple_dataset_report.report
    File report_samples = samples_report.report
    File? report_reference = reference_report.report
    Array[File]? txt_out = echo_txt.txt_out
    File? genome_txt = echo_genome_size.txt_out
  }
}
