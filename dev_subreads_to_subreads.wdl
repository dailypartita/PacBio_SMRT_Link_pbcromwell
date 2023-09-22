# Mock workflow that simulates demultiplexing with a single output, which
# may have a pre-defined UUID

version 1.0

import "dev/tasks.wdl" as dev_tasks

workflow dev_subreads_to_subreads {
  input {
    File eid_subread
    String use_uuid = ""

    Int nproc = 1
    Int max_nchunks = 100
    String log_level = "INFO"
    String? tmp_dir
  }

  call dev_tasks.subreads_to_subreads {
    input:
      subreads = eid_subread,
      use_uuid = use_uuid,
      nproc = nproc,
      log_level = log_level
  }

  output {
    File subreads_out = subreads_to_subreads.subreads_out
    File report_subreads = subreads_to_subreads.report_subreads
    File report_misc = subreads_to_subreads.report_misc
  }
}
