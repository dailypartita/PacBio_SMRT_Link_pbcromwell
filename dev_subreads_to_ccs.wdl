# Mock workflow that simulates generating a CCS dataset, which may have a
# pre-defined UUID

version 1.0

import "dev/tasks.wdl" as dev_tasks

workflow dev_subreads_to_ccs {
  input {
    File eid_subread
    String use_uuid = ""

    Int nproc = 1
    Int max_nchunks = 100
    String log_level = "INFO"
    String? tmp_dir
  }

  call dev_tasks.subreads_to_ccs {
    input:
      subreads = eid_subread,
      use_uuid = use_uuid,
      nproc = nproc,
      log_level = log_level
  }

  output {
    File ccs_out = subreads_to_ccs.ccs_out
    File report_ccs = subreads_to_ccs.report_ccs
  }
}
