# Mock workflow that simulates CCS demultiplexing with a single output, which
# may have a pre-defined UUID

version 1.0

import "dev/tasks.wdl" as dev_tasks

workflow dev_ccs_to_ccs {
  input {
    File eid_ccs
    String use_uuid = ""

    Int nproc = 1
    Int max_nchunks = 100
    String log_level = "INFO"
    String? tmp_dir
  }

  call dev_tasks.ccs_to_ccs {
    input:
      ccsxml = eid_ccs,
      use_uuid = use_uuid,
      nproc = nproc,
      log_level = log_level
  }

  output {
    File ccs_out = ccs_to_ccs.ccs_out
    File report_ccs = ccs_to_ccs.report_ccs
  }
}
