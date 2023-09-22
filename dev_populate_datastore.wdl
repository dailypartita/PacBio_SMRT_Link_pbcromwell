version 1.0

import "dev/tasks.wdl" as dev_tasks

workflow dev_populate_datastore {
  input {
    File eid_subread
    Int num_subreadsets = 25
    Int sleep_multiplier = 0

    Int nproc = 1
    Int max_nchunks = 100
    String log_level = "INFO"
    String? tmp_dir
  }

  call dev_tasks.populate_datastore {
    input:
      subreads = eid_subread,
      num_subreadsets = num_subreadsets,
      sleep_multiplier = sleep_multiplier,
      nproc = 1,
      log_level = log_level
  }

  output {
    File datastore = populate_datastore.datastore
  }
}
