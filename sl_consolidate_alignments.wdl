# FIXME this is obsolete now

version 1.0

import "tasks/pbcoretools.wdl" as pbcoretools

workflow sl_consolidate_alignments {
  input {
    File? eid_alignment
    File? eid_ccs_alignment

    Int nproc = 1
    Int max_nchunks = 1
    String log_level = "INFO"
    String? tmp_dir
  }

  # workaround for input dataset type polymorphism
  Array[File?] datasets_ = [eid_alignment, eid_ccs_alignment]
  File mapped_dataset = select_first(datasets_)

  call pbcoretools.auto_consolidate_alignments {
    input:
      mapped = mapped_dataset,
      force_consolidate = true,
      tmp_dir = tmp_dir,
      nproc = nproc,
      log_level = log_level
  }

  output {
    # these will always be generated because force_consolidate = true above,
    # but need to be optionally typed because they're optional in the task
    File? mapped_bam_datastore = auto_consolidate_alignments.datastore
  }
}
