# workflow version of SMRT Link import-fasta job

version 1.0

import "tasks/pbcoretools.wdl" as pbcoretools
import "tasks/sl_tasks.wdl" as sl_tasks

workflow sl_import_fasta {
  input {
    File eid_fasta
    String reference_name
    String organism = "unknown"
    String ploidy = "haploid"

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
    call sl_tasks.fasta_to_reference {
      input:
        reference_fasta = eid_fasta,
        reference_name = reference_name,
        organism = organism,
        ploidy = ploidy,
        log_level = log_level
    }

    # This is here mostly just for testing purposes
    call sl_tasks.simple_dataset_report {
      input:
        dataset = fasta_to_reference.referenceset,
        skip_counts = true,
        nproc = 1,
        log_level = log_level
    }
  }

  output {
    File? referenceset = fasta_to_reference.referenceset
    File? report_dataset = simple_dataset_report.report
  }
}
