# CCS BAM/FASTX export workflow.  Note that the behavior here is significantly
# different from the subreads version, as it uses tasks from the CCS workflow
# instead of calling bam2fastx directly.

version 1.0

import "wf_prepare_input.wdl"
import "tasks/pbcoretools.wdl" as pbcoretools

workflow pb_export_ccs {
  input {
    File eid_ccs
    Boolean output_fasta = true
    Boolean output_bam = false
    String dataset_filters = ""
    # Phred score, so default is Q20.  Outputs will be 
    Int filter_min_qv = 20

    # ignored
    Int max_nchunks = 0
    Int nproc = 1
    String log_level = "INFO"
    String? tmp_dir
  }

  call wf_prepare_input.prepare_input {
    input:
      dataset_xml = eid_ccs,
      dataset_filters = dataset_filters,
      filter_min_qv = filter_min_qv,
      nproc = nproc,
      log_level = log_level
  }

  call pbcoretools.auto_ccs_outputs as export_fastq {
    input:
      ccsxml = prepare_input.reads_file,
      min_qv = filter_min_qv,
      mode = "fastq",
      nproc = nproc,
      log_level = log_level
  }

  if (output_fasta) {
    call pbcoretools.auto_ccs_outputs as export_fasta {
      input:
        ccsxml = prepare_input.reads_file,
        min_qv = filter_min_qv,
        mode = "fasta",
        nproc = nproc,
        log_level = log_level
    }
  }

  if (output_bam) {
    call pbcoretools.auto_ccs_outputs as export_bam {
      input:
        ccsxml = prepare_input.reads_file,
        min_qv = filter_min_qv,
        mode = "consolidate",
        nproc = nproc,
        log_level = log_level
    }
  }

  output {
    File fastq_datastore = export_fastq.datastore
    File? fasta_datastore = export_fasta.datastore
    File? bam_datastore = export_bam.datastore
  }
}
