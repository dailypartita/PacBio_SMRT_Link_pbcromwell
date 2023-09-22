version 1.0

import "wf_prepare_input.wdl"

task lima_trim {
  input {
    File reads
    File adapters
    Boolean symmetric_barcodes = true
    Boolean peek_guess = true
    Boolean write_unbarcoded = true
    Int min_score = 0
    String? other_args

    # XXX this is excessive for most use cases, but will need to be increased
    # for large asymmetrically barcoded datasets
    Int mem_gb = 2
    Int nproc
    String log_level = "INFO"
  }
  command {
    set -e
    lima \
      -j ${nproc} \
      ${true="--peek-guess" false="" peek_guess} \
      --preset ${true="HIFI-SYMMETRIC" false="HIFI-ASYMMETRIC" symmetric_barcodes} \
      ${"--min-score " + min_score} \
      ${true="--dump-removed" false="" write_unbarcoded} \
      ${other_args} \
      --ignore-xml-biosamples \
      --split-named \
      --guess-file-json lima_guesses_report.json \
      --alarms alarms.json \
      `readlink -f ${reads}` \
      `readlink -f ${adapters}` \
      trimmed.consensusreadset.xml
    sed -i 's/ (filtered)/ (trimmed)/;' trimmed.consensusreadset.xml
    #ln -s trimmed.json trimmed.datastore.json
    if [ -f "trimmed.lima.summary" ]; then
      ln -s trimmed.lima.summary trimmed.lima.summary.txt
    fi
    if [ -f "trimmed.lima.guess" ]; then
      ln -s trimmed.lima.guess trimmed.lima.guess.txt
    fi
  }
  runtime {
    cpu: nproc
    memory: "${mem_gb}GB"
  }
  output {
    File dataset = "trimmed.consensusreadset.xml"
    File? unbarcoded = "trimmed.unbarcoded.bam"
    File summary = "trimmed.lima.summary.txt"
    File? infer_log = "trimmed.lima.guess.txt"
    File? infer_report = "lima_guesses_report.json"
    File output_dir = "."
  }
}

task require_single_primer {
  input {
    File barcodeset

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 -m pbcoretools.tasks.require_single_primer \
      --log-level ${log_level} \
      `readlink -f ${barcodeset}`
  }
  runtime {
    backend: "Local"
    memory: "1GB"
  }
  output {
    # This is only necessary to force this task to run before everything
    # else
    File barcodes_out = barcodeset
  }
}

task pbreports_trim_adapters {
  input {
    File trimmed
    File untrimmed
    File adapters
    File? lima_infer_report

    Boolean test_mode = false
    Int lima_min_score = 26

    Int nproc = 1
    String log_level = "INFO"
  }

  # FIXME the order here is odd due to the way the Python code is written
  command {
    python3 \
      -m pbreports.report.trim_adapters \
      --log-level ${log_level} \
      --nproc ${nproc} \
      --report-csv trim_adapters_summary.csv \
      `readlink -f ${trimmed}` \
      trim_adapters.report.json
  }
  Int total_mem_mb = 4096 + (nproc * 512)
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File report = "trim_adapters.report.json"
    File summary_csv = "trim_adapters_summary.csv"
    Array[File?] plot_pngs = glob("*.png")
    Array[File?] plot_jsons = glob("*.json.gz")
  }
}

workflow pb_trim_adapters {
  input {
    File eid_ccs
    File eid_barcode
    Int lima_min_score = 80
    Int filter_min_qv = 20

    Int nproc
    Int max_nchunks = 0
    String log_level = "INFO"
    String? tmp_dir
  }

  # option not exposed in SMRTLink
  Boolean lima_write_unbarcoded = true
  Boolean lima_symmetric_barcodes = true

  call wf_prepare_input.prepare_input {
    input:
      dataset_xml = eid_ccs,
      filter_min_qv = filter_min_qv,
      nproc = nproc,
      log_level = log_level
  }

  call require_single_primer {
    input:
      barcodeset = eid_barcode,
      log_level = log_level
  }

  call lima_trim as lima {
    input:
      reads = prepare_input.reads_file,
      adapters = require_single_primer.barcodes_out,
      symmetric_barcodes = lima_symmetric_barcodes,
      min_score = lima_min_score,
      write_unbarcoded = lima_write_unbarcoded,
      nproc = nproc,
      log_level = log_level
  }

  call pbreports_trim_adapters {
    input:
      trimmed = lima.dataset,
      untrimmed = prepare_input.reads_file,
      adapters = eid_barcode,
      nproc = nproc,
      log_level = log_level
  }

  output {
    File trimmed = lima.dataset
    File report_trim_adapters = pbreports_trim_adapters.report
    File summary_csv = pbreports_trim_adapters.summary_csv
#    File reports_datastore = pbreports_trim_adapters.reports_datastore
    File lima_summary = lima.summary
    File? unbarcoded = lima.unbarcoded
  }
}
