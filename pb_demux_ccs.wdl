# CCS demultiplexing, synced to on-instrument implementation when run with
# use_barcode_uuids=true

version 1.0

task lima {
  input {
    File reads
    File barcodes
    String prefix
    Boolean symmetric_barcodes = true
    Boolean peek_guess = true
    Boolean use_barcode_uuids = false
    Boolean write_unbarcoded = true
    Int? min_score
    Int? min_qv
    String? other_args
    String? new_dataset_name_arg
    File? biosamples_csv

    # XXX this is excessive for most use cases, but will need to be increased
    # for large asymmetrically barcoded datasets
    Int mem_gb = 2
    Int base_memory_mb = 0
    Int nproc
    String log_level = "INFO"
  }
  Boolean ignore_biosamples = defined(biosamples_csv) && !use_barcode_uuids
  command {
    set -e
    lima \
      -j ${nproc} \
      ${true="--peek-guess" false="" peek_guess} \
      --hifi-preset ${true="SYMMETRIC-ADAPTERS" false="ASYMMETRIC" symmetric_barcodes} \
      ${"--min-score " + min_score} \
      ${"--min-qv " + min_qv} \
      ${true="--dump-removed" false="" write_unbarcoded} \
      ${true="--reuse-uuids" false="" use_barcode_uuids} \
      ${true="--output-missing-pairs" false="" use_barcode_uuids} \
      ${true="--ignore-xml-biosamples" false="" ignore_biosamples} \
      ${new_dataset_name_arg} \
      ${"--biosample-csv " + biosamples_csv} \
      ${other_args} \
      --split-bam-named \
      --split-subdirs \
      --guess-file-json lima_guesses_report.json \
      --alarms alarms.json \
      `readlink -f ${reads}` \
      `readlink -f ${barcodes}` \
      ${prefix}.consensusreadset.xml
    if [ -f "${prefix}.lima.summary" ]; then
      ln -s ${prefix}.lima.summary ${prefix}.lima.summary.txt
    fi
    if [ -f "${prefix}.lima.guess" ]; then
      ln -s ${prefix}.lima.guess ${prefix}.lima.guess.txt
    fi
  }
  Int total_mem_mb = mem_gb * 1024 + base_memory_mb
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File dataset_xml = "${prefix}.consensusreadset.xml"
    File? unbarcoded = "${prefix}.unbarcoded.bam"
    File summary = "${prefix}.lima.summary.txt"
    File? infer_log = "${prefix}.lima.guess.txt"
    File? infer_report = "lima_guesses_report.json"
    File output_dir = "."
    File? alarms = "alarms.json"
  }
}

task setup_demux {
  input {
    File data
    File barcodes
    Boolean symmetric_barcodes
    File? biosamples_csv
    String? biosamples_csv_str
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  # either way we will provide a file to the Python tool
  Boolean have_csv = ((defined(biosamples_csv) && (biosamples_csv != "")) ||
                      (defined(biosamples_csv_str) && (biosamples_csv_str != "")))
  command {
    set -e
    if [ ! -z "${biosamples_csv_str}" ]; then
      echo "${biosamples_csv_str}" > biosamples_in.csv
    elif [ ! -z "${biosamples_csv}" ] && [ -e "${biosamples_csv}" ]; then
      ln -s "${biosamples_csv}" biosamples_in.csv
    fi
    python3 -m pbcoretools.tasks.memory.estimate_lima_memory \
      --log-level ${log_level} \
      `readlink -f ${barcodes}` \
      `readlink -f ${data}` \
      ${true="--symmetric" false="--asymmetric" symmetric_barcodes}
    python3 -m pbcoretools.tasks.generate_barcode_input \
      --log-level ${log_level} \
      ${true="--biosamples-csv biosamples_in.csv " false="" have_csv} \
      -o lima_biosamples.csv \
      `readlink -f ${data}`
    python3 -m pbcoretools.tasks.get_demux_file_prefix \
      --log-level ${log_level} \
      --output-file prefix.txt \
      `readlink -f ${data}`
  }
  # this is probably excessive but we need to test with bloated dataset XMLs
  Int MIN_MEM_MB = 200
  Int total_mem_mb = base_memory_mb + MIN_MEM_MB
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
    backend: "Local"
  }
  output {
    File lima_mem_gb_txt = "lima_mem_gb.txt"
    File? lima_biosamples_csv = "lima_biosamples.csv"
    String prefix = read_string("prefix.txt")
  }
}

task auto_ccs_outputs_barcoded {
  input {
    File demultiplexed
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }

  command {
    python3 \
      -m pbcoretools.tasks.auto_ccs_outputs_barcoded \
      --log-level ${log_level} \
      --nproc ${nproc} \
      `readlink -f ${demultiplexed}` \
      ccs_demultiplexed_outputs.datastore.json
  }
  # also probably an overestimate
  Int total_mem_mb = base_memory_mb + 1024 * (nproc + 1)
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File datastore = "ccs_demultiplexed_outputs.datastore.json"
    File output_dir = "."
  }
}

task pbreports_barcode_ccs {
  input {
    File barcoded
    String prefix
    Int lima_min_score = 80
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 \
      -m pbreports.report.barcode \
      --log-level ${log_level} \
      --nproc ${nproc} \
      --report-csv ${prefix}.barcodes_summary.csv \
      --per-barcode-reports per_barcode_reports.datastore.json \
      `readlink -f ${barcoded}` \
      ${prefix}.barcodes.report.json
  }
  Int total_mem_mb = base_memory_mb + 4096 * (nproc + 1)
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File report = "${prefix}.barcodes.report.json"
    File summary_csv = "${prefix}.barcodes_summary.csv"
    File? reports_datastore = "per_barcode_reports.datastore.json"
    Array[File?] plot_pngs = glob("*.png")
    Array[File?] plot_jsons = glob("*.json.gz")
  }
}

workflow pb_demux_ccs {
  input {
    File eid_ccs
    File eid_barcode
    File? biosamples_csv
    String? biosamples_csv_str
    String? new_dataset_name
    Int filter_min_qv = 20
    Boolean lima_symmetric_barcodes = true
    Boolean lima_peek_guess = true
    Int lima_min_score = 80
    Boolean use_barcode_uuids = false
    String? lima_overrides
    Int add_memory_mb = 0

    Int nproc = 8
    Int max_nchunks = 0
    String log_level = "INFO"
    String? tmp_dir
  }
  # this is required by the new report
  Boolean lima_write_unbarcoded = true

  call setup_demux {
    input:
      data = eid_ccs,
      barcodes = eid_barcode,
      biosamples_csv = biosamples_csv,
      biosamples_csv_str = biosamples_csv_str,
      symmetric_barcodes = lima_symmetric_barcodes,
      base_memory_mb = add_memory_mb,
      log_level = log_level
  }
  Int lima_mem_gb = read_int(setup_demux.lima_mem_gb_txt)

  # workaround for services behavior (null string becomes empty string)
  if (defined(new_dataset_name) && new_dataset_name != "") {
    String? new_dataset_name_arg = "--dataset-name '${new_dataset_name}'"
  }

  call lima {
    input:
      reads = eid_ccs,
      barcodes = eid_barcode,
      prefix = setup_demux.prefix,
      symmetric_barcodes = lima_symmetric_barcodes,
      peek_guess = lima_peek_guess,
      min_score = lima_min_score,
      min_qv = filter_min_qv,
      write_unbarcoded = lima_write_unbarcoded,
      new_dataset_name_arg = new_dataset_name_arg,
      biosamples_csv = setup_demux.lima_biosamples_csv,
      use_barcode_uuids = use_barcode_uuids,
      other_args = lima_overrides,
      mem_gb = lima_mem_gb,
      base_memory_mb = add_memory_mb,
      nproc = nproc,
      log_level = log_level
  }

  Int MAX_NPROC_REPORTS = 8
  Int nproc_reports = if (nproc > MAX_NPROC_REPORTS) then MAX_NPROC_REPORTS else nproc
  call pbreports_barcode_ccs as barcode_report {
    input:
      barcoded = lima.dataset_xml,
      prefix = setup_demux.prefix,
      lima_min_score = lima_min_score,
      base_memory_mb = add_memory_mb,
      nproc = nproc_reports,
      log_level = log_level
  }

  call auto_ccs_outputs_barcoded {
    input:
      demultiplexed = lima.dataset_xml,
      base_memory_mb = add_memory_mb,
      nproc = nproc,
      log_level = log_level
  }

  output {
    File barcoded_reads = lima.dataset_xml
    File report_barcodes = barcode_report.report
    File summary_csv = barcode_report.summary_csv
    File lima_summary = lima.summary
    File? lima_infer_log = lima.infer_log
    File? barcoded_reports_datastore = barcode_report.reports_datastore
    File demuxed_files_datastore = auto_ccs_outputs_barcoded.datastore
    File? unbarcoded = lima.unbarcoded
    File demultiplexing_files = lima.output_dir
    File fastx_files = auto_ccs_outputs_barcoded.output_dir
    File? user_barcodes_csv = setup_demux.lima_biosamples_csv
  }
}
