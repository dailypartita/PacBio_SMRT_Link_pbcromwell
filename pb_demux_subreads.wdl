# LEGACY SUBREAD DEMULTIPLEXING

version 1.0

import "tasks/memory.wdl" as memory

task lima_subreads {
  input {
    File reads
    File barcodes
    Boolean symmetric_barcodes = true
    Boolean peek_guess = true
    Boolean write_unbarcoded = true
    Int min_score = 26
    Boolean? ignore_biosamples
    String? other_args
    Boolean generate_missing_datasets = false

    # XXX this is excessive for most use cases, but will need to be increased
    # for large asymmetrically barcoded datasets
    Int mem_gb = 2
    Int nproc
    String log_level = "INFO"
  }
  command {
    lima \
      -j ${nproc} \
      ${true="--peek-guess" false="" peek_guess} \
      ${true="--same" false="--different" symmetric_barcodes} \
      ${"--min-score " + min_score} \
      ${true="--dump-removed" false="" write_unbarcoded} \
      ${true="--ignore-biosamples" false="" ignore_biosamples} \
      ${true="--output-missing-pairs" false="" generate_missing_datasets} \
      ${other_args} \
      --split-bam-named \
      --guess-file-json lima_guesses_report.json \
      --alarms alarms.json \
      `readlink -f ${reads}` \
      `readlink -f ${barcodes}` \
      demultiplex.json
    if [ -f "demultiplex.lima.summary" ]; then
      mv demultiplex.lima.summary demultiplex.lima.summary.txt
    fi
    if [ -f "demultiplex.lima.guess" ]; then
      mv demultiplex.lima.guess demultiplex.lima.guess.txt
    fi
  }
  runtime {
    cpu: nproc
    memory: "${mem_gb}GB"
  }
  output {
    File datastore = "demultiplex.json"
    File? unbarcoded = "demultiplex.unbarcoded.bam"
    File summary = "demultiplex.lima.summary.txt"
    File? infer_log = "demultiplex.lima.guess.txt"
    File? infer_report = "lima_guesses_report.json"
    File output_dir = "."
  }
}

task update_barcoded_sample_metadata {
  input {
    File lima_datastore
    File input_reads
    File barcodes
    Int min_bq_filter = 26
    Boolean isoseq_mode = false
    Boolean use_barcode_uuids = false

    Int nproc
    Int mem_gb
    String log_level = "DEBUG"
  }
  command {
    python3 -m pbcoretools.tasks.update_barcoded_sample_metadata \
      --log-level ${log_level} \
      --nproc ${nproc} \
      ${true="--isoseq-mode" false="" isoseq_mode} \
      ${true="--use-barcode-uuids" false="" use_barcode_uuids} \
      ${"--min-bq-filter " + min_bq_filter} \
      `readlink -f ${input_reads}` \
      `readlink -f ${lima_datastore}` \
      `readlink -f ${barcodes}` \
      metadata_updated.datastore.json
  }
  runtime {
    cpu: nproc
    memory: "${mem_gb}GB"
  }
  output {
    File datastore = "metadata_updated.datastore.json"
    File? empty_files_datastore = "missing_barcodes.datastore.json"
  }
}

task create_input_dataset {
  input {
    File input_reads
    String dataset_name
    String dataset_ext = ".subreadset.xml"
    File? biosamples_csv_file
    String? biosamples_csv_str

    Int nproc = 1
    String log_level = "DEBUG"
  }
  String sample_file_name = "User_Input_Barcoded_Samples.csv"
  command <<<
    set -e
    extra_args="--biosamples-csv ~{sample_file_name}"
    if [ -f "~{biosamples_csv_file}" ]; then
      echo "biosamples_csv_file exists, copying file"
      cp ~{biosamples_csv_file} ~{sample_file_name}
    elif [ ! -z "~{biosamples_csv_str}" ]; then
      echo "biosamples_csv_str is defined, creating file"
      echo "~{biosamples_csv_str}" > ~{sample_file_name}
    else
      echo "No biosamples CSV input included"
      extra_args=""
    fi
    python3 -m pbcoretools.tasks.reparent_dataset \
      --log-level ~{log_level} \
      --suffix "(copy)" \
      `readlink -f ~{input_reads}` \
      "~{dataset_name}" \
      unbarcoded~{dataset_ext} \
      $extra_args
  >>>
  runtime {
    cpu: 1
    memory: "1GB"
  }
  output {
    File new_parent = "unbarcoded${dataset_ext}"
    File? user_barcodes_csv = "${sample_file_name}"
  }
}

task datastore_to_merged_dataset {
  input {
    File datastore
    String dataset_ext = ".subreadset.xml"

    Int nproc = 1
    String log_level = "DEBUG"
    Int index_memory_gb = 4
  }
  command {
    python3 -c 'import sys; import os.path ; from pbcommand.models import DataStore ; ds = DataStore.load_from_json(os.path.realpath(sys.argv[1])) ; print("\n".join([f.path for f in ds.files.values()]))' ${datastore} > datasets.fofn
    dataset \
      --strict \
      create \
      --unique-collections \
      --no-sub-datasets \
      merged${dataset_ext} \
      datasets.fofn
  }
  runtime {
    cpu: 1
    # this might be excessive?
    memory: "${index_memory_gb}GB"
  }
  output {
    File mergedxml = "merged${dataset_ext}"
  }
}

task estimate_lima_memory {
  input {
    File data
    File barcodes
    Boolean symmetric_barcodes

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 -m pbcoretools.tasks.memory.estimate_lima_memory \
      --log-level ${log_level} \
      `readlink -f ${barcodes}` \
      `readlink -f ${data}` \
      ${true="--symmetric" false="--asymmetric" symmetric_barcodes}
  }
  runtime {
    cpu: 1
    # this is probably excessive but we need to test with bloated dataset XMLs
    memory: "200MB"
    backend: "Local"
  }
  output {
    File lima_mem_gb_txt = "lima_mem_gb.txt"
  }
}

task pbreports_barcode_subreads {
  input {
    File barcoded
    File unbarcoded
    File barcodes
    File? lima_infer_report

    Boolean test_mode = false
    Int lima_min_score = 26

    Int nproc = 1
    String log_level = "INFO"
  }

  # FIXME the order here is odd due to the way the Python code is written
  command {
    python3 \
      -m pbreports.report.barcode_subreads \
      --log-level ${log_level} \
      --nproc ${nproc} \
      ${true="--test-mode" false="" test_mode} \
      ${"--lima-min-score " + lima_min_score} \
      ${"--guess-file-json " + lima_infer_report} \
      `readlink -f ${barcoded}` \
      `readlink -f ${unbarcoded}` \
      `readlink -f ${barcodes}` \
      barcode.report.json \
      barcode_summary.csv \
      barcoded_reports.datastore.json
  }
  Int total_mem_mb = 4096 + (nproc * 512)
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File report = "barcode.report.json"
    File summary_csv = "barcode_summary.csv"
    File? reports_datastore = "barcoded_reports.datastore.json"
    Array[File?] plot_pngs = glob("*.png")
    Array[File?] plot_jsons = glob("*.json.gz")
  }
}

workflow pb_demux_subreads {
  input {
    File eid_subread
    File eid_barcode
    File? biosamples_csv
    String? biosamples_csv_str
    String new_dataset_name = ""
    Boolean lima_symmetric_barcodes = true
    Boolean lima_peek_guess = true
    Int lima_min_score = 26
    Boolean lima_write_unbarcoded = true
    Boolean use_barcode_uuids = false
    Boolean pb_test_mode = false
    String? lima_overrides

    Int nproc = 8
    # unused
    Int max_nchunks = 0
    String log_level = "INFO"
    String? tmp_dir
  }

  call memory.get_input_sizes {
    input:
      bam_dataset = eid_subread,
      fasta_dataset = eid_barcode,
      log_level = log_level
  }

  if (!use_barcode_uuids) {
    call create_input_dataset {
      input:
        input_reads = eid_subread,
        biosamples_csv_file = biosamples_csv,
        biosamples_csv_str = biosamples_csv_str,
        dataset_name = new_dataset_name,
        dataset_ext = ".subreadset.xml",
        nproc = nproc,
        log_level = log_level
    }
  }

  Array[File?] all_parents = [create_input_dataset.new_parent, eid_subread]
  File parent_dataset = select_first(all_parents)

  call estimate_lima_memory {
    input:
      data = parent_dataset,
      barcodes = eid_barcode,
      symmetric_barcodes = lima_symmetric_barcodes,
      log_level = log_level
  }
  Int lima_mem_gb = read_int(estimate_lima_memory.lima_mem_gb_txt)

  call lima_subreads as lima {
    input:
      reads = parent_dataset,
      barcodes = eid_barcode,
      symmetric_barcodes = lima_symmetric_barcodes,
      peek_guess = lima_peek_guess,
      min_score = lima_min_score,
      write_unbarcoded = lima_write_unbarcoded,
      other_args = lima_overrides,
      generate_missing_datasets = use_barcode_uuids,
      mem_gb = lima_mem_gb,
      nproc = nproc,
      log_level = log_level
  }

  call update_barcoded_sample_metadata {
    input:
      lima_datastore = lima.datastore,
      input_reads = parent_dataset,
      barcodes = eid_barcode,
      isoseq_mode = false,
      use_barcode_uuids = use_barcode_uuids,
      nproc = nproc,
      mem_gb = get_input_sizes.index_memory_gb,
      log_level = log_level
  }

  call pbreports_barcode_subreads as barcode_report {
    input:
      barcoded = update_barcoded_sample_metadata.datastore,
      unbarcoded = parent_dataset,
      barcodes = eid_barcode,
      lima_infer_report = lima.infer_report,
      test_mode = pb_test_mode,
      lima_min_score = lima_min_score,
      nproc = nproc,
      log_level = log_level
  }

  output {
    File? input_reads = create_input_dataset.new_parent
    File barcoded_reads = update_barcoded_sample_metadata.datastore
    File report_barcodes = barcode_report.report
    File summary_csv = barcode_report.summary_csv
    File lima_summary = lima.summary
    File? lima_infer_log = lima.infer_log
    File? barcoded_reports_datastore = barcode_report.reports_datastore
    File? unbarcoded = lima.unbarcoded
    File? empty_files_datastore = update_barcoded_sample_metadata.empty_files_datastore
    File demultiplexing_files = lima.output_dir
    File? user_barcodes_csv = create_input_dataset.user_barcodes_csv
  }
}
