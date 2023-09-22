version 1.0

import "tasks/pbcoretools.wdl" as pbcoretools
import "tasks/pbreports.wdl" as pbreports
import "tasks/chunking.wdl" as chunking
import "tasks/memory.wdl" as memory
import "pb_detect_methyl.wdl"

task ccs {
  input {
    File subreads
    Boolean polishing = true
    Int min_passes = 3
    Float min_read_score = 0.65
    Float min_predicted_accuracy = 0.9
    Float min_snr = 2.5
    Int min_length = 10
    Int max_length = 50000
    # this could be either of --modelPath or --modelSpec with appropriate
    # argument; we need to pass the arguments this way due to smrtlink's
    # lack of nullables
    String? model_args
    String? draft_mode_args
    String? other_args
    Boolean? disable_heuristics = false
    Boolean by_strand = false
    Boolean process_all = false
    Boolean include_kinetics = false
    Boolean split_hd = false
    Int add_memory_mb = 0

    Int nproc
    Int mem_mb_per_core = 256
    String log_level = "INFO"
    String? chunk_str
  }
  Boolean all_kinetics = include_kinetics && process_all
  Boolean hifi_kinetics = include_kinetics && !(process_all)
  command {
    ccs \
      `readlink -f ${subreads}` \
      out.consensusreadset.xml \
      --log-level ${log_level} \
      ${"--chunk " + chunk_str} \
      ${true="" false="--noPolish" polishing} \
      ${true="--all" false="" process_all} \
      ${true="--all-kinetics --subread-fallback" false="" all_kinetics} \
      ${true="--hifi-kinetics" false="" hifi_kinetics} \
      --minLength ${min_length} \
      --maxLength ${max_length} \
      --minPasses ${min_passes} \
      --minSnr ${min_snr} \
      --minPredictedAccuracy ${min_predicted_accuracy} \
      ${draft_mode_args} \
      ${true="--disable-heuristics" false="" disable_heuristics} \
      ${true="--by-strand" false="" by_strand} \
      ${true="--hd-finder" false="" split_hd} \
      --unsorted-output \
      --alarms alarms.json \
      --task-report task-report.json \
      --report-json ccs_processing.report.json \
      --zmw-metrics-json ccs_zmws.json.gz \
      ${model_args} ${other_args} \
      -j ${nproc}
  }
  # this is absurdly inflated but we can't take risks right now
  Int total_mem_mb = add_memory_mb + 16384 + (nproc * mem_mb_per_core)
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File ccsxml = "out.consensusreadset.xml"
    File ccs_zmws = "ccs_zmws.json.gz"
    File report = "ccs_processing.report.json"
    File task_report = "task-report.json"
  }
}

task get_nchunks {
  input {
    Int bam_size
    Int max_nchunks
    Int target_size = 20000
    # This assumes subreads as input
    Int zmws_scale_factor = 10
  }
  command <<<
    python3 <<EOF
    import math
    bam_size = ~{bam_size}
    max_nchunks = ~{max_nchunks}
    target_size = ~{target_size}
    zmws_scale_factor = ~{zmws_scale_factor}
    nchunks = int(math.ceil(bam_size / (target_size * zmws_scale_factor)))
    nchunks = max(1, min(nchunks, max_nchunks))
    print("nchunks={}".format(nchunks))
    with open("nchunks.txt", "wt") as txt_out:
      txt_out.write(str(nchunks))
    EOF
  >>>
  runtime {
    cpu: 1
    memory: "100MB"
  }
  output {
    File txt_out = "nchunks.txt"
  }
}

task gather_ccs_zmws {
  input {
    Array[File] ccs_zmws
    Int add_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 -m pbcoretools.tasks.gather_ccs_zmws \
      --log-level ${log_level} \
      ccs_zmws.json.gz \
      ${sep=" " ccs_zmws}
  }
  # FIXME this task needs to be optimized
  Int total_mem_mb = add_memory_mb + 12 * 1024
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File gathered = "ccs_zmws.json.gz"
  }
}

task update_consensus_reads {
  input {
    File ccs_in
    File subreads_in
    Boolean use_run_design_uuid
    Int add_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
    Int memory_gb = 2
  }
  command {
    python3 -m pbcoretools.tasks.update_consensus_reads \
      --log-level ${log_level} \
      ${true="--use-run-design-uuid" false="" use_run_design_uuid} \
      ${ccs_in} \
      `readlink -f ${subreads_in}` \
      final.consensusreadset.xml
  }
  Int total_mem_mb = add_memory_mb + memory_gb * 1024
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File ccsxml = "final.consensusreadset.xml"
  }
}

task consolidate_reads_bam {
  input {
    File ccsxml
    File? zmws_json
    File? report_ccs_processing
    Int add_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 -m pbcoretools.tasks.consolidate_reads_bam \
      --log-level ${log_level} \
      ${"--zmws-json " + zmws_json} \
      ${"--report-ccs-processing " + report_ccs_processing} \
      ~{ccsxml}
  }
  Int total_mem_mb = add_memory_mb + 4096
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File reads_bam_fofn = "reads.fofn"
    File new_xml = "final.consensusreadset.xml"
  }
}

task pbreports_ccs2 {
  input {
    File ccsxml
    Int add_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }

  command {
    python3 \
      -m pbreports.report.ccs2 \
      --log-level ${log_level} \
      `readlink -f ${ccsxml}` \
      ccs.report.json \
      -c ccs.report.csv.zip
  }
  Int total_mem_mb = add_memory_mb + 12 * 1024
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File report = "ccs.report.json"
    File csv = "ccs.report.csv.zip"
    Array[File?] plot_pngs = glob("*.png")
  }
}

task cleanup_chunked_dataset_files {
  input {
    File chunked_xml
    # These files aren't actually used here - they're just a way to connect
    # tasks so that this is called last
    Array[File?] aux_files
    Int add_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 -m pbcoretools.tasks.delete_bam_resources \
      --log-level ${log_level} \
      ${chunked_xml}
  }
  Int total_mem_mb = add_memory_mb + 1024
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File removed_files_fofn = "removed_files.fofn"
  }
}

workflow pb_ccs {
  input {
    File eid_subread
    Float ccs_min_predicted_accuracy = 0.99
    Int ccs_min_passes = 3
    Float ccs_min_snr = 2.5
    Int ccs_min_length = 10
    Int ccs_max_length = 50000
    Boolean ccs_polish = true
    Boolean ccs_use_run_design_uuid = false
    String? ccs_model_args
    String? ccs_override_args
    Boolean ccs_by_strand = false
    Boolean ccs_process_all = false
    Boolean ccs_include_kinetics = false
    Boolean ccs_split_hd = false
    Boolean detect_methyl = false
    Int add_memory_mb = 0

    Int nproc = 1
    Int max_nchunks = 300
    Int target_size = 5000
    String log_level = "INFO"
    String? tmp_dir
  }

  Boolean flag_use_kinetics = ccs_include_kinetics || detect_methyl

  call memory.get_input_sizes {
    input:
      bam_dataset = eid_subread,
      log_level = log_level
  }

  # the built-in chunking in the CCS tool only works for a single BAM file
  # with no dataset filters.  we would prefer to use this if possible since
  # it's more efficient, but fall back on pbcore chunking otherwise
  # FIXME and we should just kill off the filters too
  Boolean is_single_bam_unfiltered = ((get_input_sizes.bam_files == 1) &&
                                      (get_input_sizes.bam_filters == 0))

  if (is_single_bam_unfiltered) {
    call get_nchunks {
      input:
        bam_size = get_input_sizes.bam_size,
        max_nchunks = max_nchunks,
        target_size = target_size
    }
    Int nchunks = read_int(get_nchunks.txt_out)

    scatter (i_chunk in range(nchunks)) {
      String chunk_str = "${i_chunk+1}/${nchunks}"
      call ccs as ccs1 {
        input:
          subreads = eid_subread,
          #subreads = subreads,
          min_passes = ccs_min_passes,
          min_predicted_accuracy = ccs_min_predicted_accuracy,
          min_snr = ccs_min_snr,
          min_length = ccs_min_length,
          max_length = ccs_max_length,
          model_args = ccs_model_args,
          polishing = ccs_polish,
          by_strand = ccs_by_strand,
          process_all = ccs_process_all,
          include_kinetics = flag_use_kinetics,
          split_hd = ccs_split_hd,
          other_args = ccs_override_args,
          add_memory_mb = add_memory_mb,
          nproc = nproc,
          chunk_str = chunk_str,
          log_level = log_level
      }
    }
  }
  # else...
  # FIXME maybe we should scatter over BAM files instead
  if (!is_single_bam_unfiltered) {
    call chunking.split_dataset as split_subreads {
      input:
        ds_in = eid_subread,
        max_nchunks = max_nchunks,
        target_size = target_size,
        split_args = "--keepReadGroups",
        memory_gb = get_input_sizes.index_memory_gb,
        base_memory_mb = add_memory_mb,
        log_level = log_level
    }

    scatter (subreads in split_subreads.chunks) {
      call ccs as ccs2 {
        input:
          subreads = subreads,
          min_passes = ccs_min_passes,
          min_predicted_accuracy = ccs_min_predicted_accuracy,
          min_snr = ccs_min_snr,
          min_length = ccs_min_length,
          max_length = ccs_max_length,
          model_args = ccs_model_args,
          polishing = ccs_polish,
          by_strand = ccs_by_strand,
          process_all = ccs_process_all,
          include_kinetics = flag_use_kinetics,
          split_hd = ccs_split_hd,
          other_args = ccs_override_args,
          add_memory_mb = add_memory_mb,
          nproc = nproc,
          log_level = log_level
      }
    }
  }

  Array[File] chunk_ccs_xmls = select_first([ccs1.ccsxml, ccs2.ccsxml])
  Array[File] chunk_ccs_zmws = select_first([ccs1.ccs_zmws, ccs2.ccs_zmws])
  Array[File] chunk_ccs_reports = select_first([ccs1.report, ccs2.report])
  Array[File] chunk_ccs_task_reports = select_first([ccs1.task_report, ccs2.task_report])

  call chunking.gather_datasets_lite as gather_ccsxml {
    input:
      chunks = chunk_ccs_xmls,
      dataset_name = "'CCS Data'",
      dataset_ext = ".consensusreadset.xml",
      prefix = "ccs_reads"
  }

  call chunking.gather_generic as gather_processing_reports {
    input:
      chunks = chunk_ccs_reports,
      output_file_name = "ccs_processing.report.json",
      gather_args = "--dataset ${update_consensus_reads.ccsxml}",
      log_level = log_level
  }

  call gather_ccs_zmws {
    input:
      ccs_zmws = chunk_ccs_zmws,
      add_memory_mb = add_memory_mb,
      log_level = log_level
  }

  call pbreports.gather_task_reports {
    input:
      task_reports = chunk_ccs_task_reports,
      output_file_name = "ccs_tasks_report.json",
      log_level = log_level
  }

  call update_consensus_reads {
    input:
      ccs_in = gather_ccsxml.gathered,
      subreads_in = eid_subread,
      use_run_design_uuid = ccs_use_run_design_uuid,
      add_memory_mb = add_memory_mb,
      nproc = nproc,
      log_level = log_level
  }

  call pbcoretools.auto_ccs_outputs as export_bam {
    input:
      ccsxml = update_consensus_reads.ccsxml,
      mode = "consolidate",
      nproc = nproc,
      add_memory_mb = add_memory_mb,
      log_level = log_level
  }

  call pbcoretools.auto_ccs_outputs as export_fasta {
    input:
      ccsxml = update_consensus_reads.ccsxml,
      mode = "fasta",
      nproc = nproc,
      add_memory_mb = add_memory_mb,
      log_level = log_level
  }

  call pbcoretools.auto_ccs_outputs as export_fastq {
    input:
      ccsxml = update_consensus_reads.ccsxml,
      mode = "fastq",
      nproc = nproc,
      add_memory_mb = add_memory_mb,
      log_level = log_level
  }

  call consolidate_reads_bam {
    input:
      ccsxml = update_consensus_reads.ccsxml,
      zmws_json = gather_ccs_zmws.gathered,
      report_ccs_processing = gather_processing_reports.gathered,
      add_memory_mb = add_memory_mb,
  }
  # this avoids reading in the BAM for call caching
  String reads_bam_file = read_string(consolidate_reads_bam.reads_bam_fofn)

  # making the final cleanup task dependent on these files forces it to be
  # run last
  Array[File?] aux_files = [pbreports_ccs2.report,
                            export_bam.datastore,
                            export_fasta.datastore,
                            export_fastq.datastore,
                            consolidate_reads_bam.reads_bam_fofn]

  call cleanup_chunked_dataset_files {
    input:
      chunked_xml = update_consensus_reads.ccsxml,
      aux_files = aux_files,
  }

  if (detect_methyl) {
    call pb_detect_methyl.pbprimrose {
      input:
        reads = consolidate_reads_bam.new_xml,
        keep_kinetics = ccs_include_kinetics,
        reuse_source_uuid = true,
        nproc = nproc,
        base_memory_mb = add_memory_mb,
        log_level = log_level
    }
    String primrose_bam = read_string(pbprimrose.ccsbam_fofn)

    call pb_detect_methyl.pbreports_detect_cpg_methyl {
      input:
        ccsxml = pbprimrose.ccsxml,
        nproc = nproc,
        base_memory_mb = add_memory_mb,
        log_level = log_level
    }
  }

  File final_xml = select_first([pbprimrose.ccsxml,
                                 consolidate_reads_bam.new_xml])
  String final_reads_bam_file = select_first([primrose_bam, reads_bam_file])

  call pbreports_ccs2 {
    input:
      ccsxml = final_xml,
      add_memory_mb = add_memory_mb,
      nproc = nproc,
      log_level = log_level
  }

  output {
    File ccsxml = final_xml
    File report_ccs2 = pbreports_ccs2.report
    File report_ccs_processing = gather_processing_reports.gathered
    File ccs_zmws = gather_ccs_zmws.gathered
    File ccs_csv = pbreports_ccs2.csv
    File bam_datastore = export_bam.datastore
    File fasta_datastore = export_fasta.datastore
    File fastq_datastore = export_fastq.datastore
    String reads_bam = final_reads_bam_file
    File report_tasks = gather_task_reports.report
    File? report_detect_cpg_methyl = pbreports_detect_cpg_methyl.report
  }
}
