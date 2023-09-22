# Standalone methylation detection

version 1.0

import "wf_prepare_input.wdl"

task pbprimrose {
  input {
    File reads
    Boolean keep_kinetics = false
    Boolean reuse_source_uuid = false
    # XXX this would be a string ' (5mC)', but there is a giant bug in the
    # Python wdl_parser module that we use to generate the pipeline JSONs for
    # SMRT Link.  since wdl_parser is unmaintained we probably need to switch
    # to miniWDL but this will require build modifications.  so instead we
    # construct the argument in Bash in the command block
    Boolean? append_dataset_name
    Int base_memory_mb = 0

    Int nproc
    String log_level = "INFO"
  }
  command <<<
    set -e
    append_dataset_args=()
    if [ "~{append_dataset_name}" = "true" ]; then
      append_dataset_args=(--append-dataset-name " (5mC)")
    fi
    primrose \
      --log-level DEBUG \
      --log-file primrose.log \
      --alarms alarms.json \
      --num-threads ~{nproc} \
      "${append_dataset_args[@]}" \
      ~{true="--reuse-source-uuid" false="" reuse_source_uuid} \
      ~{true="--keep-kinetics" false="" keep_kinetics} \
      --qv-histogram-report primrose.report.json \
      `readlink -f ~{reads}` \
      with_5mC.consensusreadset.xml
    echo `pwd`/with_5mC.bam > ccsbam.fofn
  >>>
  Int total_mem_mb = base_memory_mb + 8192
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File ccsxml = "with_5mC.consensusreadset.xml"
    File ccsbam_fofn = "ccsbam.fofn"
    File report = "primrose.report.json"
    File? alarms = "alarms.json"
  }
}

task pbreports_detect_cpg_methyl {
  input {
    File ccsxml
    Int base_memory_mb = 0

    Int nproc
    String log_level = "INFO"
  }
  command {
    python3 -m pbreports.report.detect_cpg_methyl \
      --log-level DEBUG \
      --log-file pbreports.log \
      `readlink -f ${ccsxml}` \
      detect_cpg_methyl.report.json
  }
  Int total_mem_mb = 512 + base_memory_mb
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File report = "detect_cpg_methyl.report.json"
    Array[File?] plot_pngs = glob("*.png")
  }
}

workflow pb_detect_methyl {
  input {
    File eid_ccs
    Boolean keep_kinetics = false
    Int filter_min_qv = 20
    String dataset_filters = ""
    Int add_memory_mb = 0

    Int nproc = 1
    Int max_nchunks = 300
    Int target_size = 5000
    String log_level = "INFO"
    String? tmp_dir
  }

  call wf_prepare_input.prepare_input {
    input:
      dataset_xml = eid_ccs,
      dataset_filters = dataset_filters,
      filter_min_qv = filter_min_qv,
      nproc = 1,
      log_level = log_level
  }

  call pbprimrose {
    input:
      reads = prepare_input.reads_file,
      keep_kinetics = keep_kinetics,
      append_dataset_name = true,
      base_memory_mb = add_memory_mb,
      nproc = nproc,
      log_level = log_level
  }

  call pbreports_detect_cpg_methyl {
    input:
      ccsxml = pbprimrose.ccsxml,
      base_memory_mb = add_memory_mb,
      nproc = nproc,
      log_level = log_level
  }

  output {
    File ccsxml = pbprimrose.ccsxml
    String ccsbam = read_string(pbprimrose.ccsbam_fofn)
    File report_detect_cpg_methyl = pbreports_detect_cpg_methyl.report
    File? alarms = pbprimrose.alarms
  }
}
