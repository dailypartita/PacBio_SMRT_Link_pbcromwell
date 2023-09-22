# Undo Demultiplexing

version 1.0

import "pb_ccs.wdl"

task lima_undo {
  input {
    File ccsxml
    String prefix

    String log_level = "INFO"
    Int nproc = 1
  }
  String ofn = "${prefix}.consensusreadset.xml"
  command {
    set -e

    dataset \
      --log-level ${log_level} \
      --log-file dataset-absolutize.log \
      absolutize \
      --output-xml absolutized.consensusreadset.xml \
      `readlink -f ${ccsxml}`

    lima-undo \
      --log-level ${log_level} \
      --log-file lima-undo.log \
      --num-threads ${nproc} \
      --alarms alarms.json \
      absolutized.consensusreadset.xml \
      ${ofn}

    dataset newuuid --random ${ofn}
  }
  runtime {
    cpu: nproc
    memory: "4GB"
  }
  output {
    File undemuxed = "${ofn}"
    File? alarms = "alarms.json"
  }
}

workflow pb_undo_demux {
  input {
    File eid_ccs
    String output_file_prefix = "all_samples_with_barcodes"

    String log_level = "INFO"
    Int nproc = 1
    # not used
    Int max_nchunks = 1
    String? tmp_dir
  }

  call lima_undo {
    input:
      ccsxml = eid_ccs,
      prefix = output_file_prefix,
      log_level = log_level,
      nproc = nproc
  }

  call pb_ccs.pbreports_ccs2 {
    input:
      ccsxml = lima_undo.undemuxed,
      log_level = log_level
  }

  output {
    File undemuxed = lima_undo.undemuxed
    File report_ccs2 = pbreports_ccs2.report
  }
}
