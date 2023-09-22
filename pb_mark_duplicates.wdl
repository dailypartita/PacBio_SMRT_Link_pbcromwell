version 1.0

import "wf_prepare_input.wdl"
import "tasks/pbreports.wdl"

task pbmarkdup {
  input {
    File reads_in
    Boolean cross_library = true

    Int nproc
    String log_level = "INFO"
  }
  command {
    pbmarkdup \
      --log-level ${log_level} \
      -j ${nproc} \
      ${true="--cross-library" false="" cross_library} \
      ${reads_in} \
      deduplicated.ccs.bam \
      --dup-file duplicates.ccs.bam
    pbindex duplicates.ccs.bam
    dataset \
      --strict \
      --log-level ${log_level} \
      create \
      --type ConsensusReadSet \
      --name "PCR duplicates" \
      duplicates.consensusreadset.xml \
      duplicates.ccs.bam
  }
  runtime {
    cpu: nproc
  }
  output {
    File deduplicated = "deduplicated.ccs.bam"
    File duplicates = "duplicates.ccs.bam"
    # this file is not actually exposed, but it is used to generate the report
    # (or rather the underlying .pbi index is)
    File duplicates_xml = "duplicates.consensusreadset.xml"
  }
}

# FIXME this should probably live in C++
task create_deduplicated_dataset {
  input {
    File reads_in
    File ccs_bam

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 <<EOF
    import os.path as op
    from pbcore.io import ConsensusReadSet
    ds_in = ConsensusReadSet("~{reads_in}", skipCounts=True)
    ds_out = ConsensusReadSet(op.realpath("~{ccs_bam}"))
    ds_out.name = ds_in.name + " (deduplicated)"
    ds_out.tags = ",".join(ds_in.tags.split(",") + ["deduplicated"])
    ds_out.write("deduplicated.consensusreadset.xml")
    EOF
  }
  runtime {
    cpu: 1
    memory: "4GB"
  }
  output {
    File ccsxml = "deduplicated.consensusreadset.xml"
  }
}

task pbreports_duplicates {
  input {
    File reads_out
    File duplicates

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 \
      -m pbreports.report.mark_duplicates \
      --log-level ${log_level} \
      ${reads_out} ${duplicates} \
      duplicates.report.json
  }
  runtime {
    cpu: 1
    memory: "4GB"
  }
  output {
    File report = "duplicates.report.json"
  }
}

workflow pb_mark_duplicates {
  input {
    File eid_ccs
    Boolean cross_library = true
    Int filter_min_qv = 20

    Int nproc
    String log_level = "INFO"
    Int max_nchunks = 0
    String? tmp_dir
  }

  call wf_prepare_input.prepare_input {
    input:
      dataset_xml = eid_ccs,
      filter_min_qv = filter_min_qv,
      nproc = nproc,
      log_level = log_level
  }

  call pbmarkdup {
    input:
      reads_in = prepare_input.reads_file,
      cross_library = cross_library,
      log_level = log_level,
      nproc = nproc
  }

  call create_deduplicated_dataset {
    input:
      reads_in = prepare_input.reads_file,
      ccs_bam = pbmarkdup.deduplicated,
      log_level = log_level,
      nproc = 1
  }

  call pbreports_duplicates {
    input:
      reads_out = create_deduplicated_dataset.ccsxml,
      duplicates = pbmarkdup.duplicates_xml,
      log_level = log_level,
      nproc = 1
  }

  output {
    File deduplicated = create_deduplicated_dataset.ccsxml
    File duplicates = pbmarkdup.duplicates
    File report_duplicates = pbreports_duplicates.report
  }
}
