# Mock version of CCS Analysis - requires PacBioTestData

version 1.0

task mock_ccs {
  input {
    File eid_subread
    Boolean use_run_design_uuid = false
    Boolean include_kinetics = false
    Boolean process_all = false
    Boolean detect_methyl = false
    Boolean split_hd = false

    Int nproc = 1
    String log_level = "INFO"
  }
  command <<<
    python3 <<EOF
    import json
    import os
    from pbcore.io import SubreadSet, ConsensusReadSet
    import pbcore.data
    from pbcommand.models.report import Report, Attribute
    ccs_bam = pbcore.data.getCCSBAM()
    ds_ccs = ConsensusReadSet(ccs_bam)
    ds_clr = SubreadSet(os.path.realpath("~{eid_subread}"))
    ds_ccs.metadata.collections = ds_clr.metadata.collections
    # FIXME this is an obnoxious bifurcation in our data model due to lima
    for bs in ds_clr.metadata.bioSamples:
      ds_ccs.metadata.bioSamples.addSample(bs.name)
      for bc in bs.DNABarcodes:
        ds_ccs.metadata.bioSamples[-1].DNABarcodes.addBarcode(bc.name)
    ds_ccs.name = ds_clr.name + " (CCS)"
    if "~{use_run_design_uuid}" == "true":
      ds_ccs.uuid = ds_ccs.metadata.collections[0].consensusReadSetRef.uuid
    else:
      ds_ccs.newUuid()
    ds_ccs.write("mock.consensusreadset.xml")
    Report(
      "ccs2",
      attributes=[
        Attribute("number_of_ccs_reads", value=2),
        Attribute("total_number_of_ccs_bases", value=2000),
        Attribute("include_kinetics", value="~{include_kinetics}" == True),
        Attribute("process_all", value="~{process_all}" == True),
        Attribute("detect_methyl", value="~{detect_methyl}" == True),
        Attribute("split_hd", value="~{split_hd}" == True)
      ],
      dataset_uuids=(ds_ccs.uuid,)).write_json("mock_ccs.report.json")
    EOF
  >>>
  runtime {
    memory: "1GB"
    cpu: 1
  }
  output {
    File xml = "mock.consensusreadset.xml"
    File report = "mock_ccs.report.json"
  }
}

workflow dev_mock_ccs {
  input {
    File eid_subread
    Boolean ccs_use_run_design_uuid = false
    Boolean ccs_include_kinetics = false
    Boolean ccs_process_all = false
    Boolean detect_methyl = false
    Boolean ccs_split_hd = false

    Int nproc = 1
    Int max_nchunks = 100
    Int target_size = 20000
    String log_level = "INFO"
    String? tmp_dir
  }

  call mock_ccs {
    input:
      eid_subread = eid_subread,
      use_run_design_uuid = ccs_use_run_design_uuid,
      include_kinetics = ccs_include_kinetics,
      process_all = ccs_process_all,
      detect_methyl = detect_methyl,
      split_hd = ccs_split_hd,
      log_level = log_level,
      nproc = 1
  }

  output {
    File ccsxml = mock_ccs.xml
    File report_ccs2 = mock_ccs.report
  }
}
