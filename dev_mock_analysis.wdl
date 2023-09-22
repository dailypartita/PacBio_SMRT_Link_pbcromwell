# Mock analysis workflow for testing auto-analysis behavior

version 1.0

task mock_analysis {
  input {
    File eid_ccs
    File eid_ref_dataset
    Boolean param_a = false
    String param_b = ""
    Int param_c = 0
    Float param_d = 0.99

    String log_level = "INFO"
    Int nproc = 1
  }
  command <<<
    python3 <<EOF
    import os.path
    from pbcommand.models.report import Report, Attribute
    from pbcore.io import ConsensusReadSet, ReferenceSet
    ds = ConsensusReadSet(os.path.realpath("~{eid_ccs}"), skipCounts=True)
    ref = ReferenceSet(os.path.realpath("~{eid_ref_dataset}"), skipCounts=True)
    Report(
      "mock_analysis",
      attributes=[
        Attribute("metric_a", value="~{param_a}" == "true"),
        Attribute("metric_b", value="~{param_b}"),
        Attribute("metric_c", value=~{param_c}),
        Attribute("metric_d", value=~{param_d}),
      ],
      title="Mock Analysis Report",
      dataset_uuids=(ds.uuid,ref.uuid)).write_json("mock_analysis.report.json")
    EOF
  >>>
  runtime {
    cpu: 1
    memory: "500MB"
    backend: "Local"
  }
  output {
    File report = "mock_analysis.report.json"
  }
}

workflow dev_mock_analysis {
  input {
    File eid_ccs
    File eid_ref_dataset
    Boolean param_a = false
    String param_b = ""
    Int param_c = 0
    Float param_d = 0.99

    Int nproc = 1
    Int max_nchunks = 100
    String log_level = "INFO"
    String? tmp_dir
  }

  call mock_analysis {
    input:
      eid_ccs = eid_ccs,
      eid_ref_dataset = eid_ref_dataset,
      param_a = param_a,
      param_b = param_b,
      param_c = param_c,
      param_d = param_d,
      log_level = log_level
  }

  output {
    File report_mock_analysis = mock_analysis.report
  }
}
