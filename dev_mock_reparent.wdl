# dev workflow for testing behavior of dataset parentage
version 1.0

task make_new_parent {
  input {
    File subreads
  }
  command {
    dataset absolutize \
      `readlink -f ${subreads}` \
      --output-xml new_parent.subreadset.xml
    dataset newuuid new_parent.subreadset.xml --random
  }
  runtime {
    backend: "Local"
    cpu: 1
  }
  output {
    File parent = "new_parent.subreadset.xml"
  }
}

task make_child {
  input {
    File parent
  }
  command <<<
    python3 <<EOF
    from pbcore.io import SubreadSet
    ds = SubreadSet("~{parent}")
    ds.metadata.addParentDataSet(ds.uuid,
                                 ds.datasetType,
                                 createdBy="AnalysisJob",
                                 timeStampedName="")
    ds.newUuid(random=True)
    ds.name = ds.name + " (child)"
    ds.write("child.subreadset.xml")
    EOF
  >>>
  runtime {
    backend: "Local"
    cpu: 1
  }
  output {
    File child = "child.subreadset.xml"
  }
}

workflow dev_mock_reparent {
  input {
    File eid_subread

    Int nproc = 1
    Int max_nchunks = 1
    String log_level = "INFO"
    String? tmp_dir
  }

  call make_new_parent {
    input:
      subreads = eid_subread
  }

  call make_child as make_child_1 {
    input:
      parent = eid_subread
  }

  call make_child as make_child_2 {
    input:
      parent = make_new_parent.parent
  }

  output {
    # this file will have a parent UUID known to smrtlink
    File child_1 = make_child_1.child
    # this file will have a randomly generated parent UUID
    File child_2 = make_child_2.child
  }
}
