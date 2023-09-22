# workflow for testing SL-4200.  This will succeed, but SMRT Link will fail
# when it tries to import the datastore

version 1.0

task emit_bad_datastore {
  input {
    File subreads

    Int nproc = 1
    String log_level = "INFO"
  }
  command <<<
    python3 <<EOF
    import uuid
    from pbcommand.models.common import DataStoreFile, DataStore, FileTypes
    files = [
      DataStoreFile(
        uuid.uuid4(),
        "emit-bad-datastore-0",
        FileTypes.DS_SUBREADS.file_type_id,
        "/bad/path/dataset.subreadset.xml")
    ]
    datastore = DataStore(files)
    datastore.write_json("mock.datastore.json")
    EOF
  >>>
  runtime {
    cpu: 1
  }
  output {
    File datastore = "mock.datastore.json"
  }
}

workflow dev_bad_datastore {
  input {
    File eid_subread

    Int nproc = 1
    Int max_nchunks = 100
    String log_level = "INFO"
    String? tmp_dir
  }

  call emit_bad_datastore {
    input:
      subreads = eid_subread,
      nproc = nproc,
      log_level = log_level
  }

  output {
    File mock_datastore = emit_bad_datastore.datastore
  }
}
