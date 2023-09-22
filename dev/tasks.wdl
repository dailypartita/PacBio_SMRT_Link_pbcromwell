version 1.0

task dev_failing_task {
  input {
    Int? exit_code
    Int sleep_time = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    echo "Sleeping for ${sleep_time}s..."
    sleep ${sleep_time}
    echo "Exiting with exit code ${exit_code}" >/dev/stderr
    echo "Exiting with exit code ${exit_code}" > status.txt
    exit ${exit_code}
  }
  runtime {
    backend: "Local"
    cpu: 1
    memory: "1GB"
  }
  output {
    File status = "status.txt"
  }
}

task samples_report {
  input {
    File dataset

    Int nproc = 1
    String log_level = "INFO"
  }
  command <<<
    python3 <<EOF
    from pbcore.io import openDataSet
    from pbcommand.models.report import Report, Attribute
    import os.path
    def _to_samples_str(s): return ";".join(sorted(list(s)))
    ds_file = os.path.realpath("~{dataset}")
    ds = openDataSet(ds_file, skipCounts=True)
    well_samples = {c.wellSample.name for c in ds.metadata.collections}
    bio_samples = {bs.name for bs in ds.metadata.bioSamples}
    attributes = [
      Attribute("dataset_name", value=ds.name),
      Attribute("well_samples", value=_to_samples_str(well_samples)),
      Attribute("bio_samples", value=_to_samples_str(bio_samples))
    ]
    r = Report("samples_report", title="Dataset Samples",
               attributes=attributes)
    r.write_json("samples_info.report.json")
    EOF
  >>>
  runtime {
    backend: "Local"
    cpu: 1
    memory: "1GB"
  }
  output {
    File report = "samples_info.report.json"
  }
}

task get_chem_bundle_version {
  input {
    File dataset

    Int nproc = 1
    String log_level = "INFO"
  }
  command <<<
    python3 <<EOF
    from pbcommand.models.report import Report, Attribute
    import xml.dom.minidom
    import os
    CHEM_DIR = os.environ["SMRT_CHEMISTRY_BUNDLE_DIR"]
    manifest = os.path.join(CHEM_DIR, "manifest.xml")
    dom = xml.dom.minidom.parse(manifest)
    version = dom.getElementsByTagName("Version")[0].lastChild.data
    attributes = [
      Attribute("bundle_version", value=version),
      Attribute("bundle_dir", value=version)
    ]
    r = Report("bundle_report", title="Bundle Report", attributes=attributes)
    r.write_json("bundle_info.report.json")
    EOF
  >>>
  runtime {
    backend: "Local"
    cpu: 1
    memory: "1GB"
  }
  output {
    File report = "bundle_info.report.json"
  }
}

task populate_datastore {
  input {
    File subreads
    Int num_subreadsets = 25
    Int sleep_multiplier = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 \
      -m pbcromwell.testkit.dev_tasks.populate_datastore \
      --log-level ${log_level} \
      --num-subreadsets ${num_subreadsets} \
      --sleep-multiplier ${sleep_multiplier} \
      `readlink -f ${subreads}` \
      mock.datastore.json
  }
  runtime {
    cpu: 1
    cpu: 1
    memory: "1GB"
  }
  output {
    File datastore = "mock.datastore.json"
  }
}

task subreads_to_subreads {
  input {
    File subreads
    String use_uuid = ""

    Int nproc = 1
    String log_level = "INFO"
  }
  command <<<
    python3 <<EOF
    from pbcore.io import SubreadSet
    from pbcommand.models.report import Report, Attribute
    import os.path
    import time
    import uuid
    subreads_file = os.path.realpath("~{subreads}")
    use_uuid = "~{use_uuid}"
    with SubreadSet(subreads_file, skipCounts=False) as ds:
        ds.name = ds.name + " (demultiplexed) (MOCK)"
        ds.tags = "internal,testdata"
        if use_uuid:
            uuid_ = uuid.UUID(use_uuid)
            ds.uuid = str(uuid_)
        else:
            ds.newUuid(random=True)
        ds.write("mock.subreadset.xml")
        rpt = Report("dev_subread_report",
                     attributes=[Attribute("n_reads", value=len(ds))],
                     dataset_uuids=[ds.uuid])
        rpt.write_json("subreads.report.json")
        rpt2 = Report("dev_misc_report",
                      attributes=[Attribute("time", value=time.time())])
        rpt2.write_json("dev_misc.report.json")
    EOF
  >>>
  runtime {
    backend: "Local"
    cpu: 1
    memory: "1GB"
  }
  output {
    File subreads_out = "mock.subreadset.xml"
    File report_subreads = "subreads.report.json"
    File report_misc = "dev_misc.report.json"
  }
}

task subreads_to_ccs {
  input {
    File subreads
    String use_uuid = ""

    Int nproc = 1
    String log_level = "INFO"
  }
  command <<<
    python3 <<EOF
    from pbcore.io import SubreadSet, ConsensusReadSet
    import pbcore.data
    from pbcommand.models.report import Report, Attribute
    import os.path
    import time
    import uuid
    subreads_file = os.path.realpath("~{subreads}")
    use_uuid = "~{use_uuid}"
    ccs_file_name = pbcore.data.getCCSBAM()
    with SubreadSet(subreads_file, skipCounts=True) as ds:
        with ConsensusReadSet(ccs_file_name) as ccs:
            ccs.name = ds.name + " (CCS) (MOCK)"
            ccs.tags = "internal,testdata,pbcore"
            if use_uuid:
                uuid_ = uuid.UUID(use_uuid)
                ccs.uuid = str(uuid_)
            else:
                ds.newUuid(random=True)
            ccs.write("mock.consensusreadset.xml")
            rpt = Report("dev_ccs_report",
                         attributes=[Attribute("n_reads", value=len(ccs))],
                         dataset_uuids=[ccs.uuid])
            rpt.write_json("ccs.report.json")
            rpt2 = Report("dev_misc_report",
                          attributes=[Attribute("time", value=time.time())])
            rpt2.write_json("dev_misc.report.json")
    EOF
  >>>
  runtime {
    backend: "Local"
    cpu: 1
    memory: "1GB"
  }
  output {
    File ccs_out = "mock.consensusreadset.xml"
    File report_ccs = "ccs.report.json"
    File report_misc = "dev_misc.report.json"
  }
}

task ccs_to_ccs {
  input {
    File ccsxml
    String use_uuid = ""

    Int nproc = 1
    String log_level = "INFO"
  }
  command <<<
    python3 <<EOF
    from pbcore.io import ConsensusReadSet
    from pbcommand.models.report import Report, Attribute
    import os.path
    import time
    import uuid
    ccs_file = os.path.realpath("~{ccsxml}")
    use_uuid = "~{use_uuid}"
    with ConsensusReadSet(ccs_file) as ds:
        if use_uuid:
            uuid_ = uuid.UUID(use_uuid)
            ds.uuid = str(uuid_)
        else:
            ds.newUuid(random=True)
        ds.name += " (copy) (MOCK)"
        ds.tags = "internal,testdata"
        ds.write("mock.consensusreadset.xml")
        rpt = Report("dev_ccs_report",
                     attributes=[Attribute("n_reads", value=len(ds))],
                     dataset_uuids=[ds.uuid])
        rpt.write_json("ccs.report.json")
        rpt2 = Report("dev_misc_report",
                      attributes=[Attribute("time", value=time.time())])
        rpt2.write_json("dev_misc.report.json")
    EOF
  >>>
  runtime {
    backend: "Local"
    cpu: 1
    memory: "1GB"
  }
  output {
    File ccs_out = "mock.consensusreadset.xml"
    File report_ccs = "ccs.report.json"
    File report_misc = "dev_misc.report.json"
  }
}

task dump_environment {
  input {
    Int sleep_time = 0
    Int nproc = 1
    String log_level = "INFO"
  }
  command <<<
    sleep ~{sleep_time}
    env | sort > environment.txt
    python3 <<EOF
    import os, json
    with open("environment.json", "w") as env_out:
      env_out.write(json.dumps(dict(os.environ),
                               indent=2,
                               separators=(',', ':'),
                               sort_keys=True))
    EOF
  >>>
  runtime {
    cpu: 1
    memory: "100MB"
  }
  output {
    File env_out = "environment.txt"
    File env_json = "environment.json"
  }
}

# this does nothing except generate an alarms.json file and exit 0
task emit_alarm {
  input {
    File dataset

    Int nproc = 1
    String log_level = "INFO"
  }
  command <<<
    python3 <<EOF
    from pbcommand.models.common import PacBioAlarm
    import logging
    alarm = PacBioAlarm(
      exception=None,
      info=None,
      message="This is a job alarm!",
      name="Test Alarm",
      severity=logging.WARN,
      owner="python")
    alarm.to_json("alarms.json")
    open("test_out.txt", "w").write("Hello, world!")
    EOF
  >>>
  runtime {
    backend: "Local"
    cpu: 1
    memory: "1GB"
  }
  output {
    File test_out = "test_out.txt"
  }
}

task raise_alarm {
  input {
    File dataset

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 -m pbcromwell.testkit.dev_tasks.raise_alarm \
      --log-level ${log_level} \
      `readlink -f ${dataset}` \
      out.txt
  }
  runtime {
    backend: "Local"
    cpu: 1
    memory: "1GB"
  }
  output {
    File test_out = "out.txt"
  }
}

task mock_demux {
  input {
    File dataset
    String dataset_ext = ".subreadset.xml"

    Int nproc
    String log_level = "DEBUG"
  }
  command <<<
    python3 <<EOF
    from pbcore.io import openDataSet
    from pbcommand.models.common import DataStore, DataStoreFile
    import os.path as op
    import uuid
    dataset_ext = "~{dataset_ext}"
    ds = openDataSet(op.realpath("~{dataset}"))
    ds_name = ds.name
    files = []
    for bioSample in ds.metadata.bioSamples:
      for dna_bc in bioSample.DNABarcodes:
        ds.name = ds_name + " ({b})".format(b=dna_bc.name)
        ds.uuid = str(uuid.uuid4())
        ds_out = "lima.{b}{e}".format(b=dna_bc.name, e=dataset_ext)
        ds.write(ds_out)
        files.append(DataStoreFile(ds.uuid,
                                   "barcoding.tasks.lima-0",
                                   ds.datasetType,
                                   op.abspath(ds_out)))
    datastore = DataStore(files)
    datastore.write_json("mock_demux.datastore.json")
    EOF
  >>>
  runtime {
    cpu: 1
    cpu: 1
    memory: "4GB"
  }
  output {
    File datastore = "mock_demux.datastore.json"
  }
}

task mock_update_barcoded_sample_metadata {
  input {
    File lima_datastore
    File input_reads
    File barcodes
    Boolean use_barcode_uuids = false
    String extension = ".subreadset.xml"

    Int nproc
    String log_level = "DEBUG"
  }
  command <<<
    python3 <<EOF
    from pbcoretools.file_utils import mock_update_barcoded_sample_metadata
    import os.path as op
    use_barcode_uuids = "~{use_barcode_uuids}" == "true"
    datastore = op.realpath("~{lima_datastore}")
    datastore = mock_update_barcoded_sample_metadata(
        base_dir=op.dirname(datastore),
        datastore_file=datastore,
        input_reads=op.realpath("~{input_reads}"),
        barcode_set=op.realpath("~{barcodes}"),
        use_barcode_uuids=use_barcode_uuids,
        extension="~{extension}")
    datastore.write_json("mock_updated.datastore.json")
    EOF
  >>>
  runtime {
    cpu: 1
    memory: "1GB"
  }
  output {
    File datastore = "mock_updated.datastore.json"
  }
}
