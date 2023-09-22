version 1.0

task fasta_to_reference {
  input {
    File reference_fasta
    String reference_name
    String organism = "unknown"
    String ploidy = "haploid"

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    cp -H ${reference_fasta} reference.fasta
    samtools faidx reference.fasta
    dataset \
      --strict \
      --log-level ${log_level} \
      create \
      --type ReferenceSet \
      --name "${reference_name}" \
      --organism "${organism}" \
      --ploidy "${ploidy}" \
      smrtlink_reference.referenceset.xml \
      reference.fasta
  }
  runtime {
    cpu: 1
    memory: "4GB"
  }
  output {
    File referenceset = "smrtlink_reference.referenceset.xml"
  }
}

task fasta_to_barcodes {
  input {
    File barcodes_fasta
    String dataset_name

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    cp -H ${barcodes_fasta} barcodes.fasta
    samtools faidx barcodes.fasta
    dataset \
      --strict \
      --log-level ${log_level} \
      create \
      --type BarcodeSet \
      --name "${dataset_name}" \
      smrtlink_barcodes.barcodeset.xml \
      barcodes.fasta
  }
  runtime {
    cpu: 1
    backend: "Local"
    memory: "1GB"
  }
  output {
    File barcodeset = "smrtlink_barcodes.barcodeset.xml"
  }
}

# This task combines multiple validation functions used in SMRT Link
# integration, but it leaves testing of the results to the runner code.
task simple_dataset_report {
  input {
    File? dataset
    Boolean skip_counts = false
    File? skip_if_report_exists

    Int nproc = 1
    String log_level = "INFO"
  }
  command <<<
    python3 <<EOF
    import sys
    skip_if_report_exists = "~{skip_if_report_exists}"
    if skip_if_report_exists:
      sys.exit(0)
    from pbcore.io import openDataSet
    from pbcommand.models.report import Report, Attribute
    ds_file = "~{dataset}"
    ds = openDataSet(ds_file, skipCounts=~{true="True" false="False" skip_counts})
    ds.updateCounts()
    attributes = [
      Attribute("dataset_name", ds.name, name="Dataset Name"),
      Attribute("number_of_records", ds.numRecords, name="Number of Records"),
      Attribute("total_length", ds.totalLength, name="Total Length")
    ]
    r = Report("simple_dataset_report",
               title="Dataset Details",
               attributes=attributes,
               dataset_uuids=[ds.uuid])
    r.write_json("dataset_info.report.json")
    EOF
  >>>
  runtime {
    backend: "Local"
    cpu: 1
    memory: "1GB"
  }
  output {
    File? report = "dataset_info.report.json"
  }
}
