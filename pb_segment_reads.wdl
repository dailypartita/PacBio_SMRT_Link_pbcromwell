# Read segmentation workflow for single-cell-Iso-Seq (and other future apps)

version 1.0

task prepare_inputs {
  input {
    File ccs_xml
    File? eid_barcode
    File? adapters_fasta

    Int nproc = 1
    String log_level = "INFO"
  }
  String prefix = "segmented"
  command <<<
    set -e
    python3 <<EOF
    import os.path
    import sys
    from pbcore.io import BarcodeSet, ConsensusReadSet
    from pbcommand.models.common import PacBioAlarm
    if len("~{adapters_fasta}") > 0:
      os.symlink(os.path.realpath("~{adapters_fasta}"), "adapters.fasta")
    elif len("~{eid_barcode}") > 0:
      bcs = BarcodeSet(os.path.realpath("~{eid_barcode}"))
      os.symlink(bcs.externalResources[0].resourceId, "adapters.fasta")
    else:
      PacBioAlarm.dump_error(
        "alarms.json",
        "MissingAdaptersError",
        "Either a BarcodeSet XML or an adapters FASTA is required"
        "Either a BarcodeSet XML or an adapters FASTA is required"
        "MissingAdaptersError",
        "ERROR")
      sys.exit(1)
    ccs = ConsensusReadSet(os.path.realpath("~{ccs_xml}"), trustCounts=True)
    for i, resource in enumerate(ccs.externalResources, start=1):
      # we can get away with skipCounts=True here because skera doesn't care
      chunk = ConsensusReadSet(resource.resourceId, skipCounts=True)
      chunk.write(f"chunk.{i}.consensusreadset.xml")
    EOF
  >>>
  runtime {
    cpu: 1
    # not loading CCS BAM indexes keeps this low
    memory: "1GB"
  }
  output {
    File adapters = "adapters.fasta"
    Array[File] chunks = glob("chunk.*.consensusreadset.xml")
  }
}

task skera {
  input {
    File ccs_xml
    File adapters_fasta
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  String prefix = "segmented"
  command <<<
    skera split \
      --log-level ~{log_level} \
      --log-file skera.log \
      --alarms alarms.json \
      -j ~{nproc} \
      $(readlink -f "~{ccs_xml}") \
      "~{adapters_fasta}" \
      "~{prefix}.consensusreadset.xml"
    echo "`pwd`/~{prefix}.bam" > bam.fofn
    echo "`pwd`/~{prefix}.non_passing.bam" > non_passing_bam.fofn
  >>>
  # usually much less than this, but can get >500M when highly parallel
  Int total_mem_mb = base_memory_mb + 4096
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File segmented = "${prefix}.consensusreadset.xml"
    File read_lengths_csv = "${prefix}.read_lengths.csv"
    File ligations_csv = "${prefix}.ligations.csv"
    File report = "${prefix}.summary.json"
    String segmented_bam_fn = read_string("bam.fofn")
    String non_passing_bam_fn = read_string("non_passing_bam.fofn")
    File? alarms = "alarms.json"
  }
}

task gather_sreads {
  input {
    File eid_ccs
    Array[File] segmented
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    set -e

    python3 -m pbcoretools.tasks.gather_segmented_reads \
      --log-level DEBUG \
      --log-file gather_segmented_reads.log \
      segmented.consensusreadset.xml \
      `readlink -f ${eid_ccs}` \
      ${sep=" " segmented}

    python3 -m pbcoretools.tasks.memory.get_dataset_size \
      --skip-counts \
      segmented.consensusreadset.xml
  }
  Int total_mem_mb = base_memory_mb + 1024
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File segmented = "segmented.consensusreadset.xml"
    File num_records_txt = "numrecords.txt"
  }
}

task pbreports_read_segmentation {
  input {
    File segmented
    # This task also gathers the remaining chunked outputs
    Array[File] skera_report
    Array[File] read_lengths_csv
    Array[File] ligations_csv
    # this is required for memory allocation
    Int number_of_sreads_m
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    set -e

    python3 -m pbreports.tasks.gather_skera_reports \
      --log-level DEBUG \
      --log-file gather_skara_outputs.log \
      skera_summary.report.json \
      ${sep=" " skera_report}

    python3 -m pbcoretools.tasks.gather \
      --log-level DEBUG \
      --log-file gather_read_lengths_csv.log \
      segmented.read_lengths.csv ${sep=" " read_lengths_csv}
    gzip segmented.read_lengths.csv

    python3 -m pbcoretools.tasks.gather \
      --log-level DEBUG \
      --log-file gather_ligations_csv.log \
      segmented.ligations.csv ${sep=" " ligations_csv}

    python3 -m pbreports.report.read_segmentation \
      --log-level ${log_level} \
      --log-file pbreports.log \
      `readlink -f ${segmented}` \
      skera_summary.report.json \
      segmented.read_lengths.csv.gz \
      segmented.ligations.csv \
      read_segmentation.report.json
  }
  # this is dominated by reading segmented.read_lengths.csv into a Pandas
  # dataframe, which contains 4 integers per S-read, but we use shorter types
  # internally to cut down memory overhead
  Int SCALE_FACTOR = 40
  Int total_mem_mb = base_memory_mb + 1024 + (number_of_sreads_m * SCALE_FACTOR)
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File report = "read_segmentation.report.json"
    Array[File?] plot_pngs = glob("*.png")
  }
}

workflow pb_segment_reads {
  input {
    File eid_ccs
    File? eid_barcode
    File? adapters_fasta
    Int add_memory_mb = 0

    String log_level = "INFO"
    Int nproc = 1
    # ignored
    Int max_nchunks = 1
    File? tmp_dir
  }

  call prepare_inputs {
    input:
      ccs_xml = eid_ccs,
      eid_barcode = eid_barcode,
      adapters_fasta = adapters_fasta,
      nproc = 1,
      log_level = log_level
  }

  scatter (chunk_xml in prepare_inputs.chunks) {
    call skera {
      input:
        ccs_xml = chunk_xml,
        adapters_fasta = prepare_inputs.adapters,
        nproc = nproc,
        base_memory_mb = add_memory_mb,
        log_level = log_level
    }
  }

  call gather_sreads {
    input:
      eid_ccs = eid_ccs,
      segmented = skera.segmented,
      base_memory_mb = add_memory_mb,
      log_level = log_level
  }
  Int number_of_sreads = read_int(gather_sreads.num_records_txt)
  Int number_of_sreads_m = number_of_sreads / 1048576

  call pbreports_read_segmentation {
    input:
      segmented = gather_sreads.segmented,
      skera_report = skera.report,
      read_lengths_csv = skera.read_lengths_csv,
      ligations_csv = skera.ligations_csv,
      number_of_sreads_m = number_of_sreads_m,
      base_memory_mb = add_memory_mb,
      log_level = log_level
  }

  output {
    File segmented_reads = gather_sreads.segmented
    Array[File] segmented_reads_bam = skera.segmented_bam_fn
    Array[File] non_passing_bam = skera.non_passing_bam_fn
    File report_read_segmentation = pbreports_read_segmentation.report
  }
}
