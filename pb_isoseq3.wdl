# Parallelized Iso-Seq 3 workflow with optional minimap2 alignment

version 1.0

import "tasks/pbmm2.wdl" as mapping
import "tasks/memory.wdl" as memory

# Inspect the input dataset to determine whether it has been demultiplexed
# already (this is assumed to include 3' and 5' primers), and if so write
# the datastore and dataset files equivalent to lima output.  This assumes
# that datasets have gone through the "garden path" demultiplexing workflow
# that produces split BAM files.
task dataset_demux {
  input {
    File ccsxml
    Int memory_gb
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command <<<
    python3 <<EOF
    import subprocess
    import os.path as op
    import sys
    import pysam
    from pbcore.io import ConsensusReadSet
    from pbcommand.models.common import DataStoreFile, DataStore, PacBioAlarm
    ds = ConsensusReadSet(op.realpath("~{ccsxml}"), strict=True)
    if len(ds.externalResources) == 0:
      PacBioAlarm.dump_error(
        "alarms.json",
        "EmptyDatasetError",
        "The input ConsensusReadSet does not contain any BAM files, please select a valid PacBio dataset with a non-zero read count.",
        "The input ConsensusReadSet does not contain any BAM files, please select a valid PacBio dataset with a non-zero read count.",
        "EmptyDatasetError",
        "ERROR")
      sys.exit(1)
    elif ds.isBarcoded:
      files = []
      for i, resource in enumerate(ds.externalResources, start=1):
        ds_bc = ConsensusReadSet(resource.bam)
        ds_bc.name = ds.name + f" (sample {i})"
        ofn = f"sample_{i}.consensusreadset.xml"
        ds_bc.write(ofn)
        print(f"INFO: wrote sample {i} to {ofn}")
        files.append(DataStoreFile(ds_bc.uuid, "dataset_demux", ds.datasetType,
                                   op.abspath(ofn)))
      DataStore(files).write_json("samples.datastore.json")
    else:
      print("INFO: dataset is not barcoded")
      sample_names = set([rg.SampleName for rg in ds.readGroupTable])
      if len(sample_names) == 1 and sample_names != {"UnnamedSample"}:
        print("WARNING: SM tag is present, will reheader BAM")
        for ext_res, bam in zip(ds.externalResources, ds.resourceReaders()):
          print(f"INFO: reheadering {ext_res.bam}")
          header = dict(bam.peer.header)
          for rg in header["RG"]:
            if "SM" in rg:
              sample_name = rg.pop("SM")
              print(f"WARNING: removed sample tag '{sample_name}'")
          new_bam = op.basename(ext_res.bam)[:-3] + "reheadered.bam"
          with pysam.AlignmentFile(new_bam, "wb", header=header) as bam_out:
            for rec in bam.peer:
              bam_out.write(rec)
          subprocess.check_call(["samtools", "index", new_bam])
          subprocess.check_call(["pbindex", new_bam])
          ext_res.resourceId = new_bam
          ext_res.pbi = new_bam + ".pbi"
          ext_res.bai = new_bam + ".bai"
        ds.write("reheadered.consensusreadset.xml")
    EOF
  >>>
  Int total_mem_mb = base_memory_mb + 1024 * (memory_gb + 1)
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File? datastore = "samples.datastore.json"
    File? reheadered = "reheadered.consensusreadset.xml"
    Array[File?] datasets = glob("sample*.consensusreadset.xml")
    File? alarms = "alarms.json"
  }
}

task lima_isoseq {
  input {
    File ccsxml
    File barcodes
    Int base_memory_mb = 0

    Int nproc
    String log_level = "INFO"
  }
  command {
    set -e
    lima \
      -j ${nproc} \
      --isoseq \
      --peek-guess \
      --ignore-biosamples \
      --alarms alarms.json \
      `readlink -f ${ccsxml}` \
      `readlink -f ${barcodes}` \
      fl_transcripts.json
    mv fl_transcripts.json fl_transcripts.datastore.json
    if [ -f "fl_transcripts.lima.summary" ]; then
      mv fl_transcripts.lima.summary fl_transcripts.lima.summary.txt
    fi
    if [ -f "fl_transcripts.lima.guess" ]; then
      mv fl_transcripts.lima.guess fl_transcripts.lima.guess.txt
    fi
  }
  # in practice it looks even lower than this
  Int total_mem_mb = base_memory_mb + 2048
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File datastore = "fl_transcripts.datastore.json"
    Array[File] datasets = glob("fl_transcripts.*.consensusreadset.xml")
    File summary = "fl_transcripts.lima.summary.txt"
    File? infer_log = "fl_transcripts.lima.guess.txt"
  }
}

task refine {
  input {
    File fl_ccs
    File barcodes
    Boolean require_polya = true
    Boolean single_sample = false
    Int base_memory_mb = 0

    Int nproc
    String log_level = "INFO"
  }

  command {
    set -e

    isoseq3 \
      refine \
      --log-level ${log_level} \
      -j ${nproc} \
      --alarms alarms.json \
      ${true="--require-polya" false="" require_polya} \
      `readlink -f ${fl_ccs}` \
      `readlink -f ${barcodes}` \
      flnc.consensusreadset.xml \
      flnc.filter_summary.json \
      flnc.report.csv

    python3 -m pbcoretools.tasks.isoseq.collect_files \
      flnc.bam \
      --flnc-bam flnc.bam \
      --flnc-report flnc.report.csv \
      ${true="--single-sample" false="" single_sample} \
      --datastore isoseq_refine.datastore.json
  }
  Int total_mem_mb = 200 + base_memory_mb
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File flnc_ccs = "flnc.consensusreadset.xml"
    #File flnc_bam = "flnc.bam"
    #File report_csv = "flnc.report.csv"
    File datastore = "isoseq_refine.datastore.json"
    File report_filter_summary = "flnc.filter_summary.json"
  }
}

task cluster {
  input {
    File flnc_ccs
    # XXX unused?
    File barcodeset
    Int n_flnc_reads
    Int base_memory_mb = 0

    Int nproc
    String log_level = "INFO"
  }

  command {
    set -e
    isoseq3 \
      cluster \
      --alarms alarms.json \
      --log-level ${log_level} \
      --log-file isoseq_cluster.log \
      `readlink -f ${flnc_ccs}` \
      --use-qvs \
      unpolished.transcriptset.xml \
      -j ${nproc}
    dataset create --type TranscriptSet \
       --generateIndices \
       unpolished.hq.transcriptset.xml \
       unpolished.hq.bam
    dataset create --type TranscriptSet \
       --generateIndices \
       unpolished.lq.transcriptset.xml \
       unpolished.lq.bam
    python3 -m pbreports.tasks.isoseq_summarize \
      --log-level ${log_level} \
      unpolished.cluster_report.csv \
      unpolished.hq.transcriptset.xml \
      `readlink -f ${flnc_ccs}` \
      --csv-out hq_transcripts.fl_counts.csv
  }
  # best guess at empirical formula, needs more testing
  Int total_mem_kb = 1024 * (base_memory_mb + 2048 + (1024 * nproc)) + (n_flnc_reads * 40)
  runtime {
    cpu: nproc
    memory: "${total_mem_kb}KB"
  }
  output {
    File unpolished = "unpolished.transcriptset.xml"
    File report_csv = "unpolished.cluster_report.csv"
    File hq_bam = "unpolished.hq.bam"
    File lq_bam = "unpolished.lq.bam"
    File hq_transcripts = "unpolished.hq.transcriptset.xml"
    File lq_transcripts = "unpolished.lq.transcriptset.xml"
    File summary_csv = "hq_transcripts.fl_counts.csv"
  }
}

task bam2fasta_transcripts {
  input {
    File hq_transcripts
    File lq_transcripts
    File subreads
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }

  command {
    python3 -m pbcoretools.tasks.isoseq.bam2fasta_transcripts \
      --log-level ${log_level} \
      `readlink -f ${hq_transcripts}` \
      `readlink -f ${lq_transcripts}` \
      `readlink -f ${subreads}` \
      hq_transcripts.fasta \
      lq_transcripts.fasta
  }
  Int total_mem_mb = base_memory_mb + 1024
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File hq_fasta = "hq_transcripts.fasta"
    File lq_fasta = "lq_transcripts.fasta"
  }
}

task collapse {
  input {
    File hq_transcript_aln
    File flnc_ccs
    Int base_memory_mb = 0

    Float min_aln_coverage = 0.99
    Float min_aln_identity = 0.95
    Int max_fuzzy_junction = 5

    Int nproc
    String log_level = "INFO"
  }
  command <<<
    set -e
    DIRNAME=$(dirname $(readlink -f ~{hq_transcript_aln}))
    ln -s $DIRNAME/mapped.bam
    ln -s $DIRNAME/mapped.bam.bai
    isoseq3 collapse \
      --log-level ~{log_level} \
      --alarms alarms.json \
      -j ~{nproc} \
      --do-not-collapse-extra-5exons \
      --min-aln-coverage ~{min_aln_coverage} \
      --min-aln-identity ~{min_aln_identity} \
      --max-fuzzy-junction ~{max_fuzzy_junction} \
      $(readlink -f ~{hq_transcript_aln}) \
      $(readlink -f ~{flnc_ccs}) \
      collapse_isoforms.gff
  >>>
  Int total_mem_mb = base_memory_mb + 8192
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File fasta = "collapse_isoforms.fasta"  # CCSONLY transcripts do not contain QVs
    File gff = "collapse_isoforms.gff"
    File abundance = "collapse_isoforms.abundance.txt"
    File group = "collapse_isoforms.group.txt"
    File readstat = "collapse_isoforms.read_stat.txt"
    File report = "collapse_isoforms.report.json"
    File mapped_bam = "mapped.bam"
    File mapped_bam_bai = "mapped.bam.bai"
  }
}

task consolidate_transcripts {
  input {
    File subreads
    File hq_transcripts
    File lq_transcripts
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    set -e
    python3 -m pbcoretools.tasks.isoseq.consolidate_transcripts \
      --log-level ${log_level} \
      `readlink -f ${subreads}` \
      ${hq_transcripts} \
      ${lq_transcripts} \
      combined.hq.transcriptset.xml \
      combined.lq.transcriptset.xml \
      resources.hq.json \
      resources.lq.json
    dataset \
      --log-level ${log_level} \
      absolutize \
      combined.hq.transcriptset.xml
    dataset \
      --log-level ${log_level} \
      absolutize \
      combined.lq.transcriptset.xml
  }
  # again 2GB is probably an overestimate, since transcripts will never be
  # anywhere near as numerous as sequencing reads
  Int total_mem_mb = base_memory_mb + 2048
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File out_hq_transcripts = "combined.hq.transcriptset.xml"
    File out_lq_transcripts = "combined.lq.transcriptset.xml"
  }
}

# This is an appalling bit of glue code necessary to integrate this workflow
# with SMRT Link
task collect_cluster_outputs {
  input {
    File sample_bam
    File cluster_report_csv
    File barcode_overview_report
    File hq_fasta
    File lq_fasta
    File? hq_aln_bam
    File? hq_aln_bai
    File? collapse_fasta
    File? collapse_gff
    File? collapse_group
    File? collapse_abundance
    File? collapse_readstat
    Boolean single_sample = false
    Boolean all_samples = false
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 -m pbcoretools.tasks.isoseq.collect_files \
      --log-level ${log_level} \
      ${sample_bam} \
      --cluster-report-csv ${cluster_report_csv} \
      --barcode-overview-report ${barcode_overview_report} \
      --hq-fasta ${hq_fasta} \
      --lq-fasta ${lq_fasta} \
      ${"--hq-aln-bam " + hq_aln_bam} \
      ${"--hq-aln-bai " + hq_aln_bai} \
      ${"--collapse-fasta " + collapse_fasta} \
      ${"--collapse-gff " + collapse_gff} \
      ${"--collapse-group " + collapse_group} \
      ${"--collapse-abundance " + collapse_abundance} \
      ${"--collapse-readstat " + collapse_readstat} \
      ${true="--all-samples" false="" all_samples} \
      ${true="--single-sample" false="" single_sample} \
      --datastore isoseq_cluster.datastore.json
  }
  Int total_mem_mb = base_memory_mb + 200
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File datastore = "isoseq_cluster.datastore.json"
  }
}

# FIXME this is an awful python mushroom cap that should be pushed into the
# C++ layer whenever we have time
task demux_files {
  input {
    File flnc_ccs
    File? cluster_csv
    File? hq_transcripts
    File? lq_transcripts
    File? read_stats
    File? collapse_fasta
    String? datastore_prefix
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 -m pbcoretools.tasks.isoseq.joint_demux \
      --log-level ${log_level} \
      ${flnc_ccs} \
      ${"--cluster-csv " + cluster_csv} \
      ${"--transcripts-hq " + hq_transcripts} \
      ${"--transcripts-lq " + lq_transcripts} \
      ${"--collapse-fasta " + collapse_fasta} \
      ${"--read-stats " + read_stats} \
      ${"--datastore-prefix " + datastore_prefix}
  }
  Int total_mem_mb = base_memory_mb + 4096
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    Array[File] datastores = glob("*.datastore.json")
  }
}

task pbreports_isoseq_primers {
  input {
    # these are all *.datastore.json
    Array[File] barcoded
    File unbarcoded
    File barcodes
    Array[File] report_filter_summary
    Array[File] flnc_ccs
    Int base_memory_mb = 0

    # NOTE this doesn't really benefit from parallelization because there
    # should be only one barcoded dataset output right now
    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    set -e

    python3 -m pbcoretools.tasks.gather \
      --log-level DEBUG \
      lima_out.datastore.json \
      ${sep=" " barcoded}

    dataset \
      --log-level DEBUG \
      --strict \
      create \
      --trustCounts \
      flnc_ccs_all_samples.consensusreadset.xml \
      ${sep=" " flnc_ccs}

    python3 -m pbcoretools.tasks.gather \
      flnc_summary_all.report.json \
      ${sep=" " report_filter_summary}

    python3 -m pbreports.report.isoseq_primers \
      --log-level ${log_level} \
      -o isoseq_primers.report.json \
      lima_out.datastore.json \
      `readlink -f ${unbarcoded}` \
      `readlink -f ${barcodes}` \
      flnc_summary_all.report.json \
      flnc_ccs_all_samples.consensusreadset.xml
  }
  Int total_mem_mb = base_memory_mb + 4096 * 4
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File report = "isoseq_primers.report.json"
    Array[File?] plot_pngs = glob("*.png")
    Array[File?] plot_jsons = glob("*.json.gz")
    File all_flnc_ccs = "flnc_ccs_all_samples.consensusreadset.xml"
  }
}

task pbreports_isoseq3 {
  input {
    Array[File] hq_transcripts
    Array[File] lq_transcripts
    # these are optional for a reason, see comment in main workflow body
    File? cluster_csv
    File? flnc_ccs
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }

  command {
    set -e
    dataset --skipCounts create \
      --type TranscriptSet \
      hq_transcripts.transcriptset.xml \
      ${sep=" " hq_transcripts}
    dataset --skipCounts create \
      --type TranscriptSet \
      lq_transcripts.transcriptset.xml \
      ${sep=" " lq_transcripts}
    python3 -m pbreports.report.isoseq3 \
      --log-level ${log_level} \
      hq_transcripts.transcriptset.xml \
      lq_transcripts.transcriptset.xml \
      isoseq3.report.json \
      ${"--cluster-csv " + cluster_csv} \
      ${"--flnc-ccs " + flnc_ccs}
  }
  # being super paranoid here
  Int total_mem_mb = base_memory_mb + 8192
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File report = "isoseq3.report.json"
    Array[File?] plot_pngs = glob("*.png")
  }
}

task pbreports_isoseq_mapping {
  input {
    Array[File] collapse_fasta_or_fastq
    Array[File] collapse_report_json
    File? read_stats
    File? flnc_ccs
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    set -e
    python3 -m pbcoretools.tasks.gather \
      collapsed_isoforms.fasta \
      ${sep=' ' collapse_fasta_or_fastq}
    python3 -m pbcoretools.tasks.gather \
      all_samples_collapse.report.json \
      ${sep=' ' collapse_report_json}
    python3 -m pbreports.report.isoseq_mapping \
      --log-level ${log_level} \
      collapsed_isoforms.fasta \
      --collapse-report all_samples_collapse.report.json \
      ${"--read-stats " + read_stats} \
      ${"--flnc-ccs " + flnc_ccs} \
      --report-json isoseq_mapping.report.json
  }
  # FIXME this is again excessive, but safest for now
  Int total_mem_mb = base_memory_mb + 8192
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File report = "isoseq_mapping.report.json"
    Array[File?] plot_pngs = glob("*.png")
    File? abundance = "collapse_isoforms.fl_count.txt"
  }
}

workflow pb_isoseq3 {
  input {
    File eid_ccs
    File eid_barcode
    File? eid_ref_dataset

    File? biosamples_csv
    String? biosamples_csv_str
    Int filter_min_qv = 20
    String dataset_filters = ""

    # main window options
    Boolean run_clustering = true
    Boolean cluster_separately = false

    # Options exposed in smrtlink:
    Boolean refine_require_polya = true

    Int mapping_min_length = 50
    Float mapping_min_concordance = 95
    Float mapping_min_coverage = 99
    String? mapping_pbmm2_overrides
    Int add_memory_mb = 0

    Int isocollapse_max_fuzzy_junction = 5

    Int nproc = 1
    String log_level = "INFO"
    # this is completely ignored now
    Int max_nchunks = 0
    String? tmp_dir
  }

  call memory.get_input_sizes {
    input:
      bam_dataset = eid_ccs,
      fasta_dataset = eid_ref_dataset
  }

  call dataset_demux {
    input:
      ccsxml = eid_ccs,
      memory_gb = get_input_sizes.index_memory_gb,
      base_memory_mb = add_memory_mb,
      log_level = log_level
  }

  Boolean was_demuxed = defined(dataset_demux.datastore)
  if (was_demuxed) {
    Array[File] datasets_in_demuxed = select_all(dataset_demux.datasets)
  }
  if (defined(dataset_demux.reheadered)) {
    Array[File] datasets_reheadered = select_all([dataset_demux.reheadered])
  }
  Array[File] datasets_in = select_first([datasets_in_demuxed,
                                          datasets_reheadered,
                                          [eid_ccs]])

  scatter (sample in datasets_in) {
    call lima_isoseq {
      input:
        ccsxml = sample,
        barcodes = eid_barcode,
        base_memory_mb = add_memory_mb,
        nproc = nproc
    }
  }

  Array[File] processed_data = flatten(lima_isoseq.datasets)
  # the behavior of the clustering and mapping steps is very different when
  # there are multiple samples
  Boolean is_single_sample = length(processed_data) == 1

  scatter (sample in processed_data) {
    call refine {
      input:
        fl_ccs = sample,
        barcodes = eid_barcode,
        require_polya = refine_require_polya,
        base_memory_mb = add_memory_mb,
        nproc = nproc,
        single_sample = is_single_sample
    }
  }

  call pbreports_isoseq_primers {
    input:
      barcoded = lima_isoseq.datastore,
      unbarcoded = eid_ccs,
      barcodes = eid_barcode,
      report_filter_summary = refine.report_filter_summary,
      flnc_ccs = refine.flnc_ccs,
      base_memory_mb = add_memory_mb
  }

  if (run_clustering) {

    # We can do clustering either separately (e.g. for separate species) or
    # jointly (e.g. tissue samples).  If the latter, we use the merged FLNC
    # CCS output that was used as input for the "isoseq_primers" report
    if (cluster_separately) {
      Array[File] cluster_samples_1 = refine.flnc_ccs
    } # there is no 'else'
    if (!cluster_separately) {
      Array[File] cluster_samples_2 = [pbreports_isoseq_primers.all_flnc_ccs]
    }
    # in this scope, both _1 and _2 become Array[File]?
    Array[File] cluster_samples = select_first([cluster_samples_1,
                                                cluster_samples_2])

    scatter (flnc_ccs in cluster_samples) {
      # we need the number of FLNC reads for memory allocation
      call memory.get_input_sizes as get_flnc_size {
        input:
          bam_dataset = flnc_ccs
      }

      call cluster {
        input:
          flnc_ccs = flnc_ccs,
          barcodeset = eid_barcode,
          n_flnc_reads = get_flnc_size.bam_size,
          base_memory_mb = add_memory_mb,
          nproc = nproc
      }

      call bam2fasta_transcripts {
        input:
          hq_transcripts = cluster.hq_transcripts,
          lq_transcripts = cluster.lq_transcripts,
          subreads = eid_ccs,
      }

      call consolidate_transcripts {
        input:
          subreads = eid_ccs,
          hq_transcripts = cluster.hq_transcripts,
          lq_transcripts = cluster.lq_transcripts
      }

      # Optional mapping if reference is defined
      if (defined(eid_ref_dataset)) {
        File reference = select_first([eid_ref_dataset])
        Int genome_length_mb = select_first([get_input_sizes.genome_length_mb])

        call mapping.pbmm2_align {
          input:
            unmapped = consolidate_transcripts.out_hq_transcripts,
            reference = reference,
            aln_ext = ".transcriptalignmentset.xml",
            min_concordance = mapping_min_concordance,
            min_length = mapping_min_length,
            pbmm2_overrides = mapping_pbmm2_overrides,
            genome_length_mb = genome_length_mb,
            mem_scale_factor = 6,
            base_memory_mb = add_memory_mb,
            preset_mode = "ISOSEQ",
            nproc = nproc
        }

        call collapse {
          input:
            hq_transcript_aln = pbmm2_align.mapped,
            flnc_ccs = flnc_ccs,
            max_fuzzy_junction = isocollapse_max_fuzzy_junction,
            min_aln_identity = mapping_min_concordance / 100.0,
            min_aln_coverage = mapping_min_coverage / 100.0,
            base_memory_mb = add_memory_mb,
            nproc = nproc,
            log_level = log_level
        }
      }

      call collect_cluster_outputs {
        input:
          sample_bam = flnc_ccs,
          cluster_report_csv = cluster.report_csv,
          barcode_overview_report = cluster.summary_csv,
          hq_fasta = bam2fasta_transcripts.hq_fasta,
          lq_fasta = bam2fasta_transcripts.lq_fasta,
          hq_aln_bam = collapse.mapped_bam,
          hq_aln_bai = collapse.mapped_bam_bai,
          collapse_fasta = collapse.fasta,
          collapse_gff = collapse.gff,
          collapse_group = collapse.group,
          collapse_abundance = collapse.abundance,
          collapse_readstat = collapse.readstat,
          single_sample = is_single_sample,
          all_samples = !cluster_separately,
          base_memory_mb = add_memory_mb
      }
    }

    # This report behaves differently when the cluster report CSV and FLNC
    # CCS files are given as input - we only want that behavior in joint
    # clustering experiments.  We also need to "demux" the outputs again,
    # which unfortunately requires yet another python task
    Boolean is_joint_clustering = !(cluster_separately || is_single_sample)
    if (is_joint_clustering) {
      File cluster_csv_file = cluster.report_csv[0]
      File flnc_ccs_file = pbreports_isoseq_primers.all_flnc_ccs

      call demux_files as demux_cluster_files {
        input:
          flnc_ccs = flnc_ccs_file,
          cluster_csv = cluster_csv_file,
          hq_transcripts = cluster.hq_transcripts[0],
          lq_transcripts = cluster.lq_transcripts[0],
          datastore_prefix = "cluster_files",
          base_memory_mb = add_memory_mb,
          log_level = log_level
      }
    }
    call pbreports_isoseq3 {
      input:
        hq_transcripts = cluster.hq_transcripts,
        lq_transcripts = cluster.lq_transcripts,
        cluster_csv = cluster_csv_file,
        flnc_ccs = flnc_ccs_file,
        base_memory_mb = add_memory_mb,
        log_level = log_level
    }

    if (defined(eid_ref_dataset)) {
      if (is_joint_clustering) {
        File read_stats_file = select_first(collapse.readstat)

        call demux_files as demux_collapse_files {
          input:
            flnc_ccs = pbreports_isoseq_primers.all_flnc_ccs,
            read_stats = read_stats_file,
            collapse_fasta = collapse.fasta[0],
            datastore_prefix = "collapse_files",
            base_memory_mb = add_memory_mb,
            log_level = log_level
        }
      }
      call pbreports_isoseq_mapping {
        input:
          collapse_fasta_or_fastq = select_all(collapse.fasta),
          collapse_report_json = select_all(collapse.report),
          read_stats = read_stats_file,
          flnc_ccs = flnc_ccs_file,
          base_memory_mb = add_memory_mb,
          log_level = log_level
      }
    }
  }

  output {
    File report_isoseq_primers = pbreports_isoseq_primers.report
    File? report_isoseq = pbreports_isoseq3.report
    # this has the per-sample FLNC BAM and report CSVs
    Array[File] datastore_refine = refine.datastore
    # everything else per-sample
    Array[File]? datastore_cluster = collect_cluster_outputs.datastore
    Array[File]? datastore_demux_cluster = demux_cluster_files.datastores
    Array[File]? datastore_demux_collapse = demux_collapse_files.datastores
    File? report_isoseq_mapping = pbreports_isoseq_mapping.report
    File? collapse_abundance_joint = pbreports_isoseq_mapping.abundance
    #File? lima_summary = lima_isoseq.summary
    #File? lima_infer_log = lima_isoseq.infer_log
  }
}
