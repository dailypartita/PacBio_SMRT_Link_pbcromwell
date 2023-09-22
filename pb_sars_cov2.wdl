# Workflow for SARS-CoV-2 variant calling on HiFi reads.  The input MUST be
# separate per-sample CCS BAM files processed using the SMRT Link
# demultiplexing workflow or equivalent lima command.

version 1.0

# To avoid any confusion over DataSet XML filters this workflow is very strict
# about splitting the input dataset by BAM file, each of which is assumed to
# represent a single sample (as output by the pb_demux_* workflows).  The
# output will be a series of chunk XMLs, which will include a barcode quality
# filter if the BAMs contain barcoding annotations (the standard use case).
# We also need to collect FASTA files for tools that don't read our XML.
# (The barcodes are fine as is but since we need to resolve the absolute path
# of the FASTA file this is handled here as well.)
task collect_inputs {
  input {
    File eid_ccs
    File eid_ref_dataset
    File eid_ref_dataset_2
    File eid_barcode
    Int min_bq = 80

    Int nproc = 1
    String log_level = "INFO"
    File? tmp_dir
  }
  command {
    python3 -m pbcoretools.tasks.sars_cov2_inputs \
      --log-level DEBUG \
      --min-bq ${min_bq} \
      ${eid_ccs} ${eid_barcode} ${eid_ref_dataset} ${eid_ref_dataset_2}
  }
  runtime {
    cpu: 1
    memory: "1GB"
  }
  output {
    Array[File] chunks = glob("chunk.*.consensusreadset.xml")
    File genome_fasta = "genome.fasta"
    File genome_fasta_fai = "genome.fasta.fai"
    File pbaa_guide_fasta = "pbaa_guide.fasta"
    File pbaa_guide_fasta_fai = "pbaa_guide.fasta.fai"
    File primers_fasta = "primers.fasta"
  }
}

# This task is run in parallel over sample BAMs.  Because there are many
# ways for a sample to fail (including but not limited to sample contamination,
# poor coverage, or bugs), the command script is more fault-tolerant than
# usual.
task run_sample {
  input {
    File ccs_reads
    File primers_fasta
    File pbaa_guide_fasta
    File pbaa_guide_fasta_fai
    File genome_fasta
    File genome_fasta_fai
    Boolean neighbor_mode = true
    Int min_coverage = 4
    Float min_alt_freq = 0.5
    Int min_bq = 80

    Int nproc = 1
    String log_level = "INFO"
    File? tmp_dir
  }

  # These are constants
  Int min_score_lead = 10
  Int min_cluster_read_count = 2
  Int min_amplicon_length = 100
  String output_prefix = "output"

  # FIXME cdunn pointed out that the tests for empty or non-existent output
  # files is too tightly coupled to tool behavior and we should solve this
  # with traps instead.  it's not clear what we can do in a Cromwell script
  # but there is room for improvement here
  command {
    set -e
    set -x
    echo `pwd`

    # to remove any ambiguity we generate a new BAM file first; this will
    # have already been filtered by BQ >= 80
    dataset \
      --log-level DEBUG \
      --log-file dataset_consolidate.log \
      consolidate \
      ${ccs_reads} \
      input.ccs.bam input.consensusreadset.xml

    # the input BAM will be used to set the sample name in the gathered outputs
    echo `pwd`/input.ccs.bam > outputs.fofn

    # FIXME these tools can't read biosample info directly
    python3 <<EOF
    from pbcore.io import BamReader
    bam = BamReader("input.ccs.bam")
    biosample = bam.readGroupTable[0].SampleName.replace(" ", "_")
    with open("biosample.txt", "wt") as txt_out:
      txt_out.write(biosample)
    EOF

    ln -s ${pbaa_guide_fasta} pbaa_guide.fasta
    ln -s ${pbaa_guide_fasta_fai} pbaa_guide.fasta.fai
    ln -s ${genome_fasta} genome.fasta
    ln -s ${genome_fasta_fai} genome.fasta.fai
    ln -s ${primers_fasta} primers.fasta

    # remove and record amplification primers
    (lima \
      --log-level ${log_level} \
      --log-file trimmed_amplicons.log \
      -j ${nproc} \
      --ccs \
      --ignore-biosamples \
      --min-length ${min_amplicon_length} \
      --min-score-lead ${min_score_lead} \
      --min-score ${min_bq} \
      ${true="--neighbors" false="--different" neighbor_mode} \
      input.consensusreadset.xml primers.fasta \
      trimmed_amplicons.consensusreadset.xml) || \
      echo "ERROR: lima failed"

    if [ -s "trimmed_amplicons.bam" ]; then
      bam2fastq -u -o trimmed_amplicons trimmed_amplicons.bam
      ln -s trimmed_amplicons.bam output.trimmed_amplicons.bam
      samtools index output.trimmed_amplicons.bam
      echo `pwd`/output.trimmed_amplicons.bam >> outputs.fofn
    fi
    if [ -s "trimmed_amplicons.lima.summary" ]; then
      ln -s trimmed_amplicons.lima.summary output.trimmed_amplicons.lima.summary
      echo `pwd`/output.trimmed_amplicons.lima.summary >> outputs.fofn
    fi

    if [ -s "trimmed_amplicons.fastq" ]; then
      samtools faidx trimmed_amplicons.fastq

      # amplicon clustering - note that this can run successfully but
      # produce empty sequences e.g. for blank samples
      (pbaa cluster \
        -j ${nproc} \
        --log-level ${log_level} \
        --log-file pbaa.log \
        --trim-ends 5 \
        --min-cluster-read-count ${min_cluster_read_count} \
        pbaa_guide.fasta trimmed_amplicons.fastq pbaa_out) || \
        echo "ERROR: pbaa failed"

      # Map the primer-trimmed amplicon reads to the reference genome
      pbmm2 align \
        -j ${nproc} \
        --log-level ${log_level} \
        --log-file mapped_amplicons.log \
        --sort --preset HIFI \
        genome.fasta trimmed_amplicons.bam mapped_amplicons.bam

      # Get coverage metrics from samtools
      samtools mpileup \
        --min-BQ 1 \
        -f genome.fasta \
        -s mapped_amplicons.bam > mapped_amplicons.bam.mpileup
      samtools depth \
        -q 0 -Q 0 \
        mapped_amplicons.bam > mapped_amplicons.bam.depth
    else
      echo "ERROR: empty amplicons FASTQ"
    fi

    # If pbaa was successful, process the results
    if [ -s "pbaa_out_passed_cluster_sequences.fasta" ]; then
      (consensusVariants \
        --runName "`cat biosample.txt`" \
        --prefix pbaa_consensus \
        --read_info pbaa_out_read_info.txt \
        --hifiSupport \
        trimmed_amplicons.fastq genome.fasta \
        pbaa_out_passed_cluster_sequences.fasta > consensusVariants.log) || \
        echo "ERROR: consensusVariants failed"

      # convert pbaa output to VCF format
      (pbaa2vcf \
        --passOnly \
        --sampleCol runName \
        -o pbaa.vcf \
        pbaa_consensus_alleles.csv pbaa_consensus_variants.csv \
        genome.fasta > pbaa2vcf.log) || \
        echo "ERROR: pbaa2vcf failed"
    else
      echo "ERROR: no passed clustered sequences"
    fi

    if [ -s "pbaa_out_read_info.txt" ]; then
      # now add pretty colors to the mapped amplicons
      pbaa bampaint \
        --log-level ${log_level} \
        --log-file bampaint.log \
        pbaa_out_read_info.txt mapped_amplicons.bam \
        ${output_prefix}.aligned_reads.bam

      echo `pwd`/${output_prefix}.aligned_reads.bam >> outputs.fofn
    fi

    # get the final variants VCF
    if [ -s "pbaa.vcf" ]; then
      VCFCons \
        genome.fasta ${output_prefix} \
        --sample-name "`cat biosample.txt`" \
        --min_coverage ${min_coverage} \
        --min_alt_freq ${min_alt_freq} \
        --vcf_type pbaa \
        --input_depth mapped_amplicons.bam.depth \
        --input_vcf pbaa.vcf > vcfcons.log

      # Now map the consensus sequences to the reference genome
      pbmm2 align \
        -j ${nproc} \
        --log-level ${log_level} \
        --log-file ${output_prefix}.vcfcons.frag_aligned.log \
        --sort --preset HIFI \
        genome.fasta ${output_prefix}.vcfcons.frag.fasta \
        ${output_prefix}.vcfcons.frag_aligned.bam

      # this report essentially just annotates the existing outputs with the
      # sample name pulled from the CCS BAM header.  Inputs are the info.csv
      # and variants.csv from VCFCons and the lima counts file.  Outputs are
      # a report.json (corresponding to the info.csv) and two CSV files.
      python3 -m pbreports.report.sars_cov2_sample \
        --log-level ${log_level} \
        --log-file pbreports_sars_cov2_sample.log \
        --min-coverage ${min_coverage} \
        input.ccs.bam \
        --summary-csv ${output_prefix}.vcfcons.info.csv \
        --variants-tsv ${output_prefix}.vcfcons.variants.csv \
        --counts-tsv trimmed_amplicons.lima.counts \
        --report-out sample_variants.report.json \
        --counts-out trimmed_amplicons.lima_counts.csv \
        --variants-out sample_variants.csv

      # collect VCFCons outputs
      echo `pwd`/${output_prefix}.vcfcons.vcf >> outputs.fofn
      echo `pwd`/${output_prefix}.vcfcons.fasta >> outputs.fofn
      echo `pwd`/${output_prefix}.vcfcons.frag.fasta >> outputs.fofn
      echo `pwd`/${output_prefix}.vcfcons.frag_aligned.bam >> outputs.fofn
    else
      echo "biosample,path" > FAILED
      echo "`cat biosample.txt`,`pwd`" >> FAILED
    fi
  }
  runtime {
    cpu: nproc
    memory: "4GB"
  }
  output {
    File outputs_fofn = "outputs.fofn"
    # These may not exist if one of the steps failed
    File? counts = "trimmed_amplicons.lima_counts.csv"
    File? report = "sample_variants.report.json"
    File? variants = "sample_variants.csv"
    File? failure = "FAILED"
  }
}

task pbreports_sars_cov2_summary {
  input {
    Array[File] sample_reports
    Array[File] sample_counts
    Array[File] sample_variants

    Int nproc = 1
    String log_level = "INFO"
    File? tmp_dir
  }
  command {
    set -e

    # collect per-sample report JSON to make a table
    python3 -m pbcoretools.tasks.gather \
      --log-level ${log_level} \
      merged_samples.report.json \
      ${sep=" " sample_reports}

    # collect per-sample primer counts
    python3 -m pbcoretools.tasks.gather \
      --log-level ${log_level} \
      merged_sample_counts.csv \
      ${sep=" " sample_counts}

    # collect per-sample variants
    python3 -m pbcoretools.tasks.gather \
      --log-level ${log_level} \
      all.vcfcons.variants.csv \
      ${sep=" " sample_variants}

    # generate the summary report
    python3 -m pbreports.report.sars_cov2_summary \
      --log-level ${log_level} \
      merged_samples.report.json \
      merged_sample_counts.csv \
      summary.report.json \
      amplicon_counts.tsv.txt
  }
  runtime {
    cpu: 1
    memory: "2GB"
  }
  output {
    File report = "summary.report.json"
    File counts_tsv = "amplicon_counts.tsv.txt"
    File variants_csv = "all.vcfcons.variants.csv"
  }
}

task collect_failures {
  input {
    Array[File] failures

    Int nproc = 1
    String log_level = "INFO"
    File? tmp_dir
  }
  command {
    python3 -m pbcoretools.tasks.gather \
      --log-level ${log_level} \
      sample_failures.csv \
      ${sep=" " failures}
  }
  runtime {
    cpu: 1
    memory: "200MB"
  }
  output {
    File failures_csv = "sample_failures.csv"
  }
}

task collect_outputs {
  input {
    Array[File] outputs_fofns
    File? sample_failures_csv

    Int nproc = 1
    String log_level = "INFO"
    File? tmp_dir
  }
  String input_prefix = "output"
  String output_prefix = "all"
  command {
    set -e

    # outputs are optional if all samples fail
    touch ${output_prefix}.aligned_reads.bam.zip
    touch ${output_prefix}.vcfcons.vcf.zip
    touch ${output_prefix}.vcfcons.fasta.zip
    touch ${output_prefix}.vcfcons.frag.fasta.zip
    touch ${output_prefix}.vcfcons.frag_aligned.bam.zip
    touch ${output_prefix}.trimmed_amplicons.bam.zip
    touch ${output_prefix}.trimmed_amplicons.lima.summary.zip

    python3 -m pbcoretools.tasks.sars_cov2_outputs \
      --log-level ${log_level} \
      --input-prefix ${input_prefix} \
      --output-prefix ${output_prefix} \
      ${"--failures " + sample_failures_csv } \
      ${sep=" " outputs_fofns}

    # this is an optimization for Cromwell call caching - it's simpler to
    # declare the zip files as output File objects, but by just reading the
    # name we bypass blocking I/O on large runs
    echo `pwd`/${output_prefix}.aligned_reads.bam.zip > aligned_reads_zip.fofn
    echo `pwd`/${output_prefix}.vcfcons.vcf.zip > vcf_zip.fofn
    echo `pwd`/${output_prefix}.vcfcons.fasta.zip > fasta_zip.fofn
    echo `pwd`/${output_prefix}.vcfcons.frag.fasta.zip > frag_fasta_zip.fofn
    echo `pwd`/${output_prefix}.vcfcons.frag_aligned.bam.zip > aligned_frag_zip.fofn
    echo `pwd`/${output_prefix}.trimmed_amplicons.bam.zip > trimmed_amplicons_zip.fofn
    echo `pwd`/${output_prefix}.trimmed_amplicons.lima.summary.zip > lima_summary_zip.fofn
  }
  runtime {
    cpu: 1
    memory: "1GB"
  }
  output {
    String vcf_zip_fn = read_string("vcf_zip.fofn")
    String fasta_zip_fn = read_string("fasta_zip.fofn")
    String frag_fasta_zip_fn = read_string("frag_fasta_zip.fofn")
    String aligned_reads_zip_fn = read_string("aligned_reads_zip.fofn")
    String aligned_frag_zip_fn = read_string("aligned_frag_zip.fofn")
    String trimmed_amplicons_zip_fn = read_string("trimmed_amplicons_zip.fofn")
    String lima_summary_zip_fn = read_string("lima_summary_zip.fofn")
    # this should be relatively small
    File? errors_zip = "error_logs.zip"
  }
}

workflow pb_sars_cov2 {
  input {
    File eid_ccs
    # pbaa guide
    File eid_ref_dataset
    # genome
    File eid_ref_dataset_2
    File eid_barcode
    Boolean lima_neighbor_mode = true
    Int min_coverage = 4
    Float min_alt_freq = 0.5
    Int min_bq = 80

    Int nproc = 1
    String log_level = "INFO"
    File? tmp_dir
    # this is not actually respected
    Int max_nchunks = 1
  }

  call collect_inputs {
    input:
      eid_ccs = eid_ccs,
      eid_ref_dataset = eid_ref_dataset,
      eid_ref_dataset_2 = eid_ref_dataset_2,
      eid_barcode = eid_barcode,
      min_bq = min_bq,
      log_level = log_level
  }

  scatter (barcoded_ds in collect_inputs.chunks) {
    call run_sample {
      input:
        ccs_reads = barcoded_ds,
        primers_fasta = collect_inputs.primers_fasta,
        pbaa_guide_fasta = collect_inputs.pbaa_guide_fasta,
        pbaa_guide_fasta_fai = collect_inputs.pbaa_guide_fasta_fai,
        genome_fasta = collect_inputs.genome_fasta,
        genome_fasta_fai = collect_inputs.genome_fasta_fai,
        neighbor_mode = lima_neighbor_mode,
        min_coverage = min_coverage,
        min_alt_freq = min_alt_freq,
        min_bq = min_bq,
        log_level = log_level,
        nproc = 1
    }
  }

  # If there are one ore more FAILED sentinel files, gather them as CSV for
  # download and reporting
  Array[File] failures = select_all(run_sample.failure)
  if (length(failures) > 0) {
    call collect_failures {
      input:
        failures = failures,
        log_level = log_level
    }
  }

  # Only run the report if at least one sample succeeded
  Array[File] sample_reports = select_all(run_sample.report)
  if (length(sample_reports) > 0) {
    call pbreports_sars_cov2_summary {
      input:
        sample_reports = sample_reports,
        sample_counts = select_all(run_sample.counts),
        sample_variants = select_all(run_sample.variants),
        log_level = log_level
    }
  }

  call collect_outputs {
    input:
      outputs_fofns = run_sample.outputs_fofn,
      sample_failures_csv = collect_failures.failures_csv,
      log_level = log_level
  }

  output {
    File? report_sars_cov2 = pbreports_sars_cov2_summary.report
    File? counts_tsv = pbreports_sars_cov2_summary.counts_tsv
    File? variants_csv = pbreports_sars_cov2_summary.variants_csv
    # to SMRT Link these will just look like any other path.  note that they
    # may be completely empty files, but they will still be created
    String vcf_zip = collect_outputs.vcf_zip_fn
    String fasta_zip = collect_outputs.fasta_zip_fn
    String frag_fasta_zip = collect_outputs.frag_fasta_zip_fn
    String aligned_reads_zip = collect_outputs.aligned_reads_zip_fn
    String aligned_frag_zip = collect_outputs.aligned_frag_zip_fn
    String trimmed_amplicons_zip = collect_outputs.trimmed_amplicons_zip_fn
    String lima_summary_zip = collect_outputs.lima_summary_zip_fn
    File? sample_failures_csv = collect_failures.failures_csv
    File? errors_zip = collect_outputs.errors_zip
  }
}
