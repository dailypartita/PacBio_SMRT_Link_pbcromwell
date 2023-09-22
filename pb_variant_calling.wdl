# HiFi Mapping and Variant Calling using DeepVariant 1.4 and pbsv

version 1.0

import "tasks/chunking.wdl"
import "tasks/pbreports.wdl"
import "wf_coverage_reports.wdl"
import "wf_deep_variant.wdl"
import "wf_pbsv_intl.wdl"
import "wf_whatshap.wdl"
import "wf_prepare_input.wdl"

task pbmm2_align_wgs {
  input {
    File unmapped
    File reference
    Int i_chunk = 0
    Float min_concordance = 70
    Int min_length = 50
    String? biosample_name
    String? pbmm2_overrides
    Int base_memory_mb

    Int nproc = 1
    String log_level = "DEBUG"
    String? tmp_dir

    Int genome_length_mb = 3300
  }
  String output_prefix = "mapped-${i_chunk}"
  command {
    set -vex

    pbmm2 \
      align \
      --log-level ${log_level} \
      --log-file pbmm2.log \
      --alarms alarms.json \
      -j ${nproc} \
      --sort \
      --preset HiFi \
      --min-gap-comp-id-perc ${min_concordance} \
      --min-length ${min_length} \
      --sample "${biosample_name}" \
      ${pbmm2_overrides} \
      `readlink -f ${reference}` \
      `readlink -f ${unmapped}` \
      ${output_prefix}.bam

    dataset create \
      --generateIndices \
      --type ConsensusAlignmentSet \
      --name "Mapped Reads (chunk ${i_chunk})" \
      ${output_prefix}.consensusalignmentset.xml \
      ${output_prefix}.bam
  }
  # we do not use precalculated indices in our references, so we need to index
  # on the fly which consumes more memory
  Int mem_scale_factor = 6
  Int total_mem_mb = base_memory_mb + 4096 + (genome_length_mb * mem_scale_factor) + (nproc * 256)
  runtime {
    cpu: nproc
    memory: "${total_mem_mb}MB"
  }
  output {
    File mapped_bam = "${output_prefix}.bam"
    File mapped_bam_bai = "${output_prefix}.bam.bai"
    File mapped_xml = "${output_prefix}.consensusalignmentset.xml"
  }
}

task get_output_prefix {
  input {
    File mapped_xml
    File reference_xml
    String? mapping_biosample_name
  }
  command {
    python3 <<EOF
    import os.path as op
    from pbcore.io import openDataSet, BamReader
    # saving the base reference dataset file name for later
    ref = openDataSet(op.realpath("${reference_xml}"))
    ref_name = op.basename(ref.externalResources[0].resourceId).split(".")[0]
    with open("reference.txt", "wt") as ref_out:
      ref_out.write(ref_name)
    # determine the actual biosample name
    mapped = openDataSet(op.realpath("${mapped_xml}"))
    if "${mapping_biosample_name}" != "":
      sample_names = ["${mapping_biosample_name}"]
    else:
      sample_names = set([])
      for r in mapped.externalResources:
        bam = BamReader(r.resourceId)
        for rg in bam.readGroupTable:
          sample_names.add(rg.SampleName)
      if len(sample_names) > 1:
        pass  # TODO
    with open("biosample.txt", "wt") as biosample_out:
      biosample_out.write(sorted(list(sample_names))[0].replace(" ", "_"))
    EOF
  }
  runtime {
    cpu: 1
    memory: "500MB"
  }
  output {
    String sample_name = read_string("biosample.txt")
    String reference_name = read_string("reference.txt")
  }
}

task collect_outputs {
  input {
    String sample_name
    String reference_name
    File haplotagged_bam
    File haplotagged_bam_bai
    File phased_variants_vcf
    File phased_variants_vcf_tbi
    File variants_gvcf
    File structural_variants_vcf
    File structural_variants_vcf_tbi
    File phasing_stats_tsv
    File bcftools_stats_txt

    Int nproc = 1
  }
  String prefix = "${sample_name}.${reference_name}"
  command {
    ln -s "${variants_gvcf}" "${prefix}.deepvariant.g.vcf.gz"
    ln -s "${phased_variants_vcf}" "${prefix}.deepvariant.phased.vcf.gz"
    ln -s "${phased_variants_vcf_tbi}" "${prefix}.deepvariant.phased.vcf.gz.tbi"
    ln -s "${haplotagged_bam}" "${prefix}.haplotagged.bam"
    ln -s "${haplotagged_bam_bai}" "${prefix}.haplotagged.bam.bai"
    ln -s "${structural_variants_vcf}" "${prefix}.pbsv.vcf.gz"
    ln -s "${structural_variants_vcf_tbi}" "${prefix}.pbsv.vcf.gz.tbi"
    ln -s "${phasing_stats_tsv}" "${prefix}.phasing_stats.tsv"
    ln -s "${bcftools_stats_txt}" "${prefix}.deepvariant.stats.txt"
  }
  runtime {
    cpu: 1
    memory: "50MB"
  }
  output {
    File final_haplotagged_bam = "${prefix}.haplotagged.bam"
    File final_haplotagged_bam_bai = "${prefix}.haplotagged.bam.bai"
    File final_variants_gvcf = "${prefix}.deepvariant.g.vcf.gz"
    File final_phased_variants_vcf = "${prefix}.deepvariant.phased.vcf.gz"
    File final_phased_variants_vcf_tbi = "${prefix}.deepvariant.phased.vcf.gz.tbi"
    File final_structural_variants_vcf = "${prefix}.pbsv.vcf.gz"
    File final_structural_variants_vcf_tbi = "${prefix}.pbsv.vcf.gz.tbi"
    File final_phasing_stats_tsv = "${prefix}.phasing_stats.tsv"
    File final_bcftools_stats_txt = "${prefix}.deepvariant.stats.txt"
  }
}

task merge_haplotagged_bams {
  input {
    Array[File] chunked_haplotagged_bams

    String log_level = "INFO"
    Int nproc = 1
  }
  command {
    pbmerge \
      --num-threads ${nproc} \
      -o deepvariant.haplotagged.bam \
      ${sep=" " chunked_haplotagged_bams}
    samtools index deepvariant.haplotagged.bam
  }
  runtime {
    cpu: nproc
    memory: "1GB"
  }
  output {
    File haplotagged_bam = "deepvariant.haplotagged.bam"
    File haplotagged_bam_bai = "deepvariant.haplotagged.bam.bai"
  }
}

task pbreports_variant_calling {
  input {
    File mapped_xml
    File variants_vcf
    File sv_vcf
    File sv_report
    File? sv_plot_png

    String log_level = "INFO"
    Int nproc = 1
  }
  command {
    python3 -m pbreports.report.variant_calling \
      --log-level ${log_level} \
      --log-file pbreports.log \
      -o variant_calling.report.json \
      ${variants_vcf} \
      ${sv_vcf} \
      ${sv_report} \
      ${mapped_xml} \
      ${"--sv-plot-png " + sv_plot_png}
  }
  runtime {
    cpu: 1
    memory: "512MB"
  }
  output {
    File report = "variant_calling.report.json"
    Array[File?] plot_pngs = glob("*.png")
  }
}

workflow pb_variant_calling {
  input {
    File eid_ccs
    File eid_ref_dataset
    Int filter_min_qv = 20
    # mapping parameters
    String? mapping_biosample_name
    String? mapping_pbmm2_overrides
    # too many BAM files causes problems later
    Int mapping_max_nchunks = 4

    # pbsv call parameters
    String pbsv_chunk_length = "1M"
    Int min_sv_length = 20
    String? pbsv_override_args
    # pbsv v2.8.0, equivalent to --hifi
    Int min_percent_reads = 10                # -P 10
    Int min_reads_one_sample = 3              # -O 3
    Int min_reads_all_samples = 3             # -A 3
    Int min_reads_per_strand_all_samples = 0  # -S 0

    String deepvariant_docker_image = "google/deepvariant:1.5.0"
    Boolean enable_gpu = false
    String whatshap_docker_image = "quay.io/biocontainers/whatshap:1.4--py39hc16433a_1"
    Int add_memory_mb = 0

    Int nproc = 1
    # in practice we should use as many chunks as possible without straining
    # cromwell, because the make_examples step is very slow and single-core;
    # however other sub-workflows like mapping are throttled
    Int max_nchunks = 96
    String log_level = "INFO"
    String? tmp_dir
  }
  Int BASE_MEMORY_APP = 1024
  Int base_memory_mb = BASE_MEMORY_APP + add_memory_mb

  Int NPROC_MAX_MERGE = 16

  Float mapping_min_concordance = 70
  Int mapping_min_length = 50

  call wf_prepare_input.prepare_input {
    input:
      dataset_xml = eid_ccs,
      dataset_ext = ".consensusreadset.xml",
      filter_min_qv = filter_min_qv,
      nproc = 1,
      log_level = log_level
  }

  Int max_nchunks_split = if (max_nchunks > mapping_max_nchunks) then mapping_max_nchunks else max_nchunks
  call chunking.split_dataset {
    input:
      ds_in = prepare_input.reads_file,
      max_nchunks = max_nchunks_split,
      split_mode = "zmws",
      log_level = log_level
  }

  scatter (i_chunk in range(split_dataset.nchunks)) {
    call pbmm2_align_wgs {
      input:
        unmapped = split_dataset.chunks[i_chunk],
        reference = eid_ref_dataset,
        min_concordance = mapping_min_concordance,
        min_length = mapping_min_length,
        biosample_name = mapping_biosample_name,
        pbmm2_overrides = mapping_pbmm2_overrides,
        i_chunk = i_chunk,
        base_memory_mb = add_memory_mb,
        log_level = log_level,
        nproc = nproc,
    }
  }

  call chunking.gather_datasets_lite as gather_mapped_xml {
    input:
      chunks = pbmm2_align_wgs.mapped_xml,
      dataset_ext = ".consensusalignmentset.xml",
      dataset_name = "'Mapped BAMs'",
      log_level = log_level
  }

  call pbreports.mapping_stats {
    input:
      mapped = gather_mapped_xml.gathered,
      all_reads = eid_ccs,
      min_concordance = mapping_min_concordance,
      index_memory_gb = 4,
      base_memory_mb = add_memory_mb,
      nproc = 1,
      log_level = log_level
  }

  call wf_coverage_reports.coverage_reports {
    input:
      mapped = gather_mapped_xml.gathered,
      reference = eid_ref_dataset,
      index_memory_gb = 4,
      base_memory_mb = add_memory_mb,
      run_gc_coverage_plot = true,
      nproc = nproc,
      log_level = log_level
  }

  call wf_deep_variant.wf_deep_variant as deep_variant {
    input:
      mapped_bams = pbmm2_align_wgs.mapped_bam,
      mapped_bams_bai = pbmm2_align_wgs.mapped_bam_bai,
      eid_ref_dataset = eid_ref_dataset,
      docker_image = deepvariant_docker_image,
      enable_gpu = enable_gpu,
      base_memory_mb = base_memory_mb,
      nproc = nproc,
      max_nchunks = max_nchunks,
      log_level = log_level
  }

  call wf_pbsv_intl.pbsv_intl as pbsv {
    input:
      mapped = gather_mapped_xml.gathered,
      reference = eid_ref_dataset,
      chunk_length = pbsv_chunk_length,
      min_sv_length = min_sv_length,
      min_percent_reads = min_percent_reads,
      min_reads_one_sample = min_reads_one_sample,
      min_reads_all_samples = min_reads_all_samples,
      min_reads_per_strand_all_samples = min_reads_per_strand_all_samples,
      pbsv_override_args = pbsv_override_args,
      hifi_mode = true,
      nproc = nproc,
      log_level = log_level,
      max_nchunks = max_nchunks
  }

  call wf_whatshap.wf_whatshap as whatshap {
    input:
      mapped_bams = pbmm2_align_wgs.mapped_bam,
      mapped_bams_bai = pbmm2_align_wgs.mapped_bam_bai,
      variants_vcf = deep_variant.variants_vcf,
      variants_vcf_tbi = deep_variant.variants_vcf_index,
      reference_fasta = deep_variant.reference_fasta,
      reference_fasta_fai = deep_variant.reference_fasta_fai,
      docker_image = whatshap_docker_image,
      log_level = log_level,
      nproc = nproc
  }

  Int nproc_merge = if (nproc > NPROC_MAX_MERGE) then NPROC_MAX_MERGE else nproc
  call merge_haplotagged_bams {
    input:
      chunked_haplotagged_bams = whatshap.haplotagged_bams,
      nproc = nproc_merge
  }

  call get_output_prefix {
    input:
      mapped_xml = gather_mapped_xml.gathered,
      reference_xml = eid_ref_dataset,
      mapping_biosample_name = mapping_biosample_name
  }

  call collect_outputs {
    input:
      sample_name = get_output_prefix.sample_name,
      reference_name = get_output_prefix.reference_name,
      phased_variants_vcf = whatshap.phased_variants_vcf,
      phased_variants_vcf_tbi = whatshap.phased_variants_vcf_tbi,
      variants_gvcf = deep_variant.variants_gvcf,
      structural_variants_vcf = pbsv.variants,
      structural_variants_vcf_tbi = pbsv.variants_index,
      haplotagged_bam = merge_haplotagged_bams.haplotagged_bam,
      haplotagged_bam_bai = merge_haplotagged_bams.haplotagged_bam_bai,
      phasing_stats_tsv = whatshap.phasing_stats_tsv,
      bcftools_stats_txt = deep_variant.bcftools_stats_txt,
      nproc = nproc
  }

  Array[File] sv_plot_pngs = select_all(pbsv.report_plot_pngs)
  if (length(sv_plot_pngs) > 0) {
    File sv_plot_png = select_first(pbsv.report_plot_pngs)
  }

  call pbreports_variant_calling {
    input:
      mapped_xml = gather_mapped_xml.gathered,
      variants_vcf = whatshap.phased_variants_vcf,
      sv_vcf = pbsv.variants,
      sv_report = pbsv.report,
      sv_plot_png = sv_plot_png,
      log_level = log_level
  }

  output {
    File mapped = gather_mapped_xml.gathered
    File report_mapping_stats = mapping_stats.report
    File report_coverage = coverage_reports.report_coverage
    File report_variant_calling = pbreports_variant_calling.report
    # final named outputs
    File structural_variants_vcf = collect_outputs.final_structural_variants_vcf
    File structural_variants_vcf_tbi = collect_outputs.final_structural_variants_vcf_tbi
    File bcftools_stats_txt = collect_outputs.final_bcftools_stats_txt
    File variants_gvcf = collect_outputs.final_variants_gvcf
    File phased_variants_vcf = collect_outputs.final_phased_variants_vcf
    File phased_variants_vcf_tbi = collect_outputs.final_phased_variants_vcf_tbi
    File haplotagged_bam = collect_outputs.final_haplotagged_bam
    File haplotagged_bam_bai = collect_outputs.final_haplotagged_bam_bai
    File phasing_stats_tsv = collect_outputs.final_phasing_stats_tsv
  }
}
