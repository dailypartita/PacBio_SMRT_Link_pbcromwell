# Repeat Analysis application

version 1.0

import "wf_mapping.wdl"
import "wf_prepare_input.wdl"

task process_target_regions {
  input {
    File mapped
    File reference
    File target_regions_bed
    Int target_ploidy = 2
    Int extract_flank = 100

    String log_level = "INFO"
    # this is ignored
    Int nproc = 1
  }
  command <<<
    set -ex
    # reminder: dollar-brace is bash interpolation, tilde-brace is cromwell

    #FOR EACH target in TARGETBED
    #
    #TARGETBED has 6 columns
    #chr start stop + three more:
    #target: name of target (str)
    #motifs: comma-separated list of motifs, e.g. CCG,CAG
    #        use first motif as "primary" for histograms
    #revcomp: boolean 0/1. sometimes the ref is the opposite
    #         strand compared to what is typically reported
    prefix="output"
    while read -r chr start stop target motifs revcomp; do
        echo "chr:$chr start:$start stop:$stop target:$target motifs:$motifs revcomp:$revcomp"
        #first motif/primary
        pMotif=$(echo ${motifs} | cut -d, -f1)
        #revcomp flag, opt in
        if [ "$revcomp" == 1 ]; then
            rc="-r"
        else
            rc=""
        fi

        region="${chr}:${start}-${stop}"
        extract_fq="${prefix}.${target}.fastq"
        #extract sequence of just the repeat expansion
        python3 -m pbrepeatanalysis.extractRegion \
          --log-level ~{log_level} --log-file extractRegion.${target}.log \
          ~{mapped} \
          `readlink -f ~{reference}` \
          "${region}" ${rc} \
          -f ~{extract_flank} \
          -o ${extract_fq}
        #waterfall plot
        python3 -m pbrepeatanalysis.waterfall \
          --log-level ~{log_level} --log-file waterfall.${target}.log \
          -f png \
          -i ${extract_fq} \
          -m ${motifs} \
          -o "${prefix}.${target}.waterfall.png"
        #motif count csv
        python3 -m pbrepeatanalysis.countMotifs \
          --log-level ~{log_level} --log-file countMotifs.${target}.log \
          -b \
          -i ${extract_fq} \
          -m ${motifs} \
          -o "${prefix}.${target}.counts.csv"
        #primary count histogram
        python3 -m pbrepeatanalysis.plotCounts \
          --log-level ~{log_level} --log-file plotCounts.${target}.log \
          -i ${extract_fq} \
          -m ${pMotif} \
          -o "${prefix}.${target}"
        #cluster
        python3 -m pbrepeatanalysis.clusterByRegion \
          --log-level ~{log_level} --log-file cluster.${target}.log \
          ~{mapped} \
          `readlink -f ~{reference}` \
          "${region}" ${rc} \
          --drop --smrtlink \
          -c ~{target_ploidy} \
          -m ${motifs} \
          -p "${prefix}.${target}"

        echo "done processing .bed record $target"
    done < ~{target_regions_bed}
    find `pwd` -name "*.waterfall.png" > waterfall.fofn
    find `pwd` -name "*.insertSize.png" > insertSize.fofn
    find `pwd` -name "*.motifcount.png" > motifcount.fofn
    find `pwd` -name "*.counts.png" > counts.fofn
    find `pwd` -name "*.summarySL.csv" > summary.fofn
  >>>
  runtime {
    cpu: 1
    memory: "8GB"
  }
  output {
    Array[String] waterfall_plots = read_lines("waterfall.fofn")
    Array[String] length_plots = read_lines("insertSize.fofn")
    Array[String] motif_plots = read_lines("motifcount.fofn")
    Array[String] counts_csv = read_lines("counts.fofn")
    Array[String] summary_csv = read_lines("summary.fofn")
    #Array[File] counts_plots = glob(???)
  }
}

task count_on_target {
  input {
    File mapped
    File target_regions_bed

    String log_level = "INFO"
  }
  command {
    python3 -m pbrepeatanalysis.countOnTarget \
      --log-level ${log_level} \
      ${mapped} ${target_regions_bed} \
      -o .
  }
  runtime {
    cpu: 1
    memory: "8GB"
  }
  output {
    File counts_csv = "onTargetCounts.csv"
  }
}

task pbreports_repeat_analysis {
  input {
    File target_regions_bed
    File counts_csv
    File? plot_png
    Array[String] summary_csv
    Array[String] waterfall_plots
    Array[String] length_plots
    Array[String] motif_plots
    Int target_ploidy = 2

    String log_level = "INFO"
  }
  command {
    python3 -m pbreports.report.repeat_analysis \
      --log-level ${log_level} \
      --log-file pbreports.log \
      -o repeat_analysis.report.json \
      ${"-p " + plot_png} \
      -c ${target_ploidy} \
      ${target_regions_bed} \
      ${counts_csv} \
      ${sep=" " summary_csv} \
      ${sep=" " waterfall_plots} \
      ${sep=" " length_plots} \
      ${sep=" " motif_plots}
  }
  runtime {
    cpu: 1
    memory: "200MB"
  }
  output {
    # this needs to fail gracefully if summary_csv is empty
    File? report = "repeat_analysis.report.json"
    #Array[File]? plot_pngs = glob("*.png")
    File? alarms = "alarms.json"
  }
}

workflow pb_repeat_analysis {
  input {
    File eid_ccs
    File eid_ref_dataset
    File target_regions_bed
    String mapping_pbmm2_overrides = "-r 10000 -E 0 -L 0.1 -c 0"
    Int target_ploidy = 2
    Int extract_flank = 100
    String dataset_filters = ""
    Int filter_min_qv = 20

    Int nproc = 1
    # this is ignored
    Int max_nchunks = 0
    String log_level = "INFO"
    String? tmp_dir
  }

  call wf_prepare_input.prepare_input {
    input:
      dataset_xml = eid_ccs,
      dataset_ext = ".consensusreadset.xml",
      dataset_filters = dataset_filters,
      filter_min_qv = filter_min_qv,
      nproc = 1,
      log_level = log_level
  }

  call wf_mapping.mapping {
    input:
      reads = prepare_input.reads_file,
      reference = eid_ref_dataset,
      target_regions_bed = target_regions_bed,
      preset_mode = "HiFi",
      pbmm2_overrides = mapping_pbmm2_overrides,
      run_coverage = true,
      run_gc_coverage_plot = false,
      mem_scale_factor = 6,
      nproc = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir
  }

  call count_on_target {
    input:
      mapped = mapping.mapped,
      target_regions_bed = target_regions_bed,
      log_level = log_level
  }

  call process_target_regions {
    input:
      mapped = mapping.mapped,
      reference = eid_ref_dataset,
      target_regions_bed = target_regions_bed,
      target_ploidy = target_ploidy,
      extract_flank = extract_flank,
      log_level = log_level,
      nproc = nproc
  }

  call pbreports_repeat_analysis {
    input:
      target_regions_bed = target_regions_bed,
      counts_csv = count_on_target.counts_csv,
      summary_csv = process_target_regions.summary_csv,
      waterfall_plots = process_target_regions.waterfall_plots,
      length_plots = process_target_regions.length_plots,
      motif_plots = process_target_regions.motif_plots,
      plot_png = mapping.plot_target_coverage_png,
      target_ploidy = target_ploidy,
      log_level = log_level
  }

  output {
    File mapped = mapping.mapped
    File? report_mapping_stats = mapping.report_mapping_stats
    File? coverage_gff = mapping.coverage_gff
    File? report_coverage = mapping.report_coverage
    #File? report_target_coverage = mapping.report_target_coverage
    File? mapped_bam_datastore = mapping.mapped_bam_datastore
    File target_counts_csv = count_on_target.counts_csv
    File? report_repeat_analysis = pbreports_repeat_analysis.report
  }
}
