version 1.0

# can also be used for mapping_stats_ccs and mapping_stats_hgap variants
task mapping_stats {
  input {
    File mapped
    File? all_reads 
    String report_module = "mapping_stats"
    Float min_concordance = 70
    # RQ vs. concordance, CCS mapping only
    Boolean? show_calibration_plot

    # this is no longer parallel
    Int nproc = 1
    String log_level = "INFO"
    Int index_memory_gb
    Int base_memory_mb = 0
  }
  command {
    python3 \
      -m pbreports.report.${report_module} \
      --log-level ${log_level} \
      --min-concordance ${min_concordance} \
      ${true="--calibration-plot" false="" show_calibration_plot} \
      ${"--all-reads `readlink -f " + all_reads + "`"} \
      `readlink -f ${mapped}` \
      mapping_stats.report.json
  }
  Int total_mem_mb = (index_memory_gb + 1) * 1024 + base_memory_mb
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File report = "mapping_stats.report.json"
    Array[File?] plot_pngs = glob("*.png")
  }
}

task polished_assembly {
  input {
    File coverage_gff
    File consensus_fastq
    File? report_contigs

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 \
      -m pbreports.report.polished_assembly \
      --log-level ${log_level} \
      ${"--contigs-report " + report_contigs} \
      ${coverage_gff} \
      ${consensus_fastq} \
      polished_assembly.report.json
  }
  runtime {
    cpu: 1
    memory: "4GB"
  }
  output {
    File report = "polished_assembly.report.json"
  }
}

task subread_stats {
  input {
    File subreads

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 \
      -m pbreports.report.subread_stats \
      --log-level ${log_level} \
      `readlink -f ${subreads}` \
      subreads.report.json
  }
  runtime {
    cpu: 1
    memory: "12GB"
  }
  output {
    File report = "subreads.report.json"
  }
}

task gather_task_reports {
  input {
    Array[File] task_reports
    String output_file_name = "tasks_report.json"

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 -m pbreports.tasks.gather_task_reports \
      -o ${output_file_name} \
      ${sep=" " task_reports}
  }
  runtime {
    cpu: 1
    memory: "100MB"
  }
  output {
    # output file is "plain" JSON, not .report.json, to prevent auto-import
    File report = "${output_file_name}"
  }
}

task consensus_reports {
  input {
    File reference
    File coverage_gff
    File variants_gff
    Int nVariants = 100
    Int sortSize = 10000

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    set -e
    python3 \
      -m pbreports.report.summarize_consensus \
      --log-level ${log_level} \
      --variantsGff ${variants_gff} \
      -o consensus.gff \
      ${coverage_gff}
    python3 \
      -m pbreports.report.top_variants \
      --log-level ${log_level} \
      top_variants.report.json \
      ${variants_gff} \
      `readlink -f ${reference}`
    python3 \
      -m pbreports.report.variants \
      --log-level ${log_level} \
      variants.report.json \
      `readlink -f ${reference}` \
      consensus.gff ${variants_gff}
  }
  runtime {
    cpu: 1
    memory: "4GB"
  }
  output {
    File consensus_gff = "consensus.gff"
    File report_variants = "variants.report.json"
    File report_top_variants = "top_variants.report.json"
  }
}
