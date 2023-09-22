# the coverage tasks are run independently of the mapping step because they
# are blocking otherwise in workflows such as resequencing or basemods (due
# to the way Cromwell executes sub-workflows)

version 1.0

task summarize_coverage {
  input {
    File mapped
    File reference
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 \
      -m pbreports.report.summarize_coverage \
      --log-level ${log_level} \
      --log-file pbreports.log \
      ${mapped} \
      `readlink -f ${reference}` \
      coverage.gff
  }
  Int total_mem_mb = base_memory_mb + 8192
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File coverage_gff = "coverage.gff"
  }
}

task gc_coverage_plot {
  input {
    File reference
    File alignments
    Int region_size = 100

    Int nproc = 1
    String log_level = "INFO"
    Int mem_per_core = 2
    Int base_memory_mb = 0
  }
  command {
    python3 \
      -m pbreports.tasks.plot_gc_coverage \
      --log-level ${log_level} \
      --log-file pbreports.log \
      --region-size ${region_size} \
      ${alignments} \
      `readlink -f ${reference}` \
      -o gc_coverage_boxplot.png
  }
  # this task has been written to keep memory consumption to a reasonable
  # level even at human scale.  8GB is probably excessive but it's difficult
  # to predict right now
  Int total_mem_mb = base_memory_mb + 8192
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File? plot_png = "gc_coverage_boxplot.png"
    File? plot_png_thumb = "gc_coverage_boxplot_thumbnail.png"
  }
}

# can also be used for coverage_hgap report
task pbreports_coverage {
  input {
    File reference
    File coverage_gff
    File? gc_coverage_plot_png
    String report_module = "coverage"
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }

  command {
    python3 \
      -m pbreports.report.${report_module} \
      --log-level ${log_level} \
      --log-file pbreports.log \
      ${"--gc-cov-plot " + gc_coverage_plot_png} \
      `readlink -f ${reference}` \
      ${coverage_gff} \
      coverage.report.json
  }
  Int total_mem_mb = base_memory_mb + 8192
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File report = "coverage.report.json"
    Array[File?] plot_pngs = glob("*.png")
  }
}

task plot_target_coverage {
  input {
    File mapped
    File target_regions_bed
    Int base_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 -m pbreports.tasks.plot_target_coverage \
      --log-level ${log_level} \
      --log-file pbreports.log \
      --ignore-error \
      ${mapped} \
      -t ${target_regions_bed} \
      -o target_coverage_plot.png
  }
  # FIXME this is way too much
  Int total_mem_mb = base_memory_mb + 64 * 1024
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File? plot_png = "target_coverage_plot.png"
  }
}

task target_coverage {
  input {
    File mapped
    File target_regions_bed
    File? plot_png

    Int nproc = 1
    String log_level = "INFO"
    Int index_memory_gb
    Int base_memory_mb = 0
  }
  command {
    python3 -m pbreports.report.target_coverage \
      --log-level ${log_level} \
      --log-file pbreports.log \
      ${"--plot-png " + plot_png} \
      ${mapped} ${target_regions_bed} target_regions.report.json
  }
  Int total_mem_mb = index_memory_gb * 1024 + base_memory_mb
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File report = "target_regions.report.json"
  }
}

workflow coverage_reports {
  input {
    File mapped
    File reference
    File? target_regions_bed
    Boolean run_gc_coverage_plot = true
    Int index_memory_gb = 0
    Int base_memory_mb = 0

    Int nproc
    String log_level = "INFO"
  }

  call summarize_coverage {
    input:
      mapped = mapped,
      reference = reference,
      base_memory_mb = base_memory_mb,
      nproc = nproc
  }

  if (run_gc_coverage_plot) {
    call gc_coverage_plot {
      input:
        alignments = mapped,
        reference = reference,
        nproc = nproc,
        base_memory_mb = base_memory_mb,
        log_level = log_level
    }
  }

  call pbreports_coverage {
    input:
      reference = reference,
      coverage_gff = summarize_coverage.coverage_gff,
      gc_coverage_plot_png = gc_coverage_plot.plot_png,
      base_memory_mb = base_memory_mb,
      nproc = nproc
  }

  if (defined(target_regions_bed) && (target_regions_bed != "")) {
    call plot_target_coverage {
      input:
        mapped = mapped,
        target_regions_bed = select_first([target_regions_bed]),
        base_memory_mb = base_memory_mb,
        log_level = log_level
    }

    call target_coverage {
      input:
        mapped = mapped,
        target_regions_bed = select_first([target_regions_bed]),
        plot_png = plot_target_coverage.plot_png,
        index_memory_gb = index_memory_gb,
        base_memory_mb = base_memory_mb,
        log_level = log_level
    }
  }

  output {
    File coverage_gff = summarize_coverage.coverage_gff
    File report_coverage = pbreports_coverage.report
    File? report_target_coverage = target_coverage.report
    File? plot_target_coverage_png = plot_target_coverage.plot_png
  }
}
