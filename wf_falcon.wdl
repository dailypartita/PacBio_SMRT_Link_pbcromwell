version 1.0

task task__0_rawreads__build {
  input {
    File config
    Array[File] fasta_array
    String topdir = "../.."
    String pypeflow_nproc = "1"
    String pypeflow_mb = "4000"
  }
  command <<<
    files_fn="~{write_lines(fasta_array)}"
    python3 -m falcon_kit.mains.dazzler --config-fn=~{config} --db-fn=./raw_reads.db  build --input-fofn-fn=${files_fn} --length-cutoff-fn=./length_cutoff
    # TODO: Verify that db exists.
    #ln -sf ./length_cutoff length_cutoff
    echo "$(pwd)/raw_reads.db" >| dazzler_db_fn
    echo "$(pwd)/.raw_reads.idx" >| dazzler_idx_fn
    echo "$(pwd)/.raw_reads.bps" >| dazzler_bps_fn
    echo "$(pwd)/.raw_reads.dust.anno" >| dazzler_dust_anno_fn
    echo "$(pwd)/.raw_reads.dust.data" >| dazzler_dust_data_fn
  >>>
  output {
    File dazzler_db_txt = "dazzler_db_fn"
    File dazzler_idx_txt = "dazzler_idx_fn"
    File dazzler_bps_txt = "dazzler_bps_fn"
    File dazzler_dust_anno_txt = "dazzler_dust_anno_fn"
    File dazzler_dust_data_txt = "dazzler_dust_data_fn"
    File length_cutoff = "length_cutoff"
  }
}

task task__0_rawreads__tan_split {
  input {
    File config
    String dazzler_db
    String dazzler_idx
    String dazzler_bps
    #String dazzler_dust_anno = ".raw_reads.dust.anno"
    #String dazzler_dust_data = ".raw_reads.dust.data"
    #File length_cutoff
    #String pypeflow_nproc = "1"
    #String wildcards = "dal0_id"
  }
  command <<<
    python3 -m falcon_kit.mains.dazzler --config=~{config} --db=~{dazzler_db}                            tan-split --split=./tan-uows.json --bash-template=./bash_template.sh
  >>>
  output {
    File uows = 'all-units-of-work.tar'
  }
}

task task__0_rawreads__tan_scatter {
  input {
    File uows
    Int max_nchunks = 256
  }
  command <<<
    python3 -m falcon_kit.mains.generic_scatter_uows_tar --all=~{uows} --nchunks=~{max_nchunks} --pattern='./some-units-of-work.%.tar'
  >>>
  output {
    Array[File] chunks = glob("some-units-of-work.*.tar")
  }
}

task task__0_rawreads__tan_apply {
  input {
    File chunk
    File config
    String dazzler_db
    String dazzler_idx
    String dazzler_bps
    #String dazzler_dust_anno
    #String dazzler_dust_data
    Int nproc = 4
  }
  command <<<
    set -vex
    python3 -m falcon_kit.mains.cromwell_symlink ~{config} ~{dazzler_db} ~{dazzler_idx} ~{dazzler_bps}
    python3 -m falcon_kit.mains.cromwell_run_uows_tar --nproc=~{nproc} --nproc-per-uow=4 --uows-tar-fn=~{chunk} --tool=datander
    # glob seems to fail on dotfiles, so we rename them first.
    ls -la uow-*/.raw_reads.*.tan.anno
    ls -la uow-*/.raw_reads.*.tan.data
    python3 -m falcon_kit.mains.cromwell_undot --pattern 'uow-*/.*.*.tan.anno' --prefix dot
    python3 -m falcon_kit.mains.cromwell_undot --pattern 'uow-*/.*.*.tan.data' --prefix dot
    find `pwd`/uow-* -name "dot.*.*.tan.anno" > dazzler_tan_anno.fofn
    find `pwd`/uow-* -name "dot.*.*.tan.data" > dazzler_tan_data.fofn
    #find . -name '.*.*.tan.anno' | xargs -I XXX mv -sf XXX dotXXX
    #find . -name '.*.*.tan.data' | xargs -I XXX ln -sf XXX dotXXX
    #ls -larth .
    #Catrack -vdf raw_reads tan
  >>>
  runtime {
    cpu: nproc
  }
  output {
    File dazzler_tan_anno_fofn = "dazzler_tan_anno.fofn"
    File dazzler_tan_data_fofn = "dazzler_tan_data.fofn"
  }
}

task task__0_rawreads__tan_combine {
  input {
    File config
    String dazzler_db
    String dazzler_idx
    String dazzler_bps
    Array[File] tan_anno_array_of_fofn
    Array[File] tan_data_array_of_fofn
  }

  command <<<
    set -vex
    anno_files_fn=anno_files.fofn
    data_files_fn=data_files.fofn
    cat ~{sep=" " tan_anno_array_of_fofn} > ${anno_files_fn}
    cat ~{sep=" " tan_data_array_of_fofn} > ${data_files_fn}
    cat ${anno_files_fn}

    python3 -m falcon_kit.mains.dazzler --config=~{config} --db=~{dazzler_db}  track-combine --track tan --anno ${anno_files_fn} --data ${data_files_fn}
    # Note: track-combine, not tan-combine, which was for the UOW workflow
    echo "$(pwd)/.raw_reads.tan.anno" >| dazzler_tan_anno_fn
    echo "$(pwd)/.raw_reads.tan.data" >| dazzler_tan_data_fn
  >>>
  output {
    File dazzler_tan_anno_txt = "dazzler_tan_anno_fn"
    File dazzler_tan_data_txt = "dazzler_tan_data_fn"
  }
}

task task__0_rawreads__daligner_split {
  input {
    File config
    String dazzler_db = "raw_reads.db"
    String dazzler_idx = ".raw_reads.idx"
    String dazzler_bps = ".raw_reads.bps"
    String dazzler_dust_anno = ".raw_reads.dust.anno"
    String dazzler_dust_data = ".raw_reads.dust.data"
    String dazzler_tan_anno = ".raw_reads.tan.anno"
    String dazzler_tan_data = ".raw_reads.tan.data"
    File length_cutoff
    String pypeflow_nproc = "1"
    String wildcards = "dal0_id"
  }
  command <<<
    python3 -m falcon_kit.mains.dazzler --config=~{config} --db=~{dazzler_db} --nproc=~{pypeflow_nproc}  daligner-split --wildcards=~{wildcards} --length-cutoff-fn=~{length_cutoff} --split-fn=./all-units-of-work.json --bash-template-fn=./daligner_bash_template.sh
    #python3 -m falcon_kit.mains.generic_tar_uows --all=./all-units-of-work.json --nchunks=256 --pattern=./units-of-work.%.tar
    #python3 -m falcon_kit.mains.generic_scatter_one_uow --all-uow-list-fn=./all-units-of-work.json --one-uow-list-fn=some-units-of-work.json --split-idx=0
  >>>
  output {
    File uows = 'all-units-of-work.tar'
  }
}

task task__0_rawreads__daligner_scatter {
  input {
    File uows
    Int max_nchunks = 256
  }
  command <<<
    python3 -m falcon_kit.mains.generic_scatter_uows_tar --all=~{uows} --nchunks=~{max_nchunks} --pattern='./some-units-of-work.%.tar'
  >>>
  output {
    Array[File] chunks = glob("some-units-of-work.*.tar")
  }
}

task task__0_rawreads__daligner_apply {
  input {
    File chunk
    File config
    String dazzler_db
    String dazzler_idx
    String dazzler_bps
    String dazzler_dust_anno
    String dazzler_dust_data
    String dazzler_tan_anno
    String dazzler_tan_data
    Int nproc = 4
  }
  command <<<
    set -e
    python3 -m falcon_kit.mains.cromwell_symlink ~{config} ~{dazzler_db} ~{dazzler_idx} ~{dazzler_bps} ~{dazzler_dust_anno} ~{dazzler_dust_data} ~{dazzler_tan_anno} ~{dazzler_tan_data}
    #echo 'python3 -m falcon_kit.mains.dazzler --config=~{config} --db=~{dazzler_db}  daligner-apply --script=./run_daligner.sh --job-done={output.job_done}' > ./bash_template.sh

    python3 -m falcon_kit.mains.cromwell_run_uows_tar --nproc=~{nproc} --nproc-per-uow=4 --uows-tar-fn=~{chunk} --tool=daligner
    find `pwd`/uow-* -maxdepth 1 -name "*.las" > chunk_las_files.fofn
  >>>
  runtime {
    cpu: nproc
  }
  output {
    File result = "chunk_las_files.fofn"
  }
}

task task__0_rawreads__daligner_las_merge {
  input {
    Array[File] las_array_of_fofns
    File config
    Int nproc = 1
  }

  command <<<
    set -vex
    files_fn=las_files.fofn
    cat ~{sep=" " las_array_of_fofns} > ${files_fn}
    python3 -m falcon_kit.mains.cromwell_write_json --lines-fn=${files_fn} --json-fn=./gathered-las.json

    python3 -m falcon_kit.mains.cromwell_symlink ~{config}

    python3 -m falcon_kit.mains.dazzler --config=~{basename(config)}                  merge-split --db-prefix=raw_reads --las-paths=./gathered-las.json --wildcards=mer0_id --split-fn=all-units-of-work.json --bash-template-fn=las-merge-bash-template.sh

    python3 -m falcon_kit.mains.generic_run_units_of_work --nproc=~{nproc} --bash-template-fn=./las-merge-bash-template.sh --units-of-work-fn=./all-units-of-work.json --results-fn=./merged-las-paths.json
  >>>
  runtime {
    cpu: nproc
  }
  output {
    Array[File] las_array = glob("uow-*/*.las")
  }
}

task task__0_rawreads__cns_apply {
  input {
    #File chunk
    File las
    File config
    String dazzler_db
    String dazzler_idx
    String dazzler_bps
    File length_cutoff
    Int nproc = 8
  }
  command <<<
  set -vex
  python3 -m falcon_kit.mains.cromwell_symlink ~{config} ~{dazzler_db} ~{dazzler_idx} ~{dazzler_bps} ~{length_cutoff} ~{las}

  python3 -m falcon_kit.mains.consensus_task --nproc=~{nproc} --las-fn=~{basename(las)} --db-fn=~{basename(dazzler_db)} --length-cutoff-fn=~{basename(length_cutoff)} --config-fn=~{basename(config)} --fasta-fn=consensus.chunk.fasta
  >>>
  runtime {
    cpu: nproc
  }
  output {
    File fasta = "consensus.chunk.fasta"
  }
}

task task__0_rawreads__report {
  input {
    Array[File] fasta_array
    File config
    String dazzler_db
    String dazzler_idx
    String dazzler_bps
    File length_cutoff
  }
  command <<<
  set -vex
  python3 -m falcon_kit.mains.cromwell_symlink ~{config} ~{dazzler_db} ~{dazzler_idx} ~{dazzler_bps} ~{length_cutoff}

  files_fn="~{write_lines(fasta_array)}"
  ln -sf ${files_fn} input_preads.fofn

  python3 -m falcon_kit.mains.task_report_pre_assembly --config-fn=~{basename(config)} --length-cutoff-fn=~{basename(length_cutoff)} --raw-reads-db-fn=~{basename(dazzler_db)} --preads-fofn-fn=./input_preads.fofn --pre-assembly-report-fn=./pre_assembly_stats.json

  python3 -m pbreports.report.preassembly pre_assembly_stats.json preassembly.report.json
  >>>
  output {
    File stats = "pre_assembly_stats.json"
    File report = "preassembly.report.json"
  }
}

task task__1_preads_ovl__build {
  input {
    File config
    Array[File] fasta_array
    String topdir = "../.."
    String pypeflow_nproc = "1"
    String pypeflow_mb = "4000"
  }
  command <<<
    files_fn="~{write_lines(fasta_array)}"
    python3 -m falcon_kit.mains.dazzler --config-fn=~{config} --db-fn=./preads.db  build --input-fofn-fn=${files_fn} --length-cutoff-fn=./length_cutoff
    # TODO: Verify that db exists.
    #ln -sf ./length_cutoff length_cutoff
    echo "$(pwd)/preads.db" >| dazzler_db_fn
    echo "$(pwd)/.preads.idx" >| dazzler_idx_fn
    echo "$(pwd)/.preads.bps" >| dazzler_bps_fn
    echo "$(pwd)/.preads.dust.anno" >| dazzler_dust_anno_fn
    echo "$(pwd)/.preads.dust.data" >| dazzler_dust_data_fn
  >>>
  output {
    File dazzler_db_txt = "dazzler_db_fn"
    File dazzler_idx_txt = "dazzler_idx_fn"
    File dazzler_bps_txt = "dazzler_bps_fn"
    File dazzler_dust_anno_txt = "dazzler_dust_anno_fn"
    File dazzler_dust_data_txt = "dazzler_dust_data_fn"
    File length_cutoff = "length_cutoff"
  }
}

task task__1_preads_ovl__daligner_split {
  input {
    File config
    String dazzler_db
    String dazzler_idx
    String dazzler_bps
    String dazzler_dust_anno
    String dazzler_dust_data
    File length_cutoff
    String pypeflow_nproc = "1"
    String wildcards = "dal0_id"
  }
  command <<<
    python3 -m falcon_kit.mains.dazzler --config=~{config} --db=~{dazzler_db} --nproc=~{pypeflow_nproc}  daligner-split --wildcards=~{wildcards} --length-cutoff-fn=~{length_cutoff} --split-fn=./all-units-of-work.json --bash-template-fn=./daligner_bash_template.sh
    #python3 -m falcon_kit.mains.generic_tar_uows --all=./all-units-of-work.json --nchunks=256 --pattern=./units-of-work.%.tar
    #python3 -m falcon_kit.mains.generic_scatter_one_uow --all-uow-list-fn=./all-units-of-work.json --one-uow-list-fn=some-units-of-work.json --split-idx=0
  >>>
  output {
    File uows = 'all-units-of-work.tar'
  }
}

task task__1_preads_ovl__daligner_scatter {
  input {
    File uows
    Int max_nchunks = 256
  }
  command <<<
    python3 -m falcon_kit.mains.generic_scatter_uows_tar --all=~{uows} --nchunks=~{max_nchunks} --pattern='./some-units-of-work.%.tar'
  >>>
  output {
    Array[File] chunks = glob("some-units-of-work.*.tar")
  }
}

task task__1_preads_ovl__daligner_apply {
  input {
    File chunk
    File config
    String dazzler_db
    String dazzler_idx
    String dazzler_bps
    String dazzler_dust_anno
    String dazzler_dust_data
    #String dazzler_tan_anno
    #String dazzler_tan_data
    Int nproc = 4
  }
  command <<<
    set -e
    python3 -m falcon_kit.mains.cromwell_symlink ~{config} ~{dazzler_db} ~{dazzler_idx} ~{dazzler_bps} ~{dazzler_dust_anno} ~{dazzler_dust_data}
    #echo 'python3 -m falcon_kit.mains.dazzler --config=~{config} --db=~{dazzler_db}  daligner-apply --script=./run_daligner.sh --job-done={output.job_done}' > ./bash_template.sh

    python3 -m falcon_kit.mains.cromwell_run_uows_tar --nproc=~{nproc} --nproc-per-uow=4 --uows-tar-fn=~{chunk} --tool=daligner
    find `pwd`/uow-* -name "*.las" > chunk_las_files.fofn
  >>>
  runtime {
    cpu: nproc
  }
  output {
    #Array[File] result = glob("uow-*/*.las")
    File result = "chunk_las_files.fofn"
  }
}

task task__1_preads_ovl__daligner_las_merge {
  input {
    Array[File] las_array_of_fofns
    File config
    Int nproc = 1
  }

  command <<<
    set -vex
    files_fn=las_files.fofn
    cat ~{sep=" " las_array_of_fofns} > ${files_fn}
    python3 -m falcon_kit.mains.cromwell_write_json --lines-fn=${files_fn} --json-fn=./gathered-las.json

    python3 -m falcon_kit.mains.cromwell_symlink ~{config}

    python3 -m falcon_kit.mains.dazzler --config=~{basename(config)}                  merge-split --db-prefix=preads --las-paths=./gathered-las.json --wildcards=mer0_id --split-fn=all-units-of-work.json --bash-template-fn=las-merge-bash-template.sh

    python3 -m falcon_kit.mains.generic_run_units_of_work --nproc=~{nproc} --bash-template-fn=./las-merge-bash-template.sh --units-of-work-fn=./all-units-of-work.json --results-fn=./merged-las-paths.json
    find `pwd`/uow-* -name "*.las" > las_merged_files.fofn
  >>>
  runtime {
    cpu: nproc
  }
  output {
    File las_fofn = "las_merged_files.fofn"
  }
}

task task__1_preads_ovl__db2falcon {
  input {
    String dazzler_db
    String dazzler_idx
    String dazzler_bps
  }
  command <<<
    set -vex
    python3 -m falcon_kit.mains.cromwell_symlink ~{dazzler_db} ~{dazzler_idx} ~{dazzler_bps}
    time DB2Falcon -U ~{basename(dazzler_db)}
    [ -f ./preads4falcon.fasta ] || exit 1
  >>>
  output {
    File preads4falcon_fasta = './preads4falcon.fasta'
  }
}

task task__2_asm_falcon {
  input {
    String dazzler_db
    String dazzler_idx
    String dazzler_bps
    File preads4falcon_fasta
    File las_fofn
    String overlap_filtering_setting
    String length_cutoff_pr
    String fc_ovlp_to_graph_option
    Int nproc
  }
  command <<<
    set -vex
    python3 -m falcon_kit.mains.cromwell_symlink ~{dazzler_db} ~{dazzler_idx} ~{dazzler_bps} ~{preads4falcon_fasta}

    files_fn="~{las_fofn}"
    python3 -m falcon_kit.mains.cromwell_write_json --lines-fn=${files_fn} --json-fn=./las_fofn.json

    # Given, las_fofn.json,
    # write preads.m4:

    #overlap_filtering_setting='--max-diff 10000 --max-cov 100000 --min-cov 1 --min-len 1 --bestn 1000 --n-core 0'
    #length_cutoff_pr='1'

    falconc m4filt-falconRunner -n ~{nproc} --db "~{basename(dazzler_db)}" --las-json ./las_fofn.json --min-len ~{length_cutoff_pr} --out preads.m4 --filter-log filter.log
    cat filter.log

    time falconc m4filt-contained ~{fc_ovlp_to_graph_option} --in preads.m4 --out preads.filtered.m4

    cat >| multiple-python-calls.sh << EOF

    # Given preads.filtered.m4,
    # write sg_edges_list, c_path, utg_data, ctg_paths.
    time python3 -m falcon_kit.mains.ovlp_to_graph ~{fc_ovlp_to_graph_option} --overlap-file preads.filtered.m4 >| fc_ovlp_to_graph.log

    # Given sg_edges_list, utg_data, ctg_paths, preads4falcon.fasta,
    # write p_ctg.fasta and a_ctg_all.fasta,
    # plus a_ctg_base.fasta, p_ctg_tiling_path, a_ctg_tiling_path, a_ctg_base_tiling_path:
    time python3 -m falcon_kit.mains.graph_to_contig

    # Given a_ctg_all.fasta, write a_ctg.fasta:
    time python3 -m falcon_kit.mains.dedup_a_tigs >| a_ctg.fasta

    # Given a_ctg.fasta and a_ctg_all_tiling_path, write a_ctg_tiling_path:
    time python3 -m falcon_kit.mains.dedup_a_tp >| a_ctg_tiling_path

    # Collect all info needed to format the GFA-1 and GFA-2 representations of
    # the assembly graphs.
    time python3 -m falcon_kit.mains.collect_pread_gfa --preads-ovl preads.m4 >| asm.gfa.json
    time python3 -m falcon_kit.mains.collect_pread_gfa --preads-ovl preads.m4 --add-string-graph >| sg.gfa.json
    time python3 -m falcon_kit.mains.collect_contig_gfa >| contig.gfa.json

    # Output the assembly pread graph.
    time python3 -m falcon_kit.mains.gen_gfa_v1 asm.gfa.json >| asm.gfa
    time python3 -m falcon_kit.mains.gen_gfa_v2 asm.gfa.json >| asm.gfa2

    # Output the string graph.
    time python3 -m falcon_kit.mains.gen_gfa_v1 sg.gfa.json >| sg.gfa
    time python3 -m falcon_kit.mains.gen_gfa_v2 sg.gfa.json >| sg.gfa2

    # Output the contig graph with associate contigs attached to each primary contig.
    time python3 -m falcon_kit.mains.gen_gfa_v2 contig.gfa.json >| contig.gfa2

EOF

    python3 -m falcon_kit.mains.run_python_modules multiple-python-calls.sh

    #rm -f ./preads4falcon.fasta
  >>>
  runtime {
    cpu: nproc
  }
  output {
    File p_ctg_fasta = "p_ctg.fasta"
    File a_ctg_fasta = "a_ctg.fasta"
  }
}
workflow falcon {
  input {
    Array[File] array_of_fasta
    #File ifile_config = "/scratch/cdunn/pbcromwell/General_config.json"
    File ifile_config = "General_config.json"
    Int nproc = 1
    Int max_nchunks = 256
  }
  Map[String, String] cfg = read_json(ifile_config)

  call task__0_rawreads__build {
    input:
      topdir         = "../..",
      config         = ifile_config,
      fasta_array    = array_of_fasta,
      pypeflow_nproc = "~{nproc}",
      pypeflow_mb    = "4000",
  }

  String dazzler_db_0 = read_string(task__0_rawreads__build.dazzler_db_txt)
  String dazzler_idx_0 = read_string(task__0_rawreads__build.dazzler_idx_txt)
  String dazzler_bps_0 = read_string(task__0_rawreads__build.dazzler_bps_txt)
  String dazzler_dust_anno_0 = read_string(task__0_rawreads__build.dazzler_dust_anno_txt)
  String dazzler_dust_data_0 = read_string(task__0_rawreads__build.dazzler_dust_data_txt)

  call task__0_rawreads__tan_split {
    input:
      config            = ifile_config,
      dazzler_db        = dazzler_db_0,
      dazzler_idx       = dazzler_idx_0,
      dazzler_bps       = dazzler_bps_0,
      #dazzler_dust_anno = dazzler_dust_anno_0,
      #dazzler_dust_data = dazzler_dust_data_0,
  }
  call task__0_rawreads__tan_scatter {
    input:
      uows = task__0_rawreads__tan_split.uows,
      max_nchunks = max_nchunks
  }
  scatter (chunk in task__0_rawreads__tan_scatter.chunks) {
    call task__0_rawreads__tan_apply {
      input:
        chunk             = chunk,
        config            = ifile_config,
        dazzler_db        = dazzler_db_0,
        dazzler_idx       = dazzler_idx_0,
        dazzler_bps       = dazzler_bps_0,
        #dazzler_dust_anno = dazzler_dust_anno_0,
        #dazzler_dust_data = dazzler_dust_data_0,
    }
  }

  call task__0_rawreads__tan_combine {
    input:
      config            = ifile_config,
      dazzler_db        = dazzler_db_0,
      dazzler_idx       = dazzler_idx_0,
      dazzler_bps       = dazzler_bps_0,
      tan_anno_array_of_fofn = task__0_rawreads__tan_apply.dazzler_tan_anno_fofn,
      tan_data_array_of_fofn = task__0_rawreads__tan_apply.dazzler_tan_data_fofn,
  }

  String dazzler_tan_anno = read_string(task__0_rawreads__tan_combine.dazzler_tan_anno_txt)
  String dazzler_tan_data = read_string(task__0_rawreads__tan_combine.dazzler_tan_data_txt)

  call task__0_rawreads__daligner_split {
    input:
      config            = ifile_config,
      dazzler_db        = dazzler_db_0,
      dazzler_idx       = dazzler_idx_0,
      dazzler_bps       = dazzler_bps_0,
      dazzler_dust_anno = dazzler_dust_anno_0,
      dazzler_dust_data = dazzler_dust_data_0,
      length_cutoff     = task__0_rawreads__build.length_cutoff
  }
  call task__0_rawreads__daligner_scatter {
    input:
      uows = task__0_rawreads__daligner_split.uows,
      max_nchunks = max_nchunks
  }
  scatter (chunk in task__0_rawreads__daligner_scatter.chunks) {
    call task__0_rawreads__daligner_apply {
      input:
        chunk             = chunk,
        config            = ifile_config,
        dazzler_db        = dazzler_db_0,
        dazzler_idx       = dazzler_idx_0,
        dazzler_bps       = dazzler_bps_0,
        dazzler_dust_anno = dazzler_dust_anno_0,
        dazzler_dust_data = dazzler_dust_data_0,
        dazzler_tan_anno = dazzler_tan_anno,
        dazzler_tan_data = dazzler_tan_data,
    }
  }
  call task__0_rawreads__daligner_las_merge {
    input:
      las_array_of_fofns = task__0_rawreads__daligner_apply.result,
      config             = ifile_config,
      nproc              = nproc,
  }
  scatter (chunk in task__0_rawreads__daligner_las_merge.las_array) {
    call task__0_rawreads__cns_apply {
      input:
        las           = chunk,
        config        = ifile_config,
        dazzler_db    = dazzler_db_0,
        dazzler_idx   = dazzler_idx_0,
        dazzler_bps   = dazzler_bps_0,
        length_cutoff = task__0_rawreads__build.length_cutoff,
    }
  }
  call task__0_rawreads__report {
    input:
      config        = ifile_config,
      fasta_array   = task__0_rawreads__cns_apply.fasta,
      dazzler_db    = dazzler_db_0,
      dazzler_idx   = dazzler_idx_0,
      dazzler_bps   = dazzler_bps_0,
      length_cutoff = task__0_rawreads__build.length_cutoff,
  }
  call task__1_preads_ovl__build {
    input:
      topdir        = "../..",
      config        = ifile_config,
      fasta_array   = task__0_rawreads__cns_apply.fasta,
      pypeflow_nproc= "~{nproc}",
      pypeflow_mb   = "4000",
  }

  String dazzler_db_1 = read_string(task__1_preads_ovl__build.dazzler_db_txt)
  String dazzler_idx_1 = read_string(task__1_preads_ovl__build.dazzler_idx_txt)
  String dazzler_bps_1 = read_string(task__1_preads_ovl__build.dazzler_bps_txt)
  String dazzler_dust_anno_1 = read_string(task__1_preads_ovl__build.dazzler_dust_anno_txt)
  String dazzler_dust_data_1 = read_string(task__1_preads_ovl__build.dazzler_dust_data_txt)

  call task__1_preads_ovl__daligner_split {
    input:
      config=ifile_config,
      dazzler_db=dazzler_db_1,
      dazzler_idx=dazzler_idx_1,
      dazzler_bps=dazzler_bps_1,
      dazzler_dust_anno=dazzler_dust_anno_1,
      dazzler_dust_data=dazzler_dust_data_1,
      length_cutoff=task__1_preads_ovl__build.length_cutoff
  }
  call task__1_preads_ovl__daligner_scatter {
    input:
      uows = task__1_preads_ovl__daligner_split.uows,
      max_nchunks = max_nchunks
  }
  scatter (chunk in task__1_preads_ovl__daligner_scatter.chunks) {
    call task__1_preads_ovl__daligner_apply {
      input:
        chunk = chunk,
        config = ifile_config,
        dazzler_db = dazzler_db_1,
        dazzler_idx= dazzler_idx_1,
        dazzler_bps= dazzler_bps_1,
        dazzler_dust_anno=dazzler_dust_anno_1,
        dazzler_dust_data=dazzler_dust_data_1,
    }
  }
  call task__1_preads_ovl__daligner_las_merge {
    input:
      las_array_of_fofns = task__1_preads_ovl__daligner_apply.result,
      config             = ifile_config,
      nproc              = nproc,
  }
  call task__1_preads_ovl__db2falcon {
    input:
      dazzler_db = dazzler_db_1,
      dazzler_idx= dazzler_idx_1,
      dazzler_bps= dazzler_bps_1,
  }
  call task__2_asm_falcon {
    input:
      dazzler_db          = dazzler_db_1,
      dazzler_idx         = dazzler_idx_1,
      dazzler_bps         = dazzler_bps_1,
      las_fofn            = task__1_preads_ovl__daligner_las_merge.las_fofn,
      preads4falcon_fasta = task__1_preads_ovl__db2falcon.preads4falcon_fasta,
      overlap_filtering_setting = cfg["overlap_filtering_setting"],
      length_cutoff_pr          = cfg["length_cutoff_pr"],
      fc_ovlp_to_graph_option   = cfg["fc_ovlp_to_graph_option"],
      nproc               = nproc,
  }
  output {
    File ofile_pre_assembly_report = task__0_rawreads__report.stats
    File ofile_preads4falcon_fasta = task__1_preads_ovl__db2falcon.preads4falcon_fasta
    File ofile_p_ctg_fasta         = task__2_asm_falcon.p_ctg_fasta
    File ofile_a_ctg_fasta         = task__2_asm_falcon.a_ctg_fasta
    File report_preassembly        = task__0_rawreads__report.report
  }
}
