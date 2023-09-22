# Dev workflow to test scaling characteristics of SMRT Link and Cromwell

version 1.0

task make_chunks {
  input {
    Int n_chunks
  }
  command {
    python3 <<EOF
    for i in range(${n_chunks}):
      ofn = "chunk.%d.txt" % i
      open(ofn, "wt").write(str(i))
    EOF
  }
  runtime {
    cpu: 1
    memory: "100MB"
  }
  output {
    # this breaks at 13K files (if not sooner)
    Array[File] chunks = glob("chunk.*.txt")
  }
}

task run_chunk {
  input {
    File chunk_in
  }
  command {
    echo "i_chunk,time" > chunk.csv
    echo "`cat ${chunk_in}`,`time`" >> chunk.csv
  }
  runtime {
    cpu: 1
    memory: "10MB"
  }
  output {
    File csv_out = "chunk.csv"
  }
}

task gather_csv {
  input {
    Array[File] chunks
  }
  command {
    python3 -m pbcoretools.tasks.gather \
      --log-level DEBUG \
      gathered.csv ${sep=" " chunks}
  }
  runtime {
    cpu: 1
    memory: "500MB"
  }
  output {
    File csv_out = "gathered.csv"
  }
}

workflow dev_scaling {
  input {
    File eid_ccs
    # see SL-6710 for inspiration
    Int n_chunks = 836

    Int nproc = 1
    String log_level = "INFO"
    # this is ignored completely
    Int max_nchunks = 1
    File? tmp_dir
  }

  call make_chunks {
    input:
      n_chunks = n_chunks
  }

  scatter (chunk_in in make_chunks.chunks) {
    call run_chunk {
      input:
        chunk_in = chunk_in
    }
  }

  call gather_csv {
    input:
      chunks = run_chunk.csv_out
  }

  output {
    File csv_out = gather_csv.csv_out
  }
}
