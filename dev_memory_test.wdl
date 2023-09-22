# Dev workflow for testing memory parameterization

version 1.0

import "tasks/memory.wdl" as memory

task get_task_mem {
  input {
    File dataset
    Int mem_gb_min = 4
    Int nproc = 1
  }
  # This is completely arbitrary - it could be a Python program determining
  # these numbers, or even another C++ tool.
  Int mem_gb_base = 2
  Int mem_gb_per_core = 2
  Int default_mem_gb = mem_gb_base + (mem_gb_per_core * nproc)
  command {
    python3 -c "print(str(max(~{default_mem_gb}, ~{mem_gb_min})))" > mem_gb.txt
  }
  runtime {
    cpu: 1
    # The main limitation of this task is it must have a constant upper limit
    # on memory consumption
    memory: "1GB"
  }
  output {
    File mem_gb_txt = "mem_gb.txt"
    # XXX jira: SL-5656
    #Int mem_gb = read_int("mem_gb.txt")
  }
}

task mock_cpp_tool {
  input {
    File dataset
    Int nproc
    Int total_mem_gb
    Int? bam_size
    Int? index_memory_gb
  }
  command <<<
    python3 <<EOF
    from pbcommand.models.report import Report, Attribute
    #from pbcommand.utils import get_dataset_metadata
    total_mem_gb = ~{total_mem_gb}
    nproc = ~{nproc}
    ds_file = "~{dataset}"
    rpt = Report("memory_test",
                 #dataset_uuids=(get_dataset_metadata(ds_file).uuid,),
                 attributes=[
                   Attribute("nproc", value=nproc),
                   Attribute("total_mem_gb", value=total_mem_gb),
                   Attribute("bam_size", value=~{bam_size}),
                   Attribute("index_memory_gb", value=~{index_memory_gb})])
    rpt.write_json("memory.report.json")
    EOF
  >>>
  runtime {
    cpu: nproc
    memory: "${total_mem_gb}GB"
  }
  output {
    File report = "memory.report.json"
  }
}

workflow dev_memory_test {
  input {
    File eid_ccs
    Int mock_tool_mem_gb_min = 1

    Int nproc = 1
    String log_level = "INFO"
    String? tmp_dir
    Int max_nchunks = 1
  }

  call get_task_mem {
    input:
      dataset = eid_ccs,
      mem_gb_min = mock_tool_mem_gb_min,
      nproc = nproc
  }

  Int total_mem_gb = read_int(get_task_mem.mem_gb_txt)
  call mock_cpp_tool {
    input:
      dataset = eid_ccs,
      nproc = nproc,
      total_mem_gb = total_mem_gb,
      bam_size = get_input_sizes.bam_size,
      index_memory_gb = get_input_sizes.index_memory_gb
  }

  call memory.get_input_sizes {
    input:
      bam_dataset = eid_ccs,
      log_level = log_level,
      nproc = nproc
  }

  output {
    File report_memory = mock_cpp_tool.report
  }
}
