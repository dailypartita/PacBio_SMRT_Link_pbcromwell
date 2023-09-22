version 1.0

workflow test_workflow {
  input {
    File input_1
    File? input_2
    Int? nproc = 1

    Int sleep_time = 0
    Int exit_code = 0
  }

  call task1 {
    input:
      input_1 = input_1,
      sleep_time = sleep_time,
      exit_code = exit_code
  }

  if (defined(input_2)) {
    call task2 {
      input:
        input_2 = input_2
    }
  }

  output {
    File output_1 = task1.output_1
    File? output_2 = task2.output_2
  }
}

task task1 {
  input {
    File input_1
    Int sleep_time = 0
    Int exit_code = 0
  }
  command <<<
    sleep ~{sleep_time}
    echo "Hello, ${USER}!" > hello.txt
    cat ~{input_1} >> hello.txt
    exit ~{exit_code}
  >>>
  runtime {
    backend: "Local"
  }
  output {
    File output_1 = "hello.txt"
  }
}

task task2 {
  input {
    File? input_2
  }
  command {
    cat ${input_2} > goodbye.txt
  }
  runtime {
    backend: "Local"
  }
  output {
    File output_2 = "goodbye.txt"
  }
}
