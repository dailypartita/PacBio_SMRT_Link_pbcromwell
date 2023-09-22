# Workflow to download large software bundles e.g. from S3

version 1.0

task download_s3 {
  input {
    String file_uri
    String tmp_dir
  }
  command <<<
    set -e
    down_dir="~{tmp_dir}"
    file_name=$(echo "~{file_uri}" | sed 's/.*\///;')
    ( mkdir -p $down_dir && cd $down_dir && \
      aws s3 cp \
      "~{file_uri}"  .\
      --no-sign-request )
    echo `readlink -f "${down_dir}/${file_name}"` > download.fofn
  >>>
  runtime {
    cpu: 1
    memory: "100MB"
  }
  output {
    File downloaded = read_string("download.fofn")
  }
}

task download_curl {
  input {
    String file_uri
    String tmp_dir
  }
  command <<<
    set -e
    down_dir="~{tmp_dir}"
    file_name=$(echo "~{file_uri}" | sed 's/.*\///;')
    curl_exe=$(ls /usr/bin/curl || echo -n "curl")
    ( mkdir -p $down_dir && cd $down_dir && \
      ${curl_exe} -fsS \
      --retry 3 \
      -o ".${file_name}" \
      "~{file_uri}" && \
      mv ".${file_name}" "${file_name}" )
    echo `readlink -f "${down_dir}/${file_name}"` > download.fofn
  >>>
  runtime {
    cpu: 1
    memory: "100MB"
  }
  output {
    File downloaded = read_string("download.fofn")
  }
}

task sl_update {
  input {
    File downloaded
    String file_name
    String  root_dir
    String  tmp_dir
  }
  command <<<
    set -e
    runfile="~{tmp_dir}/*.run"
    gzfile="~{tmp_dir}/~{file_name}"
    if [[ $gzfile =~ [-_]autoupdate([-_.][^/]*)?\.(tar\.gz|tgz|tar\.bz2|tar|run)$ ]] ; then
      ~{root_dir}/admin/bin/smrtupdater --autoupdate-stage1 "$gzfile"
    else
      for i in {1..5}
      do
        if [ -f $gzfile ]; then
          break
        fi
        sleep 10
      done
      ( cd  "~{tmp_dir}" && \
        rm -f $runfile  && \
        tar -xvf $gzfile && \
       ~{root_dir}/admin/bin/smrtupdater --extract-bundles-only  $runfile )
    fi
  >>>
  runtime {
    cpu: 1
    memory: "100MB"
  }
  output {
    File updated = "${tmp_dir}/${file_name}"
  }
}

workflow sl_download_update {
  input {
    String file_prefix_uri
    String file_name
    String sl_root_dir = "."
    Boolean is_aws = false
    Boolean is_noop = true

    Int nproc = 1
    Int max_nchunks = 0
    String log_level = "INFO"
    String tmp_dir = "/tmp"
  }

  String download_dir = "${sl_root_dir}/userdata/downloads/sl-download"
  String file_uri = "${file_prefix_uri}/${file_name}"

  if (is_aws) {
    call download_s3 {
      input:
        file_uri = file_uri,
        tmp_dir = download_dir
    }
  }
  if (!is_aws) {
    call download_curl {
      input:
        file_uri = file_uri,
        tmp_dir = download_dir
    }
 }
 if (!is_noop) {
   call sl_update {
      input:
        downloaded =  select_first([download_s3.downloaded,
                                    download_curl.downloaded]),
        file_name = file_name,
        root_dir = sl_root_dir,
        tmp_dir = download_dir
   }
 }
 output {
    File downloaded = select_first([download_s3.downloaded,
                                    download_curl.downloaded])
 }
}
