version 1.1

task multiqc {
  input {
    File? user_config
    File? logo_png
    String? additional_commandline_args
    Array[File] input_logs
    String? report_prefix
    String? input_logs_folder
    Boolean? recursive
  }

  command <<< >>>

  output {
    File html_report = "placeholder.txt"
    Array[File]+ data_files = ["placeholder.txt"]
  }

  runtime {
    dx_app: object {
      type: "app",
      id: "app-G81jk9Q0Q1BQ8vqV3kv2pg8z",
      name: "multiqc/2.1.5"
    }
  }
}
