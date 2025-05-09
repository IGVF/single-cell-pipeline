version 1.0

# Import the tasks called by the pipeline
import "../tasks/task_submit_outputs.wdl" as task_submit

# WDL workflow for submitting single-cell processed files to IGVF portal

workflow submit_proccessed_files {

    input {
        String analysis_accession
        File igvf_credentials
    }

    call task_submit.submit as submit{
        input:
            analysis_accession = analysis_accession,
            igvf_credentials = igvf_credentials
    }

    output{    
        Array[File] submit_logs = submit.output_files
    }
}