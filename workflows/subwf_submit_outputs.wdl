version 1.0

# Import the tasks called by the pipeline
import "../tasks/task_submit_outputs.wdl" as task_submit
import "../tasks/task_submit_params.wdl" as task_params

# WDL workflow for submitting single-cell processed files to IGVF portal

workflow submit_proccessed_files {

    input {
        String analysis_accession
        File igvf_credentials
        Boolean? only_parameters = false

    }

    if (only_parameters == false) {
        # Additional logic or calls can be placed here if needed
        call task_submit.submit as submit{
        input:
            analysis_accession = analysis_accession,
            igvf_credentials = igvf_credentials
        }
    }

    call task_params.submit as params_submit{
        input:
            analysis_accession = analysis_accession,
            igvf_credentials = igvf_credentials
        }
    
    output{    
        Array[File]? submit_logs = submit.output_files
        File pipeline_params = params_submit.config_json
    }
}