version 1.0

# Import the tasks called by the pipeline
import "../tasks/task_kb_count.wdl" as task_kb

workflow wf_rna {
    meta {
        version: 'v1'
        author: 'Siddarth Wekhande (swekhand@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'IGVF Single Cell pipeline: Sub-workflow to process RNA libraries'
    }

    input {
    
        Array[File] read1
        Array[File] read2
        Array[File]? read_barcode
                
        String read_format
        File? replacement_list

        File kb_index_tar_gz
        
        File barcode_inclusion_list
        
        String? subpool
        String prefix
        
        # RNA kb runtime parameters
        String kb_strand
        String kb_mode
        Int? kb_cpus
        Float? kb_disk_factor
        Float? kb_memory_factor
        String? kb_docker_image
        
    }


    call task_kb.kb_count as kb{
        input:
            read1_fastqs = read1,
            read2_fastqs = read2,
            read_barcode_fastqs = read_barcode,
            replacement_list = replacement_list,
            strand = kb_strand,
            kb_index_tar_gz = kb_index_tar_gz,
            kb_mode = kb_mode,
            output_dir = prefix,
            barcode_inclusion_list = barcode_inclusion_list,
            read_format = read_format,
            subpool = subpool,
            cpus = kb_cpus,
            disk_factor = kb_disk_factor,
            memory_factor = kb_memory_factor,
            docker_image = kb_docker_image
    }

    output {
        # RNA kb outputs
        File rna_kb_h5ad = kb.rna_kb_h5ad
        File rna_kb_output_folder_tar_gz = kb.rna_kb_output_folder_tar_gz
        File rna_kb_run_info_json = kb.rna_kb_run_info_json
        File rna_kb_library_qc_metrics_json = kb.rna_kb_library_qc_metrics_json
        File rna_kb_parameters_json = kb.rna_kb_parameters_json
        String subwf_rna_version = "v1"
    }
}
