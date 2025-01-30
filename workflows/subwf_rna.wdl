version 1.0

# Import the tasks called by the pipeline
import "../tasks/task_seqspec_extract.wdl" as task_seqspec_extract
import "../tasks/task_kb_count.wdl" as task_kb
import "../tasks/task_qc_rna.wdl" as task_qc_rna
import "../tasks/task_log_rna.wdl" as task_log_rna

workflow wf_rna {
    meta {
        version: 'v0.1'
        author: 'Siddarth Wekhande (swekhand@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'IGVF Single Cell pipeline: Sub-workflow to process RNA libraries'
    }

    input {
    
        Array[File] read1
        Array[File] read2
        Array[File]? read_barcode
                
        Array[File] seqspecs
        String? read_format
        String chemistry
        File? replacement_list

        File kb_index_tar_gz
        
        Array[File] barcode_inclusion_list
        
        String? subpool
        String genome_name # GRCh38, mm10
        String prefix
        
        # RNA kb runtime parameters
        String kb_strand
        String kb_mode
        Int? kb_cpus
        Float? kb_disk_factor
        Float? kb_memory_factor
        String? kb_docker_image

        
        # RNA seqspec extract runtime parameters
        Int? seqspec_extract_cpus
        Float? seqspec_extract_disk_factor
        Float? seqspec_extract_memory_factor
        String? seqspec_extract_docker_image
        
        # RNA QC runtime parameters
        Int? qc_rna_cpus
        Float? qc_rna_disk_factor
        Float? qc_rna_memory_factor
        String? qc_rna_docker_image  
    }
    
    #Assuming this whitelist is applicable to all fastqs for kb task
    if (length(seqspecs) > 0 && !defined(read_format)) {
        #should implement check if length of seqspecs == length of read1 == length of read2
        scatter ( idx in range(length(seqspecs)) ) {
            call task_seqspec_extract.seqspec_extract as seqspec_extract {
                input:
                    seqspec = seqspecs[idx],
                    fastq_R1 = basename(read1[idx]),
                    fastq_R2 = basename(read2[idx]),
                    onlists = barcode_inclusion_list,
                    modality = "rna",
                    tool_format = "kb",
                    chemistry = chemistry,
                    #onlist_format = if chemistry=="shareseq" || chemistry=="parse" then "multi" else "product",
                    onlist_format = "product", #temp fix until bustools bug is fixed
                    cpus = seqspec_extract_cpus,
                    disk_factor = seqspec_extract_disk_factor,
                    memory_factor = seqspec_extract_memory_factor,
                    docker_image = seqspec_extract_docker_image
            }
        }
    }
    Array[File] barcode_inclusion_list_ = select_first([barcode_inclusion_list, seqspec_extract.onlist])
    
    String read_format_ = select_first([read_format, seqspec_extract.index_string ])

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
            barcode_inclusion_list = barcode_inclusion_list_[0],
            read_format = read_format_,
            subpool = subpool,
            cpus = kb_cpus,
            disk_factor = kb_disk_factor,
            memory_factor = kb_memory_factor,
            docker_image = kb_docker_image
    }
    
    call task_qc_rna.qc_rna as qc_rna {
        input:
            counts_h5ad = kb.rna_kb_h5ad,
            genome_name = genome_name,
            kb_workflow = kb_mode,
            prefix = prefix,
            cpus = qc_rna_cpus,
            disk_factor = qc_rna_disk_factor,
            memory_factor = qc_rna_memory_factor,
            docker_image = qc_rna_docker_image
    }
    
    #need to add duplicate logs from qc_rna in this task
    call task_log_rna.log_rna as log_rna {
        input:
            alignment_json = kb.rna_kb_alignment_metrics_json,
            barcodes_json = kb.rna_kb_barcode_metrics_json,
            genome_name = genome_name, 
            prefix = prefix          
    }

    output {
        File rna_align_log = kb.rna_kb_alignment_metrics_json
        File rna_kb_output = kb.rna_kb_compressed_output_folder
        File rna_h5ad = kb.rna_kb_h5ad
        File rna_log = log_rna.rna_logfile
        File rna_barcode_metadata = qc_rna.rna_barcode_metadata
        File? rna_umi_barcode_rank_plot = qc_rna.rna_umi_barcode_rank_plot
        File? rna_gene_barcode_rank_plot = qc_rna.rna_gene_barcode_rank_plot
        File? rna_gene_umi_scatter_plot = qc_rna.rna_gene_umi_scatter_plot
    }
}
