version 1.0

# Import the sub-workflow for preprocessing the fastqs.
import "tasks/task_check_inputs.wdl" as check_inputs
import "workflows/subwf_atac.wdl" as subwf_atac
import "workflows/subwf_rna.wdl" as subwf_rna
import "tasks/10x_create_barcode_mapping.wdl" as tenx_barcode_map

# WDL workflow for processing single-cell multiomic datasets

workflow single_cell_pipeline {

    input {
        # Common inputs
        Boolean create_onlist_mapping = false
        String prefix # Analysis set
        File? igvf_credentials
        String? subpool = "none" # To address
        File genome_tsv

        # ATAC-specific inputs
        Array[File] atac_read1
        Array[File] atac_read2
        Array[File] fastq_barcode
        File atac_barcode_inclusion_list
        File? chromap_genome_index_tar_gz
        File? genome_fasta
        String atac_read_format

        # RNA-specific inputs
        Array[File] rna_read1
        Array[File] rna_read2
        Array[File] fastq_barcode_rna = []
        File rna_barcode_inclusion_list
        String kb_mode = "nac"
        String rna_read_format
        File? kb_genome_index_tar_gz
        File? rna_replacement_list #assume this is gs link or empty
    }

    Map[String, File] annotations = read_map(genome_tsv)
    File genome_fasta_ = select_first([genome_fasta, annotations["fasta"]])
    File idx_tar_rna_ = select_first([kb_genome_index_tar_gz, annotations["kb_nac_idx_tar"]])
    File idx_tar_atac_ = select_first([chromap_genome_index_tar_gz, annotations["chromap_idx_tar"]])

    Boolean process_atac = if length(atac_read1)>0 then true else false
    Boolean process_rna = if length(rna_read1)>0 then true else false


    if(process_atac){
        if ( atac_read1[0] != "" ) {

            #ATAC Read1
            if ( (sub(atac_read1[0], "^gs:\/\/", "") == sub(atac_read1[0], "", "")) ){
                scatter(file in atac_read1){
                    call check_inputs.check_inputs as check_read1_atac{
                        input:
                            path = file,
                            igvf_credentials = igvf_credentials
                    }
                }
            }
            
            #ATAC Read2
            if ( (sub(atac_read2[0], "^gs:\/\/", "") == sub(atac_read2[0], "", "")) ){
                scatter(file in atac_read2){
                    call check_inputs.check_inputs as check_read2_atac{
                        input:
                            path = file,
                            igvf_credentials = igvf_credentials
                    }
                }
            }

            #ATAC barcode
            if ( (sub(fastq_barcode[0], "^gs:\/\/", "") == sub(fastq_barcode[0], "", "")) ){
                scatter(file in fastq_barcode){
                    call check_inputs.check_inputs as check_fastq_barcode{
                        input:
                            path = file,
                            igvf_credentials = igvf_credentials
                    }
                }
            }
            
            if (sub(genome_fasta_, "^gs:\/\/", "") == sub(genome_fasta_, "", "")){
                call check_inputs.check_inputs as check_genome_fasta{
                        input:
                            path = genome_fasta_,
                            igvf_credentials = igvf_credentials
                    }
            }

            if (sub(idx_tar_atac_, "^gs:\/\/", "") == sub(idx_tar_atac_, "", "")){
                call check_inputs.check_inputs as check_genome_index{
                    input:
                        path = idx_tar_atac_,
                        igvf_credentials = igvf_credentials
                }
            }
        }
    }
    
    Array[File] read1_atac_ = select_first([ check_read1_atac.output_file, atac_read1 ])
    Array[File] read2_atac_ = select_first([ check_read2_atac.output_file, atac_read2 ])
    Array[File] fastq_barcode_ = select_first([ check_fastq_barcode.output_file, fastq_barcode ])
    
    if(process_rna){
        if ( rna_read1[0] != "" ) {

            #RNA Read1
            if ( (sub(rna_read1[0], "^gs:\/\/", "") == sub(rna_read1[0], "", "")) ){
                scatter(file in rna_read1){
                    call check_inputs.check_inputs as check_read1_rna{
                        input:
                            path = file,
                            igvf_credentials = igvf_credentials
                    }
                }
            }

            #RNA Read2
            if ( (sub(rna_read2[0], "^gs:\/\/", "") == sub(rna_read2[0], "", "")) ){
                scatter(file in rna_read2){
                    call check_inputs.check_inputs as check_read2_rna{
                        input:
                            path = file,
                            igvf_credentials = igvf_credentials
                    }
                }
            }

            #RNA barcode
            if (length(fastq_barcode_rna) > 0){
                if ( (sub(fastq_barcode_rna[0], "^gs:\/\/", "") == sub(fastq_barcode_rna[0], "", "")) ){
                    scatter(file in fastq_barcode_rna){
                        call check_inputs.check_inputs as check_fastq_barcode_rna{
                            input:
                                path = file,
                                igvf_credentials = igvf_credentials
                        }
                    }
                }
            }

            #RNA barcode replacement list
            if (defined(rna_replacement_list)){
                File nonoptional_replacement_list = select_first([rna_replacement_list])
                if ( (sub(nonoptional_replacement_list, "^gs:\/\/", "") == sub(nonoptional_replacement_list, "", "")) ){
                    call check_inputs.check_inputs as check_rna_replacement_list{
                        input:
                            path = nonoptional_replacement_list,
                            igvf_credentials = igvf_credentials
                    }
                }
            }
            if (sub(idx_tar_rna_, "^gs:\/\/", "") == sub(idx_tar_rna_, "", "")){
                call check_inputs.check_inputs as check_transcriptome_index{
                    input:
                        path = idx_tar_rna_,
                        igvf_credentials = igvf_credentials
                }
            }
        }
    }
    
    
    Array[File] read1_rna_ = select_first([ check_read1_rna.output_file, rna_read1 ])
    Array[File] read2_rna_ = select_first([ check_read2_rna.output_file, rna_read2 ])
    Array[File] fastq_barcode_rna_ = select_first([ check_fastq_barcode_rna.output_file, fastq_barcode_rna ])
            
    if ( create_onlist_mapping && process_atac && process_rna){
        call tenx_barcode_map.mapping_tenx_barcodes as barcode_mapping{
            input:
                atac_barcode_inclusion_list = atac_barcode_inclusion_list,
                rna_barcode_inclusion_list = rna_barcode_inclusion_list
        }
    }
    
    if ( process_rna ) {
        if ( rna_read1[0] != "" ) {
            call subwf_rna.wf_rna as rna{
                input:
                    read1 = read1_rna_,
                    read2 = read2_rna_,
                    read_barcode = fastq_barcode_rna_,
                    barcode_inclusion_list = rna_barcode_inclusion_list,
                    kb_mode = kb_mode,
                    kb_index_tar_gz = select_first([check_transcriptome_index.output_file, idx_tar_rna_]),
                    prefix = prefix,
                    subpool = subpool,
                    read_format = rna_read_format,
                    replacement_list = if defined(check_rna_replacement_list.output_file) then check_rna_replacement_list.output_file else rna_replacement_list
            }
        }
    }

    if ( process_atac ) {
        if ( atac_read1[0] != "" ) {
            call subwf_atac.wf_atac as atac{
                input:
                    read1 = read1_atac_,
                    read2 = read2_atac_,
                    fastq_barcode = fastq_barcode_,
                    reference_fasta = select_first([check_genome_fasta.output_file, genome_fasta_]),
                    subpool = subpool,
                    barcode_inclusion_list = atac_barcode_inclusion_list,
                    reference_index_tar_gz = select_first([check_genome_index.output_file, idx_tar_atac_]),
                    prefix = prefix,
                    read_format = atac_read_format,
                    barcode_conversion_dict = barcode_mapping.tenx_barcode_conversion_dict,
            }
        }
    }

    output{
        # RNA outputs
        File? rna_kb_h5ad = rna.rna_kb_h5ad
        File? rna_kb_output_folder_tar_gz = rna.rna_kb_output_folder_tar_gz
        File? rna_kb_run_info_json = rna.rna_kb_run_info_json
        File? rna_kb_library_qc_metrics_json = rna.rna_kb_library_qc_metrics_json
        File? rna_kb_parameters_json = rna.rna_kb_parameters_json
    

        # ATAC outputs
        File? atac_bam = atac.atac_chromap_bam
        File? atac_bam_index = atac.atac_chromap_bam_index
        File? atac_bam_summary_stats = atac.atac_chromap_bam_summary
        File? atac_fragments = atac.atac_fragments
        File? atac_fragments_index = atac.atac_fragments_index
        File? atac_fragments_log = atac.atac_fragments_alignment_stats
        File? atac_fragments_metrics = atac.atac_fragments_qc_metrics
        File? atac_fragments_alignment_stats = atac.atac_fragments_alignment_stats
        File? atac_fragments_barcode_summary = atac.atac_fragments_barcode_summary

        String? subwf_rna_version = rna.subwf_rna_version
        String? subwf_atac_version = atac.subwf_atac_version
        String main_pipeline_version = "v1"
    }
}