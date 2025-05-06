version 1.0

# Import the tasks called by the pipeline
import "../tasks/task_submit_outputs.wdl" as task_submit

# WDL workflow for submitting single-cell processed files to IGVF portal

workflow submit_proccessed_files {

    input {
        String analysis_accession
        
        #Array[String] atac_mm_accession_list
        #Array[String] rna_mm_accession_list

        #String atac_bam
        #String atac_bam_index
        #String atac_fragment
        #String atac_fragment_index

        #String rna_h5ad
        #String rna_kb_tar

        File igvf_credentials

        File atac_bam_summary_stats
        File atac_fragment_alignment_stats
        File atac_fragment_barcode_summary
        File atac_fragment_metrics
        
        File rna_qc_kb_info
        File rna_qc_kb_parameters
        File rna_qc_inspect
    }

    call task_submit.submit as submit{
        input:
            analysis_accession = analysis_accession,
            igvf_credentials = igvf_credentials,
            atac_bam_summary_stats = atac_bam_summary_stats, 
            atac_fragment_alignment_stats = atac_fragment_alignment_stats,
            atac_fragment_barcode_summary = atac_fragment_barcode_summary,
            atac_fragment_metrics = atac_fragment_metrics,
            rna_qc_kb_info = rna_qc_kb_info,
            rna_qc_kb_parameters = rna_qc_kb_parameters,
            rna_qc_inspect = rna_qc_inspect
    }
}