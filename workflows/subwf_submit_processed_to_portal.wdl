version 1.0

# Import the tasks called by the pipeline
import "../tasks/task_submit.wdl" as task_submit

# WDL workflow for submitting single-cell processed files to IGVF portal

workflow submit_proccessed_files {

    input {
        String analysis_accession
        Array[String] atac_mm_accession_list
        Array[String] rna_mm_accession_list

        String atac_bam
        String atac_bam_index
        String atac_fragment
        String atac_fragment_index

        String rna_h5ad
        String rna_kb_tar

        File igvf_credentials

        File atac_qc
        File rna_qc
    }

    call task_submit.submit as submit{
        input:
            analysis_accession = analysis_accession,
            #atac_mm_accession_list= atac_mm_accession_list, 
            #rna_mm_accession_list= rna_mm_accession_list,
            #atac_bam= atac_bam,
            #atac_bam_index= atac_bam_index,
            #atac_fragment= atac_fragment,
            #atac_fragment_index= atac_fragment_index,
            #rna_h5ad= rna_h5ad,
            #rna_kb_tar= rna_kb_tar,
            igvf_credentials= igvf_credentials,
            atac_qc= atac_qc,
            rna_qc= rna_qc
    }

}