version 1.0

# TASK
# submit-to-portal
 
task submit {

    meta {
        version: 'v1'
        author: 'Siddarth Wekhande at Broad Institute of MIT and Harvard'
        description: 'Submit processed files to IGVF portal'
    }

    input {
        String analysis_accession
        
        #Array[String]? atac_mm_accession_list
        #Array[String]? rna_mm_accession_list

        #String? atac_bam
        String? atac_bam_index
        #String? atac_fragment
        #String? atac_fragment_index

        #String? rna_h5ad
        #String? rna_kb_tar

        File igvf_credentials

        File? atac_bam_summary_stats
        File? atac_fragment_alignment_stats
        File? atac_fragment_barcode_summary
        File? atac_fragment_metrics
        
        File? rna_qc_kb_info
        File? rna_qc_kb_parameters
        File? rna_qc_inspect

        String? docker = "swekhande/sw-dockers:submit-outputs"
        String? lab_key = "buenrostro-bernstein:"
        String? lab = "/labs/jason-buenrostro/"
        String? award = "/awards/HG011986/"
    }

    command <<<
        
        set -e

        # Export IGVF credentials if the file exists
        if [[ -f "~{igvf_credentials}" ]]; then
            while IFS= read -r line; do
                    export "$line"
            done < "~{igvf_credentials}"
        fi

        echo $IGVF_SECRET_KEY

        #Google auth

        echo "python /usr/local/bin/submit_outputs.py  \
        --atac_bam_summary_stats ~{atac_bam_summary_stats} \
        --atac_fragment_alignment_stats ~{atac_fragment_alignment_stats} \
        --atac_fragment_barcode_summary ~{atac_fragment_barcode_summary} \
        --atac_fragment_metrics ~{atac_fragment_metrics} \
        --rna_qc_kb_info ~{rna_qc_kb_info} \
        --rna_qc_kb_parameters ~{rna_qc_kb_parameters} \
        --rna_qc_inspect ~{rna_qc_inspect} \
        --lab ~{lab} \
        --lab_key ~{lab_key} \
        --award ~{award} \
        --analysis_set_acc ~{analysis_accession}"

        python3 /usr/local/bin/submit_outputs.py  \
        ~{if defined(atac_bam_summary_stats) then "--atac_bam_summary_stats " + atac_bam_summary_stats else ""} \
        ~{if defined(atac_fragment_alignment_stats) then "--atac_fragment_alignment_stats " + atac_fragment_alignment_stats else ""} \
        ~{if defined(atac_fragment_barcode_summary) then "--atac_fragment_barcode_summary " + atac_fragment_barcode_summary else ""} \
        ~{if defined(atac_fragment_metrics) then "--atac_fragment_metrics " + atac_fragment_metrics else ""} \
        ~{if defined(rna_qc_kb_info) then "--rna_qc_kb_info " + rna_qc_kb_info else ""} \
        ~{if defined(rna_qc_kb_parameters) then "--rna_qc_kb_parameters " + rna_qc_kb_parameters else ""} \
        ~{if defined(rna_qc_inspect) then "--rna_qc_inspect " + rna_qc_inspect else ""} \
        ~{if defined(atac_bam_index) then "--atac_bam_index " + atac_bam_index else ""} \
        --lab "~{lab}" \
        --lab_key "~{lab_key}" \
        --award "~{award}" \
        --analysis_set_acc "~{analysis_accession}"

    >>>

    output {
       Array[File] output_files = glob("IU_Logs/*")
    }

    runtime {
        cpu : 1
        memory : "10 GB"
        docker : "~{docker}"
    }
}