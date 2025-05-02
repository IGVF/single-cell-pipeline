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
        #Array[String] atac_mm_accession_list
        #Array[String] rna_mm_accession_list

        #String atac_bam
        #String atac_bam_index
        #String atac_fragment
        #String atac_fragment_index

        #String rna_h5ad
        #String rna_kb_tar

        File igvf_credentials

        File atac_qc
        File rna_qc

        String docker = "polumechanos/check_inputs:main"
        String lab_key = "buenrostro-bernstein:"
        String lab = "/labs/jason-buenrostro/"
        String award = "/awards/HG011986/"
    }

    command <<<
        
        set -e

        # Export IGVF credentials if the file exists
        if [[ -f "~{igvf_credentials}" ]]; then
            while IFS= read -r line; do
                    export "$line"
            done < "~{igvf_credentials}"
        fi

        #Google auth

        python submit.py \
        --atac_qc atac_qc \
        --rna_qc rna_qc \
        --lab lab \ 
        --lab_key lab_key \ 
        --award award \
        --analysis_set_acc analysis_accession

    >>>

    runtime {
        cpu : 1
        memory : "10 GB"
        docker : "~{docker}"
    }
}