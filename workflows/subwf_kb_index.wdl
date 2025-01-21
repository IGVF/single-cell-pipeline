version 1.0

import "../tasks/task_kb_index.wdl" as task_kb
import "../tasks/task_check_inputs.wdl" as task_check_inputs

workflow wf_rna {
    meta {
        version: 'v1'
        author: 'Siddarth Wekhande (swekhand@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'IGVF Single Cell pipeline: Sub-workflow to create kb reference'
    }

    input {
        File genome_fasta
        File gene_gtf
        String kb_mode
        String output_folder
    }

    call task_check_inputs.check_inputs as genome_check {
        input:
            path = genome_fasta
    }

    call task_check_inputs.check_inputs as gtf_check {
        input:
            path = gene_gtf
    }
    
    call task_kb.kb_index as kb{
        input:
            genome_fasta = genome_check.output_file,
            gene_gtf = gtf_check.output_file,
            kb_mode = kb_mode,
            output_folder = output_folder
    }

    output {
        File rna_index = kb.rna_index
    }
}