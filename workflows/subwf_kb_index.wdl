version 1.0

import "../tasks/task_kb_index.wdl" as task_kb

workflow wf_rna {
    meta {
        version: 'v1'
        author: 'Siddarth Wekhande (swekhand@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'IGVF Single Cell pipeline: Sub-workflow to create kb reference'
    }
    
    call task_kb.kb_index as kb
    
    output {
        File rna_index = kb.rna_index
    }
}