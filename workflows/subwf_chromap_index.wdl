version 1.0

import "../tasks/task_chromap_index.wdl" as chromap_index
import "../tasks/task_check_inputs.wdl" as task_check_inputs

workflow generate_chromap_index {
    meta {
        version: 'v1'
        author: 'Your Name (your.email@example.com)'
        description: 'IGVF Single Cell pipeline: Sub-workflow to create chromap index'
    }

    input {
        File genome_fasta
        String output_dir
    }

    if ( (sub(genome_fasta, "^gs:\/\/", "") == sub(genome_fasta, "", "")) ){
        call task_check_inputs.check_inputs as genome_check {
            input:
                path = genome_fasta
        }
    }

    call chromap_index.generate_chromap_index as chromap {
        input:
            genome_fasta = select_first([genome_check.output_file, genome_fasta]),
            output_dir = output_dir
    }

    output {
        File chromap_index_tar = chromap.atac_chromap_index
    }
}