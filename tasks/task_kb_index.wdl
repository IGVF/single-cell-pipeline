version 1.0

# TASK
# rna-kb

    
task kb_index {

    meta {
        version: 'v0.1'
        author: 'Siddarth Wekhande (swekhand@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'create reference using kb'
    }
    
    input {
         # This task takes in input genome and gtf and creates the index based on kb workflow
            
        File genome_fasta
        File gene_gtf   
        String kb_mode #kb_mode can be either "nac" or "standard"
        String output_folder
        
        Int? cpus = 4
        Float? disk_factor = 1
        Float? memory_factor = 1
        
        #TODO:We need to setup a docker registry.
        String? docker_image = "polumechanos/igvf-kb:dev"
        
    }
    
    
    # Determine the size of the input
    Float input_file_size_gb = size(genome_fasta, "G") + size(gene_gtf, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 48 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(50.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"


    command <<<
    
        set -e

        bash $(which monitor_script.sh) 1>&2 &
             
        mkdir ~{output_folder}

        run_kallisto index ~{kb_mode} \
            --temp_dir ./kb_tmp_dir \
            --genome_fasta ~{genome_fasta} \
            --gtf ~{gene_gtf} \
            --output_dir ~{output_folder}
        
    >>>

    output {
        File rna_index = "~{output_folder}.tar.gz"
    }

    runtime {
        cpu: cpus
        docker: "${docker_image}"
        singularity: "docker://${docker_image}"
        disks: "local-disk ${disk_gb} ${disk_type}"
        memory: "${mem_gb} GB"
    }
    
    parameter_meta {   
        genome_fasta: {
            description: 'Genome reference',
            help: 'Genome reference in .fa.gz file',
            example: 'hg38.fa.gz'
        }

        gene_gtf: {
            description: 'Gene annotations',
            help: 'Gene annotations in gtf file format',
            example: 'gencode_v47.gtf'
        }

        kb_mode: {
            description: 'KB alignment mode',
            help: 'Which strategy to use when aligning with kb',
            example: '"nac" or "standard"'
        }
        
        cpus: {
            description: 'Number of cpus.',
            help: 'Set the number of cpus used'
        }
            
        disk_factor: {
            description: 'Multiplication factor to determine disk required for task align.',
            help: 'This factor will be multiplied to the size of FASTQs to determine required disk of instance (GCP/AWS) or job (HPCs).',
            default: 8.0
        }
            
        memory_factor: {
            description: 'Multiplication factor to determine memory required for task align.',
            help: 'This factor will be multiplied to the size of FASTQs to determine required memory of instance (GCP/AWS) or job (HPCs).',
            default: 0.15
        }
             
        docker_image: {
            description: 'Docker image.',
            help: 'Docker image for the alignment step.',
            example: ["us.gcr.io/buenrostro-share-seq/share_task_bowtie2"]
        }
        
    }
    
}