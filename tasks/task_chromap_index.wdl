version 1.0

# TASK
# atac-chromap

task generate_chromap_index {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard IGVF pipeline: align ATAC task using chromap'
    }

    input {
        File genome_fasta
        String output_dir

        Int? cpus = 1
        Float? disk_factor = 1
        #TODO: With this setting it usually caps at 75%.
        Float? memory_factor = 0.15
        #TODO:We need to setup a docker registry.
        String? docker_image = "docker.io/igvf/chromap:main"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(genome_fasta, "G")
    # Determining memory size base on the size of the input files.
    Float mem_gb = 32.0 + size(genome_fasta, "G") + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    String monitor_log = "atac_align_monitor.log"

    command <<<
        set -e

        # Create index
        mkdir ~{output_dir}
        echo '------ indexing ------' 1>&2
        run_chromap index --genome_fasta ~{genome_fasta} --output_dir ~{output_dir}

    >>>

    output {
        File atac_chromap_index = "~{output_dir}.tar.gz"
    }


    runtime {
        cpu: cpus
        docker: "${docker_image}"
        singularity: "docker://${docker_image}"
        disks: "local-disk ${disk_gb} ${disk_type}"
        memory: "${mem_gb} GB"
    }

    parameter_meta {
        output_dir: {
                description: 'Output directory.',
                help: 'Output directory for the index.',
                default: "chromap_index"
            }
        genome_fasta: {
                description: 'Reference FASTA file.',
                help: 'Reference FASTA file for the alignment step.',
                example: ["hg38.fa"]
            }
        cpus: {
                description: 'Number of cpus.',
                help: 'Set the number of cpus used by bowtie2',
                default: 16
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
