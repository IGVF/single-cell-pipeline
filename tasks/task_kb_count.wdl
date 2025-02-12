version 1.0

# TASK
# rna-kb

    
task kb_count {

    meta {
        version: 'v0.1'
        author: 'Siddarth Wekhande (swekhand@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'align RNA using kb'
    }
    
    input {
         # This task takes in input the raw RNA fastqs and their associated index string and whitelist, processes the barcodes accordingly, aligns them to the genome, and outputs the count matrices.
                
        Array[File] read1_fastqs #These filenames must EXACTLY match the ones specified in seqspec
        Array[File] read2_fastqs #These filenames must EXACTLY match the ones specified in seqspec
        String kb_mode # standard or nac
        String output_dir
        String strand
        String read_format
        File kb_index_tar_gz
        File barcode_inclusion_list 
        
        Array[File]? read_barcode_fastqs
        String? subpool
        File? replacement_list

        Int threads = 4
        
        #Will these be used? Need to run tests to optimize
        Int cpus = 4
        Float disk_factor = 1.0
        Float memory_factor = 0.15
        
        #TODO:We need to setup a docker registry.
        String docker_image = "polumechanos/igvf-kb:dev"
        
    }
    
    
    # Determine the size of the input
    Float input_file_size_gb = size(read1_fastqs, "G") + size(read2_fastqs, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 24.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    # Define the output names
    String index_dir = basename(kb_index_tar_gz, ".tar.gz")
    String kb_run_info_json = "${output_dir}/run_info.json"
    String library_qc_metrics_json = "${output_dir}/inspect.json"
    String kb_info_json = "${output_dir}/kb_info.json"

    command <<<
    
        set -e

        bash $(which monitor_script.sh) 1>&2 &
        
        #set up fastq order as l1r1, l1r2, l2r1, l2r2, etc.
        interleaved_files_string=$(paste -d' ' <(printf "%s\n" ~{sep=" " read_barcode_fastqs}) <(printf "%s\n" ~{sep=" " read1_fastqs}) <(printf "%s\n" ~{sep=" " read2_fastqs}) | tr -s ' ')
           
        mkdir ~{output_dir}
        tar xvzf ~{kb_index_tar_gz} --no-same-owner -C ./
        
        if [[ '~{barcode_inclusion_list}' == *.gz ]]; then
            echo '------ Decompressing the RNA barcode inclusion list ------' 1>&2
            gunzip -c ~{barcode_inclusion_list} > barcode_inclusion_list.txt
        else
            echo '------ No decompression needed for the RNA barcode inclusion list ------' 1>&2
            cat ~{barcode_inclusion_list} > barcode_inclusion_list.txt
        fi

        run_kallisto quantify ~{kb_mode} \
            --index_dir ~{index_dir} \
            --read_format ~{read_format} \
            --output_dir ~{output_dir} \
            --strand ~{strand} \
            ~{"--subpool " + subpool} \
            ~{"--replacement_list " + replacement_list} \
            --threads ~{threads} \
            --barcode_onlist barcode_inclusion_list.txt \
            $interleaved_files_string
    >>>

    output {
        File rna_kb_output_folder_tar_gz = "~{output_dir}.tar.gz"
        File rna_kb_h5ad = "~{output_dir}.h5ad"
        File rna_kb_run_info_json = kb_run_info_json
        File rna_kb_library_qc_metrics_json = library_qc_metrics_json
        File rna_kb_parameters_json = kb_info_json
    }

    runtime {
        cpu: cpus
        docker: "${docker_image}"
        disks: "local-disk ${disk_gb} ${disk_type}"
        memory: "${mem_gb} GB"
    }
    
    parameter_meta {
    
        read1_fastqs: {
            description: 'List of input read1 fastqs.',
            help: 'List of raw read1 fastqs that will be corrected and processed using kb',
            example: ['l0.r1.fq.gz', 'l1.r1.fq.gz']
        }
        
        read2_fastqs: {
            description: 'List of input read2 fastqs.',
            help: 'List of raw read2 fastqs that will be corrected and processed using kb',
            example: ['l0.r2.fq.gz', 'l1.r2.fq.gz']
        }

        kb_mode: {
            description: 'kb mode.',
            help: 'kb mode to use for the alignment. Can be standard or nac',
            example: ['standard', 'nac']
        }

        output_dir: {
            description: 'Output directory.',
            help: 'Directory where the output files will be stored',
            example: ['output']
        }

        strand: {
            description: 'Strand.',
            help: 'Strand to use for the alignment. Can be forward, reverse, or both',
            example: ['forward', 'reverse', 'unstranded']
        }

        read_format: {
            description: 'Read format.',
            help: 'Read format to use for the alignment. Follows the kb documentation'
        }

        kb_index_tar_gz: {
            description: 'kb index tar gz.',
            help: 'Tar gz file containing the kb index',
            example: ['index.tar.gz']
        }

        barcode_inclusion_list: {
            description: 'Barcode inclusion list.',
            help: 'File containing the barcodes to include in the analysis',
            example: ['barcodes.txt']
        }

        read_barcode_fastqs: {
            description: 'List of input read barcode fastqs.',
            help: 'List of raw read barcode fastqs that will be corrected and processed using kb'
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