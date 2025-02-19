version 1.0

# TASK
# atac-chromap

task atac_align_chromap {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard IGVF pipeline: align ATAC task using chromap'
    }

    input {
        # This task takes in input the preprocessed ATAC fastqs and align them to the genome.
        Array[File] fastq_R1
        Array[File] fastq_R2
        Array[File]? fastq_barcode
        File reference_index_tar_gz
        File barcode_inclusion_list
        File? barcode_conversion_dict
        File reference_fasta        
        String read_format
        String prefix

        String? subpool

        Int? cpus = 8
        Float? disk_factor = 1
        #TODO: With this setting it usually caps at 75%.
        Float? memory_factor = 0.15
        #TODO:We need to setup a docker registry.
        String? docker_image = "docker.io/igvf/chromap:dev"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fastq_R1, "G") + size(fastq_R2, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 24.0 + size(reference_fasta, "G") + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    String monitor_log = "atac_align_monitor.log.txt"

    String index_dir = basename(reference_index_tar_gz, ".tar.gz")

    command <<<
        set -e

        bash $(which monitor_script.sh) 1>&2 &

        # Extracting index
        echo '------ Extracting indexing ------' 1>&2
        tar xvzf ~{reference_index_tar_gz} --no-same-owner -C ./

        if [[ '~{barcode_inclusion_list}' == *.gz ]]; then
            echo '------ Decompressing the barcode inclusion list ------' 1>&2
            gunzip -c ~{barcode_inclusion_list} > barcode_inclusion_list.txt
        else
            echo '------ No decompression needed for the barcode inclusion list ------' 1>&2
            cat ~{barcode_inclusion_list} > barcode_inclusion_list.txt
        fi

        # [r1|r2|bc]:start:end:strand
        # --read-format bc:0:15,r1:16:-1
        # The start and end are inclusive and -1 means the end of the read. User may use multiple fields to specify non-consecutive segments, e.g. bc:0:15,bc:32:-1.
        # The strand is presented by '+' and '-' symbol, if '-' the barcode will be reverse-complemented after extraction
        echo '------ align chromap ------' 1>&2
        run_chromap align \
            --index_dir ~{index_dir} \
            --read_format ~{read_format} \
            --reference_fasta ~{reference_fasta} \
            --prefix ~{prefix} \
            ~{"--subpool " + subpool} \
            --threads ~{cpus} \
            --barcode_onlist barcode_inclusion_list.txt \
            ~{"--barcode_translate " + barcode_conversion_dict} \
            --read1 ~{sep="," fastq_R1} \
            --read2 ~{sep="," fastq_R2} \
            --read_barcode ~{sep="," fastq_barcode}


    >>>

    output {
        File atac_fragments = "~{prefix}.fragments.tsv.gz"
        File atac_fragments_index = "~{prefix}.fragments.tsv.gz.tbi"
        File atac_barcode_summary = "~{prefix}.barcode.summary.csv"
        File atac_alignment_log = "~{prefix}.log.txt"
    }


    runtime {
        cpu: cpus
        docker: "${docker_image}"
        disks: "local-disk ${disk_gb} ${disk_type}"
        memory: "${mem_gb} GB"
    }

    parameter_meta {
        prefix: {
            description: 'Prefix for the output files.',
            help: 'Prefix for the output files.',
            example: 'output.atac'
        }
        fastq_R1: {
            description: 'Read1 fastq.',
            help: 'Processed fastq for read1.',
            example: 'input.atac.R1.fq.gz'
        }
        fastq_R2: {
            description: 'Read2 fastq.',
            help: 'Processed fastq for read2.',
            example: 'input.atac.R2.fq.gz'
        }
        fastq_barcode: {
            description: 'Barcode fastq.',
            help: 'Processed fastq for barcodes.',
            example: 'input.atac.barcode.fq.gz'
        }
        reference_index_tar_gz: {
            description: 'Reference index tar.gz file.',
            help: 'Compressed tarball containing the reference index files.',
            example: 'reference_index.tar.gz'
        }
        barcode_inclusion_list: {
            description: 'Barcode inclusion list.',
            help: 'List of barcodes to include in the analysis.',
            example: 'barcode_inclusion_list.txt'
        }
        barcode_conversion_dict: {
            description: 'Barcode conversion dictionary.',
            help: 'Dictionary for converting barcodes.',
            example: 'barcode_conversion_dict.txt'
        }
        reference_fasta: {
            description: 'Reference fasta file.',
            help: 'Reference genome fasta file.',
            example: 'reference_genome.fa'
        }

        read_format: {
            description: 'Read format.',
            help: 'Format of the reads for alignment.',
            example: 'bc:0:15,r1:16:-1'
        }
        subpool: {
            description: 'Subpool identifier.',
            help: 'Identifier for the subpool.',
            default: 'none'
        }

        cpus: {
            description: 'Number of CPUs.',
            help: 'Set the number of CPUs used by the aligner.',
            default: 8
        }
        disk_factor: {
            description: 'Disk factor.',
            help: 'Multiplication factor to determine disk required for task align.',
            default: 1.0
        }
        memory_factor: {
            description: 'Memory factor.',
            help: 'Multiplication factor to determine memory required for task align.',
            default: 0.15
        }
        docker_image: {
            description: 'Docker image.',
            help: 'Docker image for the alignment step.',
            example: 'us.gcr.io/buenrostro-share-seq/task_chromap:dev'
        }
    }
}
