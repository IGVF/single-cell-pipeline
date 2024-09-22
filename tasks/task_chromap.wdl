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
        File? barcode_inclusion_list
        File? barcode_conversion_dict
        File reference_fasta

        Boolean? trim_adapters = true
        Boolean? remove_pcr_duplicates = true
        Boolean? remove_pcr_duplicates_at_cell_level = true
        Boolean? Tn5_shift = true
        Boolean? low_mem = true
        Boolean? bed_output = true
        Int? max_insert_size = 2000
        Int? quality_filter = 0
        

        Int? multimappers = 4 # As per ENCODE pipeline
        Int? bc_error_threshold = 1
        Float? bc_probability_threshold = 0.9
        #TODO: This should come from a previous task parsing the seqspec.
        String? read_format 

        String? subpool = "none"
        String genome_name # GRCh38, mm10
        String prefix = "test-sample"

        Int? cpus = 8
        Float? disk_factor = 1
        #TODO: With this setting it usually caps at 75%.
        Float? memory_factor = 0.15
        #TODO:We need to setup a docker registry.
        String? docker_image = "us.gcr.io/buenrostro-share-seq/task_chromap:dev"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fastq_R1, "G") + size(fastq_R2, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 24.0 + size(reference_fasta, "G") + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    # Define the output names
    String fragments = '${prefix}.atac.filter.fragments.${genome_name}.tsv'
    String barcode_log = "${prefix}.atac.align.k${multimappers}.${genome_name}.barcode.summary.csv"
    String alignment_log = "${prefix}.atac.align.k${multimappers}.${genome_name}.log.txt"

    String monitor_log = "atac_align_monitor.log"

    command <<<
        set -e

        bash $(which monitor_script.sh) 1>&2 &

        # Extracting index
        echo '------ Extracting indexing ------' 1>&2
        time tar -xzf ~{reference_index_tar_gz}

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
        time chromap \
                ~{true='--trim-adapters ' false='' trim_adapters} \
                ~{true='--remove-pcr-duplicates ' false='' remove_pcr_duplicates} \
                ~{true='--remove-pcr-duplicates-at-cell-level ' false='' remove_pcr_duplicates_at_cell_level} \
                ~{true='--Tn5-shift ' false='' Tn5_shift} \
                ~{true='--low-mem ' false='' low_mem} \
                ~{true='--BED ' false='' bed_output} \
                ~{"-l " + max_insert_size} \
                ~{"--bc-error-threshold " + bc_error_threshold} \
                ~{"--bc-probability-threshold " + bc_probability_threshold} \
                ~{"--read-format " + read_format} \
                ~{"--drop-repetitive-reads " + multimappers} \
                -x chromap_index/index \
                -r ~{reference_fasta} \
                ~{"-q " + quality_filter} \
                -t ~{cpus} \
                -1 ~{sep="," fastq_R1} \
                -2 ~{sep="," fastq_R2} \
                -b ~{sep="," fastq_barcode} \
                --barcode-whitelist barcode_inclusion_list.txt \
                ~{"--barcode-translate " + barcode_conversion_dict} \
                -o ~{fragments} \
                --summary ~{barcode_log} > ~{alignment_log} 2>&1
        
        if [[ '~{subpool}' != "none" ]]; then
            echo '------  Add subpool to barcode name ------' 1>&2
            awk -v OFS="\t" -v subpool=~{subpool} '{$4=$4"_"subpool; print $0}' ~{fragments} > temp
            mv temp ~{fragments}
            awk -v FS="," -v OFS="," -v subpool=~{subpool} 'NR==1{print $0;next}{$1=$1"_"subpool; print $0}' ~{barcode_log} > temp
            mv temp ~{barcode_log}
        fi
        #TODO: We might want to correct here for the +4/-4
        bgzip -c ~{fragments} > ~{fragments}.gz
        tabix --zero-based --preset bed ~{fragments}.gz

    >>>

    output {
        File? atac_fragments = "~{fragments}.gz"
        File? atac_fragments_index = "~{fragments}.gz.tbi"
        File? atac_align_barcode_statistics = barcode_log
        File? atac_alignment_log = alignment_log
    }


    runtime {
        cpu: cpus
        docker: "${docker_image}"
        singularity: "docker://${docker_image}"
        disks: "local-disk ${disk_gb} ${disk_type}"
        memory: "${mem_gb} GB"
    }

    parameter_meta {
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
        trim_adapters: {
            description: 'Trim adapters.',
            help: 'Boolean flag to indicate if adapters should be trimmed.',
            default: true
        }
        remove_pcr_duplicates: {
            description: 'Remove PCR duplicates.',
            help: 'Boolean flag to indicate if PCR duplicates should be removed.',
            default: true
        }
        remove_pcr_duplicates_at_cell_level: {
            description: 'Remove PCR duplicates at cell level.',
            help: 'Boolean flag to indicate if PCR duplicates should be removed at the cell level.',
            default: true
        }
        Tn5_shift: {
            description: 'Tn5 shift.',
            help: 'Boolean flag to indicate if Tn5 shift should be applied.',
            default: true
        }
        low_mem: {
            description: 'Low memory mode.',
            help: 'Boolean flag to indicate if low memory mode should be used.',
            default: true
        }
        bed_output: {
            description: 'Output in BED format.',
            help: 'Boolean flag to indicate if output should be in BED format.',
            default: true
        }
        max_insert_size: {
            description: 'Maximum insert size.',
            help: 'Maximum insert size for alignment.',
            default: 2000
        }
        quality_filter: {
            description: 'Quality filter.',
            help: 'Quality filter threshold for alignment.',
            default: 0
        }
        multimappers: {
            description: 'Number of multimappers allowed.',
            help: 'This is the integer that will be passed to the -k parameter of chromap.',
            example: 4
        }
        bc_error_threshold: {
            description: 'Barcode error threshold.',
            help: 'Error threshold for barcode matching.',
            default: 1
        }
        bc_probability_threshold: {
            description: 'Barcode probability threshold.',
            help: 'Probability threshold for barcode matching.',
            default: 0.9
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
        genome_name: {
            description: 'Reference genome name.',
            help: 'The name of the reference genome used by the aligner. This is appended to the output file name.',
            examples: ['GRCh38', 'mm10']
        }
        prefix: {
            description: 'Prefix for output files.',
            help: 'Prefix that will be used to name the output files.',
            example: 'my-experiment'
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
