version 1.0

# Import the tasks called by the pipeline
import "../tasks/task_seqspec_extract.wdl" as task_seqspec_extract
import "../tasks/task_chromap.wdl" as task_align_chromap
import "../tasks/task_chromap_bam.wdl" as task_align_chromap_bam
import "../tasks/task_log_atac.wdl" as task_log_atac


workflow wf_atac {
    meta {
        version: 'v1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) and Sai Ma @ Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Sub-workflow to process the ATAC portion of SHARE-seq libraries.'
    }

    input {
        # ATAC sub-workflow inputs
        String chemistry
        String prefix = "igvf_output_atac"
        String? subpool = "none"

        File? barcode_conversion_dict # For 10X multiome
        File reference_fasta

        # Align-specific inputs
        Array[File] read1
        Array[File] read2
        Array[File] fastq_barcode
        Array[File] seqspecs
        Array[File] barcode_inclusion_list
        File reference_index_tar_gz
        String? read_format
        # Runtime parameters
        Int? align_cpus
        Float? align_disk_factor = 8.0
        Float? align_memory_factor = 0.15
        String? align_docker_image
        
        Int? align_bam_cpus
        Float? align_bam_disk_factor = 8.0
        Float? align_bam_memory_factor = 0.15

        # Merge-specific inputs
        # Runtime parameters
        Int? merge_cpus
        Float? merge_disk_factor = 8.0
        Float? merge_memory_factor = 0.15
        String? merge_docker_image

        # QC-specific inputs
        File? raw_bam
        File? raw_bam_index
        File? filtered_bam
        File? filtered_bam_index

        # RNA seqspec extract runtime parameters
        Int? seqspec_extract_cpus
        Float? seqspec_extract_disk_factor
        Float? seqspec_extract_memory_factor
        String? seqspec_extract_docker_image
    }
    # run only if seqspec array is not empty
    if (length(seqspecs) > 0 && !defined(read_format)) {
        #should implement check if length of seqspecs == length of read1 == length of read2
        scatter ( idx in range(length(seqspecs)) ) {
            call task_seqspec_extract.seqspec_extract as seqspec_extract {
                input:
                    seqspec = seqspecs[idx],
                    fastq_R1 = basename(read1[idx]),
                    fastq_R2 = basename(read2[idx]),
                    fastq_barcode = basename(fastq_barcode[idx]),
                    onlists = barcode_inclusion_list,
                    modality = "atac",
                    tool_format = "chromap",
                    onlist_format = "product",
                    chemistry = chemistry,
                    cpus = seqspec_extract_cpus,
                    disk_factor = seqspec_extract_disk_factor,
                    memory_factor = seqspec_extract_memory_factor,
                    docker_image = seqspec_extract_docker_image
            }
        }
    }
    
    #Assuming this whitelist is applicable to all fastqs for kb task
    Array[File] barcode_inclusion_list_ = select_first([barcode_inclusion_list, seqspec_extract.onlist])
    
    String index_string_ = select_first([read_format, seqspec_extract.index_string ])
        
    call task_align_chromap.chromap_generate_fragments as generate_fragments {
        input:
            fastq_R1 = read1,
            fastq_R2 = read2,
            fastq_barcode = fastq_barcode,
            reference_fasta = reference_fasta,
            reference_index_tar_gz = reference_index_tar_gz,
            subpool = subpool,
            prefix = prefix,
            barcode_inclusion_list = barcode_inclusion_list_[0],
            barcode_conversion_dict = barcode_conversion_dict,
            disk_factor = align_disk_factor,
            memory_factor = align_memory_factor,
            cpus = align_cpus,
            docker_image = align_docker_image,
            read_format = index_string_
    }

    call task_align_chromap_bam.chromap_generate_bam as generate_bam {
        input:
            fastq_R1 = read1,
            fastq_R2 = read2,
            fastq_barcode = fastq_barcode,
            reference_fasta = reference_fasta,
            reference_index_tar_gz = reference_index_tar_gz,
            subpool = subpool,
            prefix = prefix,
            barcode_inclusion_list = barcode_inclusion_list_[0],
            barcode_conversion_dict = barcode_conversion_dict,
            disk_factor = align_disk_factor,
            memory_factor = align_memory_factor,
            cpus = align_cpus,
            docker_image = align_docker_image,
            read_format = index_string_
    }

    call task_log_atac.log_atac as log_atac {
        input:
            alignment_log = generate_fragments.atac_alignment_log,
            barcode_log = generate_fragments.atac_barcode_summary,
            prefix = prefix
    }

    output {
        # Align
        File atac_fragments_alignment_stats = generate_fragments.atac_alignment_log
        File atac_fragments = generate_fragments.atac_fragments
        File atac_fragments_index = generate_fragments.atac_fragments_index
        File atac_chromap_barcode_summary = generate_fragments.atac_barcode_summary
        File? atac_qc_metrics = log_atac.atac_statistics_json

        File atac_chromap_bam = generate_bam.atac_bam
        File? atac_chromap_bam_summary = generate_bam.atac_bam_summary
    }
}