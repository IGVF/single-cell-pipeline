version 1.0

# TASK
# submit-pipeline-params-to-portal
 
task submit {

    meta {
        version: 'v1'
        author: 'Siddarth Wekhande at Broad Institute of MIT and Harvard'
        description: 'Submit pipeline parameters to IGVF portal'
    } 
    
    input {
        String analysis_accession
        File igvf_credentials
   
        String? kb_strand = "forward"
        String? atac_barcode_inclusion_list
        String? atac_read_format = "bc:15:22,bc:53:60,bc:91:98,r1:0:-1,r2:0:-1"
        String chemistry = "shareseq"
        String? chromap_genome_index_tar_gz = "gs://fc-secure-de19fd29-2253-41cd-9751-1788cf7ad1a5/submissions/intermediates/d11fbeb3-94ff-4e6c-b1e3-91c2d7eee97d/generate_chromap_index/2cac0f34-dac0-4298-b195-d5af95a3a2f5/call-chromap/chromap_IGVFFI0653VCGH.tar.gz"
        String? create_onlist_mapping = "false"
        String? prefix
        String? rna_barcode_inclusion_list
        String? rna_read_format = "0,15,23,0,53,61,0,91,99:2,0,10:1,0,0"
        String? subpool
        Array[String]? atac_read1 = "[]"
        Array[String]? atac_read2 = "[]"
        Array[String]? rna_read1 = "[]"
        Array[String]? rna_read2 = "[]"
        Array[String]? fastq_barcode = "[]"
        Array[String]? fastq_barcode_rna = "[]"
        String? genome_tsv = "gs://broad-buenrostro-pipeline-genome-annotations/IGVF_human_v43/IGVF_human_v43_Homo_sapiens_genome_files_hg38_v43.tsv"
        String? genome_fasta = "gs://fc-secure-de19fd29-2253-41cd-9751-1788cf7ad1a5/submissions/intermediates/56c3edb0-832a-4ca5-98c4-06ab82bd930a/generate_chromap_index/6c2a2005-bf1c-45d6-851f-170920e8cfde/call-genome_check/cacheCopy/glob-aae8b15f635ae9fc31e845b03c8537e4/IGVFFI0653VCGH.fasta.gz"
        String? kb_genome_index_tar_gz = "gs://fc-secure-de19fd29-2253-41cd-9751-1788cf7ad1a5/submissions/intermediates/7378e57c-d2c2-47ac-980d-b301ba9077a1/wf_rna/ff6039ff-1e1a-4471-938b-9c309770a16f/call-kb/kb_IGVFFI0653VCGH_IGVFFI7217ZMJZv2.tar.gz"
        String? docker = "swekhande/sw-dockers:submit-outputs"

        Boolean dry_run = true
    }

     command <<<

     set -e

        # Export IGVF credentials if the file exists
        if [[ -f "~{igvf_credentials}" ]]; then
            while IFS= read -r line; do
                    export "$line"
            done < "~{igvf_credentials}"
        fi

    cat > config.json <<EOF
{
  "atac_barcode_inclusion_list": "~{atac_barcode_inclusion_list}",
  "atac_read1": ~{sep=',' atac_read1},
  "atac_read2": ~{sep=',' atac_read2},
  "atac_read_format": "~{atac_read_format}",
  "chromap_genome_index_tar_gz": "~{chromap_genome_index_tar_gz}",
  "create_onlist_mapping": ~{create_onlist_mapping},
  "fastq_barcode": ~{sep=',' fastq_barcode},
  "fastq_barcode_rna": ~{sep=',' fastq_barcode_rna},
  "genome_tsv": "~{genome_tsv}",
  "genome_fasta": "~{genome_fasta}",
  "kb_genome_index_tar_gz": "~{kb_genome_index_tar_gz}",
  "kb_mode": "nac",
  "prefix": "~{prefix}",
  "rna_barcode_inclusion_list": "~{rna_barcode_inclusion_list}",
  "rna_read1": ~{sep=',' rna_read1},
  "rna_read2": ~{sep=',' rna_read2},
  "rna_read_format": "~{rna_read_format}",
  "rna_replacement_list": null,
  "single_cell_pipeline.atac.align_bam_cpus": null,
  "single_cell_pipeline.atac.align_bam_disk_factor": 12,
  "single_cell_pipeline.atac.align_bam_memory_factor": 0.15,
  "single_cell_pipeline.atac.align_cpus": null,
  "single_cell_pipeline.atac.align_disk_factor": 12,
  "single_cell_pipeline.atac.align_docker_image": null,
  "single_cell_pipeline.atac.align_memory_factor": 0.15,
  "single_cell_pipeline.atac.filtered_bam": null,
  "single_cell_pipeline.atac.filtered_bam_index": null,
  "single_cell_pipeline.atac.merge_cpus": null,
  "single_cell_pipeline.atac.merge_disk_factor": 8,
  "single_cell_pipeline.atac.merge_docker_image": null,
  "single_cell_pipeline.atac.merge_memory_factor": 0.15,
  "single_cell_pipeline.atac.raw_bam": null,
  "single_cell_pipeline.atac.raw_bam_index": null,
  "single_cell_pipeline.barcode_mapping.cpus": 16,
  "single_cell_pipeline.barcode_mapping.disk_factor": 0.5,
  "single_cell_pipeline.barcode_mapping.docker_image": "debian:bullseye-slim",
  "single_cell_pipeline.barcode_mapping.memory_factor": 0.15,
  "single_cell_pipeline.check_fastq_barcode.cpus": 1,
  "single_cell_pipeline.check_fastq_barcode.disk_factor": 1,
  "single_cell_pipeline.check_fastq_barcode.docker_image": "docker.io/igvf/check-inputs:v1",
  "single_cell_pipeline.check_fastq_barcode.memory_factor": 1,
  "single_cell_pipeline.check_fastq_barcode_rna.cpus": 1,
  "single_cell_pipeline.check_fastq_barcode_rna.disk_factor": 1,
  "single_cell_pipeline.check_fastq_barcode_rna.docker_image": "docker.io/igvf/check-inputs:v1",
  "single_cell_pipeline.check_fastq_barcode_rna.memory_factor": 1,
  "single_cell_pipeline.check_genome_fasta.cpus": 1,
  "single_cell_pipeline.check_genome_fasta.disk_factor": 1,
  "single_cell_pipeline.check_genome_fasta.docker_image": "docker.io/igvf/check-inputs:v1",
  "single_cell_pipeline.check_genome_fasta.memory_factor": 1,
  "single_cell_pipeline.check_genome_index.cpus": 1,
  "single_cell_pipeline.check_genome_index.disk_factor": 1,
  "single_cell_pipeline.check_genome_index.docker_image": "docker.io/igvf/check-inputs:v1",
  "single_cell_pipeline.check_genome_index.memory_factor": 1,
  "single_cell_pipeline.check_read1_atac.cpus": 1,
  "single_cell_pipeline.check_read1_atac.disk_factor": 1,
  "single_cell_pipeline.check_read1_atac.docker_image": "docker.io/igvf/check-inputs:v1",
  "single_cell_pipeline.check_read1_atac.memory_factor": 1,
  "single_cell_pipeline.check_read1_rna.cpus": 1,
  "single_cell_pipeline.check_read1_rna.disk_factor": 1,
  "single_cell_pipeline.check_read1_rna.docker_image": "docker.io/igvf/check-inputs:v1",
  "single_cell_pipeline.check_read1_rna.memory_factor": 1,
  "single_cell_pipeline.check_read2_atac.cpus": 1,
  "single_cell_pipeline.check_read2_atac.disk_factor": 1,
  "single_cell_pipeline.check_read2_atac.docker_image": "docker.io/igvf/check-inputs:v1",
  "single_cell_pipeline.check_read2_atac.memory_factor": 1,
  "single_cell_pipeline.check_read2_rna.cpus": 1,
  "single_cell_pipeline.check_read2_rna.disk_factor": 1,
  "single_cell_pipeline.check_read2_rna.docker_image": "docker.io/igvf/check-inputs:v1",
  "single_cell_pipeline.check_read2_rna.memory_factor": 1,
  "single_cell_pipeline.check_rna_replacement_list.cpus": 1,
  "single_cell_pipeline.check_rna_replacement_list.disk_factor": 1,
  "single_cell_pipeline.check_rna_replacement_list.docker_image": "docker.io/igvf/check-inputs:v1",
  "single_cell_pipeline.check_rna_replacement_list.memory_factor": 1,
  "single_cell_pipeline.check_transcriptome_index.cpus": 1,
  "single_cell_pipeline.check_transcriptome_index.disk_factor": 1,
  "single_cell_pipeline.check_transcriptome_index.docker_image": "docker.io/igvf/check-inputs:v1",
  "single_cell_pipeline.check_transcriptome_index.memory_factor": 1,
  "single_cell_pipeline.rna.kb_cpus": null,
  "single_cell_pipeline.rna.kb_disk_factor": 4,
  "single_cell_pipeline.rna.kb_docker_image": null,
  "single_cell_pipeline.rna.kb_memory_factor": null,
  "single_cell_pipeline.rna.kb_strand": "forward",
  "subpool": "~{subpool}"
}
EOF

    #Replace all single quotes with double quotes in config.json
    sed -i 's/'\''/"/g' config.json

    cat config.json

    cat > pipeline_parameters.json << EOF
[
    {
        "aliases": [
            "buenrostro-bernstein:~{analysis_accession}_pipeline_config"
        ],
        "lab": "jason-buenrostro",
        "award": "HG011986",
        "attachment": {
            "path": "config.json"
        },
        "description": "Terra workflow configuration for the single-cell pipeline run",
        "document_type": "pipeline parameters"
    }
]
EOF

        cat pipeline_parameters.json

        iu_register -p document -i pipeline_parameters.json -m prod

        echo "Posted document."

    if [ "~{dry_run}" = false ]; then
        echo "Not a dry run. Proceeding with patch."
        printf "record_id\tpipeline_parameters\n" > patch.tsv
        printf "~{analysis_accession}\t[\"buenrostro-bernstein:~{analysis_accession}_pipeline_config\"]\n" >> patch.tsv
        cat patch.tsv
        iu_register -p analysis_set -i patch.tsv -m prod --patch
    else
        echo "Dry run enabled. Skipping patch."
    fi

  >>>

  output {
    File config_json = "config.json"
  }

  runtime {
    cpu: 1
    memory: "10G"
    docker : "~{docker}"
  }

}