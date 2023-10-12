version 1.0

# TASK
# rna-qc

task qc_rna {
    meta {
        version: 'v0.1'
        author: 'Siddarth Wekhande (swekhand@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'calculate QC metrics from RNA'
    }

    input {
        # This function takes in input counts_unfiltered directory from kb-counts in a gzipped format
        
        File mtx_tar
        Int? umi_cutoff = 100
        Int? gene_cutoff = 100
        String genome_name
        String? subpool
        String? prefix
        
        
        #Will these be used? Need to run tests to optimize
        Int? cpus = 2
        Float? disk_factor = 1.0
        Float? memory_factor = 0.5
          
        String? docker_image = "swekhande/shareseq-prod:igvf-qc-rna"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(mtx_tar, "G")

    # Determining memory size based on the size of the input files.
    Float mem_gb = 8.0 + memory_factor * input_file_size_gb

    # Determining disk size based on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type based on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    
    String barcode_metadata = "~{default="share-seq" prefix}.qc.rna.~{genome_name}.barcode.metadata.tsv"
    String duplicates_log = "~{default="share-seq" prefix}.qc.rna.~{genome_name}.duplicates.log.txt"
    String umi_barcode_rank_plot = "~{default="share-seq" prefix}.qc.rna.~{genome_name}.umi.barcode.rank.plot.png"
    String gene_barcode_rank_plot = "~{default="share-seq" prefix}.qc.rna.~{genome_name}.gene.barcode.rank.plot.png"
    String gene_umi_scatter_plot = "~{default="share-seq" prefix}.qc.rna.~{genome_name}.gene.umi.scatter.plot.png"
    String monitor_log = "monitor.log"

    command <<<
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &

        tar -xvzf ~{mtx_tar} .C ./

        cat cells_x_genes.mtx | awk -v OFS="\t" 'NR>3{count[$2]+=$3; tot[$2]+=1}END{for (bc in count){ print bc,count[bc],tot[bc]} }' | sort -k1,1n > barcode_count_statistics_dedup.raw.tsv
        
        echo -e "barcode\tunique_umi\tgenes_final\n" > barcode_count_statistics_dedup.tsv
        
        awk -v pkr=~{if defined(subpool) then "_~{subpool}" else ""} -v OFS="\t" 'FNR==NR{bc[NR]=$1}FNR!=NR{print bc[$1]pkr,$2,$3; delete bc[$1]}END{for(idx in bc){print bc[idx]pkr,0,0}}' <(cat cells_x_genes.barcodes.txt) barcode_count_statistics_dedup.raw.tsv | sort -k1,1 >> barcode_count_statistics_dedup.tsv

        # Extract barcode metadata (total counts, genes) from anndata file
        python3 $(which qc_rna_extract_metrics.py) adata.h5ad \
                                                 ~{subpool} \
                                                 tmp_metadata.tsv

        join -t $'\t' -e 0 -j1 <(cat tmp_metadata.tsv | (sed -u 1q;sort -k1,1)) barcode_count_statistics_dedup.tsv | \
        awk -v OFS="\t" 'NR==1{print $0,"FRIG"}NR>1{printf "%s\t%4.2f\n",$0,$9/$2}' > ~{barcode_metadata}


        awk 'NR>1{total+=$2; duplicate+=$2-$9; unique+=$9} END {print "total reads:", total; print "unique reads:", unique; print "duplicate reads:", duplicate; print "FRIG:",unique/total}' ~{barcode_metadata} > ~{duplicates_log}

        # Make QC plots - need to write script
        # Rscript $(which rna_qc_plots.R) ~{barcode_metadata} ~{umi_cutoff} ~{gene_cutoff} ~{umi_barcode_rank_plot} ~{gene_barcode_rank_plot} ~{gene_umi_scatter_plot}
    >>>

    output {
        File rna_barcode_metadata = "~{barcode_metadata}"
        File rna_duplicates_log = "~{duplicates_log}"
        File? rna_umi_barcode_rank_plot = "~{umi_barcode_rank_plot}"
        File? rna_gene_barcode_rank_plot = "~{gene_barcode_rank_plot}"
        File? rna_gene_umi_scatter_plot = "~{gene_umi_scatter_plot}"
    }

    runtime {
        cpu : cpus
        memory : "~{mem_gb} GB"
        disks: "local-disk ~{disk_gb} ~{disk_type}"
        docker : "${docker_image}"
    }

    parameter_meta {
        bam: {
                description: 'Alignment bam file',
                help: 'Aligned reads in bam format.',
                example: 'hg38.aligned.bam'
            }
        umi_cutoff: {
                description: 'UMI cutoff',
                help: 'Cutoff for number of UMIs required when making UMI barcode rank plot.',
                example: 10
            }
        gene_cutoff: {
                description: 'Gene cutoff',
                help: 'Cutoff for number of genes required when making gene barcode rank plot.',
                example: 10
            }
        subpool: {
                description: 'Experiment subpool',
                help: 'Id of the sample subpool. Can be used to distinguish the 10X lanes.',
                examples: ['SS-PKR-000']
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name genome reference used to align.',
                example: ['hg38', 'mm10', 'hg19', 'mm9']
            }
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files',
                example: 'MyExperiment'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: samtools',
                example: ['put link to gcr or dockerhub']
            }
    }
}
