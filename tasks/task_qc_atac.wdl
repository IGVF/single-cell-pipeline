version 1.0

# TASK
# qc-atac

task qc_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard IGVF pipeline: ATAC qc statistics task'
    }

    input {
        # This function takes in input the raw and filtered bams
        # and compute some alignment metrics along with the TSS
        # enrichment plot.
        File? fragments
        File? fragments_index
        File? barcode_summary
        File? chrom_sizes
        File? tss
        File? barcode_conversion_dict

        Int? fragment_min_snapatac_cutoff = 1 # This is the fragment cutoff for snapatac
        Int? tsse_cutoff = 0 # This is the TSS enrichment cutoff for snapatac
        File? gtf
        String? genome_name
        String? prefix
        String? subpool="none"

        # Runtime
        Int? cpus = 60
        Float? disk_factor = 10.0
        Float? memory_factor = 0.3
        String docker_image = "docker.io/polumechanos/qc-atac-atomic:igvf"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fragments, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 32.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(100.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    # pdf string needed as required input to Picard CollectInsertSizeMetrics
    String hist_log_pdf = '${default="share-seq" prefix}.atac.qc.hist.${genome_name}.log.pdf'
    String hist_log_png = '${default="share-seq" prefix}.atac.qc.hist.${genome_name}.log.png'
    String hist_log = '${default="share-seq" prefix}.atac.qc.hist.${genome_name}.log.txt'
    String tss_pileup_prefix = '${default="share-seq" prefix}.atac.qc.tss.pileup.${genome_name}.log'
    String tss_pileup_out = '${default="share-seq" prefix}.atac.qc.tss.pileup.${genome_name}.log.png'
    String final_snapatac2_barcode_metadata = '${default="share-seq" prefix}.atac.qc.${genome_name}.snapatac2.barcode.metadata.tsv'
    String final_chromap_barcode_metadata = '${default="share-seq" prefix}.atac.qc.${genome_name}.chromap.barcode.metadata.tsv'
    String fragment_barcode_rank_plot = "${default="share-seq" prefix}.atac.qc.${genome_name}.fragment.barcode.rank.plot.png"


    String monitor_log = "atac_qc_monitor.log"


    command<<<
        set -e

        # I am not writing to a file anymore because Google keeps track of it automatically.
        bash $(which monitor_script.sh) 1>&2 &

        if [[ '~{gtf}' == *.gz ]]; then
            cp ~{gtf} gtf.gz
        else
            echo '------ Compressing GTF ------' 1>&2
            gzip -c ~{gtf} > gtf.gz   
        fi

        ln -s ~{fragments} in.fragments.tsv.gz
        ln -s ~{fragments_index} in.fragments.tsv.gz.tbi

        if [ ~{if defined(barcode_conversion_dict) then "true" else "false"} == "true" ]; then
            echo '------ There is a conversion list ------' 1>&2
            if [ '~{subpool}' != "none" ]; then
                echo '------ There is a subpool ------' 1>&2
                awk -v subpool=~{subpool} -v OFS="\t" '{print $1"_"subpool,$2"_"subpool}' ~{barcode_conversion_dict} > temp_conversion
            else
                cp ~{barcode_conversion_dict} temp_conversion
            fi
            awk -v FS='[,|\t]' -v OFS=',' 'FNR==NR{map[$2]=$1; next}FNR==1{print $0}FNR>1 && map[$1] {print map[$1],$2,$3,$4,$5}' temp_conversion ~{barcode_summary} > temp_summary
        else
            cp ~{barcode_summary} temp_summary
        fi

        cut -f1 ~{chrom_sizes} > list-names-chromosomes
        grep -wFf list-names-chromosomes ~{tss} > filtered.tss.bed

        # TSS enrichment stats
        echo '------ START: Compute TSS enrichment bulk ------' 1>&2
        time python3 /usr/local/bin/compute_tss_enrichment_bulk.py \
            -e 2000 \
            -p ~{cpus} \
            --regions filtered.tss.bed \
            --prefix "~{prefix}.atac.qc.~{genome_name}" \
            in.fragments.tsv.gz

        echo '------ START: Compute TSS enrichment snapatac2 ------' 1>&2
        echo '------ Extend TSS ------' 1>&2
        awk -v OFS="\t" '{if($2-150<0){$2=0}else{$2=$2-150};$3=$3+150; print $0}' filtered.tss.bed > tss.extended.bed
        /usr/local/bin/bedClip -verbose=2 tss.extended.bed ~{chrom_sizes} tss.extended.clipped.bed 2> tss.bedClip.log.txt
        echo '------ Promoter ------' 1>&2
        awk -v OFS="\t" '{if($2-2000<0){$2=0}else{$2=$2-2000};$3=$3+2000; print $0}' filtered.tss.bed > promoter.bed
        /usr/local/bin/bedClip -verbose=2 promoter.bed ~{chrom_sizes} promoter.clipped.bed 2> promoter.bedClip.log.txt

        echo '------ Snapatac2 ------' 1>&2
        # h5ad: f"{prefix}.snapatac.h5ad"
        # f"{prefix}.cutoff{fragment_min_snapatac_cutoff}.tsse.png"
        # f"{prefix}.cutoff500.tsse.png"
        # f"{prefix}.cutoff{fragment_min_snapatac_cutoff}.insertsize.distribution.png"
        # f"{prefix}.barcode_stats.tsv"
        #
        time python3 /usr/local/bin/snapatac2-tss-enrichment.py in.fragments.tsv.gz \
                                                                gtf.gz \
                                                                ~{chrom_sizes} \
                                                                tss.extended.clipped.bed \
                                                                promoter.clipped.bed \
                                                                ~{fragment_min_snapatac_cutoff} \
                                                                "~{prefix}.atac.qc.~{genome_name}"

        ls

        echo '------ START: Generate metadata ------' 1>&2

        awk -v FS=',' -v OFS=" " 'NR==1{$1=$1;print $0,"unique","pct_dup","pct_unmapped";next}{$1=$1;if ($2-$3-$4-$5>0){print $0,($2-$3-$4-$5),$3/($2-$4-$5),($5+$4)/$2} else { print $0,0,0,0}}' temp_summary  | sed 's/ /\t/g' > ~{final_chromap_barcode_metadata}

        cat "~{prefix}.atac.qc.~{genome_name}.barcode_stats.tsv" | sed 's/ /\t/g' > ~{final_snapatac2_barcode_metadata}

        # Barcode rank plot
        echo '------ START: Generate barcode rank plot ------' 1>&2
        time Rscript /usr/local/bin/atac_qc_plots.R ~{final_chromap_barcode_metadata} ~{fragment_min_snapatac_cutoff} ~{fragment_barcode_rank_plot}
    >>>

    output {
        File? atac_qc_final_hist_png = "~{prefix}.cutoff~{fragment_min_snapatac_cutoff}.insertsize.distribution.png"

        Float atac_qc_bulk_tsse = read_float("${prefix}.atac.qc.${genome_name}.tss_score_bulk.txt")

        File atac_qc_snapatac2_barcode_metadata = "~{final_snapatac2_barcode_metadata}"
        File atac_qc_chromap_barcode_metadata = "~{final_chromap_barcode_metadata}"
        File atac_qc_tss_enrichment_plot = "${prefix}.atac.qc.${genome_name}.tss_enrichment_bulk.png"
        File atac_qc_tsse_fragments_plot = "~{prefix}.atac.qc.~{genome_name}.cutoff~{fragment_min_snapatac_cutoff}.tsse.png"
        File atac_qc_tsse_strict_fragments_plot = "~{prefix}.atac.qc.~{genome_name}.cutoff500.tsse.png"
        File atac_qc_snapatac_h5ad = "~{prefix}.atac.qc.~{genome_name}.snapatac.h5ad"
        File? atac_qc_barcode_rank_plot = "~{fragment_barcode_rank_plot}"
    }

    runtime {
        cpu: cpus
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker: "${docker_image}"
        memory: "${mem_gb} GB"
    }

    parameter_meta {
        tss: {
                description: 'TSS bed file',
                help: 'List of TSS in bed format used for the enrichment plot.',
                example: 'refseq.tss.bed'
            }
        fragment_min_snapatac_cutoff: {
                description: 'Fragment cutoff',
                help: 'Cutoff for number of fragments required when making fragment barcode rank plot.',
                example: 10
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.',
                examples: ['hg38', 'mm10', 'both']
            }
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2',
                examples: '4'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz',
                example: ['put link to gcr or dockerhub']
            }
    }
}