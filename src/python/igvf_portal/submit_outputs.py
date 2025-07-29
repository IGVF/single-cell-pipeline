from igvf_utils.connection import Connection
import gcs_functions as gs
import sys
import re
import argparse
import json

def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1", "True", "T", "TRUE")

def main():
    parser = argparse.ArgumentParser(description="Submit single-cell pipeline outputs to IGVF.")
    
    # Adding arguments
    parser.add_argument("--atac_bam", required=False, help="Path to the ATAC BAM file.")
    parser.add_argument("--atac_bam_index", required=False, help="Path to the ATAC BAM index file.")
    parser.add_argument("--atac_fragment", required=False, help="Path to the ATAC fragment file.")
    parser.add_argument("--atac_fragment_index", required=False, help="Path to the ATAC fragment index file.")
    parser.add_argument("--atac_bam_summary_stats", required=False, help="Path to the ATAC QC metrics file.")
    parser.add_argument("--atac_fragment_alignment_stats", required=False, help="Path to the ATAC QC metrics file.")
    parser.add_argument("--atac_fragment_barcode_summary", required=False, help="Path to the ATAC QC metrics file.")
    parser.add_argument("--atac_fragment_metrics", required=False, help="Path to the ATAC QC metrics file.")
    parser.add_argument("--rna_h5ad", required=False, help="Path to the RNA H5AD file.")
    parser.add_argument("--rna_kb_tar", required=False, help="Path to the RNA KB tar file.")
    parser.add_argument("--rna_qc_kb_info", required=False, help="Path to the RNA QC metrics file.")
    parser.add_argument("--rna_qc_kb_parameters", required=False, help="Path to the RNA QC metrics file.")
    parser.add_argument("--rna_qc_inspect", required=False, help="Path to the RNA QC metrics file.")
    parser.add_argument("--lab", required=False, help="Lab name")
    parser.add_argument("--lab_key", required=False, help="Lab alias key")
    parser.add_argument("--award", required=False, help="Lab award")
    parser.add_argument("--analysis_set_acc", required=False)
    parser.add_argument("--genome", required=False, help="Genome assembly")
    parser.add_argument("--controlled_access", required=False, help="Controlled access flag")

    parser.add_argument(
        "--atac_mm_list",
        required=False,
        help="List of IGVF accessions.",
        type=lambda s: [re.search(r'IGVF\w+', part.strip()).group() for part in s.strip('[]').strip('[]').split(',') if re.search(r'IGVF\w+', part.strip())]
    )

    parser.add_argument(
        "--atac_r1_acc",
        required=False,
        help="List of IGVF accessions.",
        type=lambda s: [re.search(r'IGVF\w+', part.strip()).group() for part in s.strip('[]').strip('[]').split(',') if re.search(r'IGVF\w+', part.strip())]
    )

    parser.add_argument(
        "--atac_r2_acc",
        required=False,
        help="List of IGVF accessions.",
        type=lambda s: [re.search(r'IGVF\w+', part.strip()).group() for part in s.strip('[]').strip('[]').split(',') if re.search(r'IGVF\w+', part.strip())]
    )

    parser.add_argument(
        "--atac_bc_acc",
        required=False,
        help="List of IGVF accessions.",
        type=lambda s: [re.search(r'IGVF\w+', part.strip()).group() for part in s.strip('[]').strip('[]').split(',') if re.search(r'IGVF\w+', part.strip())]
    )

    parser.add_argument(
        "--rna_mm_list",
        required=False,
        help="List of IGVF accessions.",
        type=lambda s: [re.search(r'IGVF\w+', part.strip()).group() for part in s.strip('[]').strip('[]').split(',') if re.search(r'IGVF\w+', part.strip())]
    )

    parser.add_argument(
        "--rna_bc_acc",
        required=False,
        help="List of IGVF accessions.",
        type=lambda s: [re.search(r'IGVF\w+', part.strip()).group() for part in s.strip('[]').strip('[]').split(',') if re.search(r'IGVF\w+', part.strip())]
    )

    parser.add_argument(
        "--rna_r1_acc",
        required=False,
        help="List of IGVF accessions.",
        type=lambda s: [re.search(r'IGVF\w+', part.strip()).group() for part in s.strip('[]').strip('[]').split(',') if re.search(r'IGVF\w+', part.strip())]
    )

    parser.add_argument(
        "--rna_r2_acc",
        required=False,
        help="List of IGVF accessions.",
        type=lambda s: [re.search(r'IGVF\w+', part.strip()).group() for part in s.strip('[]').strip('[]').split(',') if re.search(r'IGVF\w+', part.strip())]
    )
    
    parser.add_argument(
        "--atac_seqspec_acc",
        required=False,
        help="List of configuration file URLs (comma-separated). Extracts IGVF accessions.",
        type=lambda s: [re.search(r'IGVF\w+', part.strip()).group() for part in s.strip('[]').strip('[]').split(',') if re.search(r'IGVF\w+', part.strip())]
    )

    parser.add_argument(
        "--rna_seqspec_acc",
        required=False,
        help="List of configuration file URLs (comma-separated). Extracts IGVF accessions.",
        type=lambda s: [re.search(r'IGVF\w+', part.strip()).group() for part in s.strip('[]').strip('[]').split(',') if re.search(r'IGVF\w+', part.strip())]
    )

    args = parser.parse_args()
    
    print("Lab:", args.lab)
    print("Lab key:", args.lab_key)
    print("Award:", args.award)
    print("Analysis Set Accession:", args.analysis_set_acc)
    
    args.controlled_access = str2bool(args.controlled_access)

    conn=Connection("prod")

    atac_reference_files = {
        "GRCm39": ["IGVFFI5593VLWB", "IGVFFI9282QLXO"],
        "GRCh38": ["IGVFFI7969JLFC", "IGVFFI0653VCGH"],
    }

    rna_reference_files = {
        "GRCm39": ["IGVFFI5078MNED"],
        "GRCh38": ["IGVFFI9561BASO"],
    }

    #ATAC BAM
    if args.atac_bam:
        print("ATAC r1 Accession:", args.atac_r1_acc)
        print("ATAC r2 Accession:", args.atac_r2_acc)
        print("ATAC bc Accession:", args.atac_bc_acc)
        print("ATAC seqspec Accession:", args.atac_seqspec_acc)
        print("ATAC MM Accession:", args.atac_mm_list)
        print("ATAC BAM:", args.atac_bam)
        payload = {}
        payload["submitted_file_name"] = args.atac_bam
        payload["md5sum"] = gs.get_md5sum(args.atac_bam.split("/")[2], "/".join(args.atac_bam.split("/")[3:]))
        payload["aliases"] = [args.lab_key + args.analysis_set_acc + "_bam"]
        payload["lab"] = args.lab 
        payload["award"] = args.award
        payload["file_format"] = "bam"
        payload["file_set"] = args.analysis_set_acc
        payload["content_type"] = "alignments"
        payload["derived_from"] = args.atac_r1_acc + args.atac_r2_acc + args.atac_bc_acc + args.atac_seqspec_acc
        payload["controlled_access"] = args.controlled_access
        payload["redacted"] = False
        payload["filtered"] = False
        payload["assembly"] = args.genome
        payload["reference_files"] = atac_reference_files[args.genome]
        payload["analysis_step_version"] = "/analysis-step-versions/39c0498d-91f6-42de-8896-2fab1403f032/"
        payload[Connection.PROFILE_KEY] = "alignment_file"
        _schema_property = conn.get_profile_from_payload(payload).properties
        print(payload)    
        
        conn.post(payload, upload_file=True, truncate_long_strings_in_payload_log=True)

    #ATAC BAM Index
    if args.atac_bam_index:
        print("ATAC BAM Index:", args.atac_bam_index)
        payload = {}
        payload["submitted_file_name"] = args.atac_bam_index
        payload["md5sum"] = gs.get_md5sum(args.atac_bam_index.split("/")[2], "/".join(args.atac_bam_index.split("/")[3:]))
        payload["aliases"] = [args.lab_key + args.analysis_set_acc + "_bam_index"]
        payload["lab"] = args.lab 
        payload["award"] = args.award
        payload["file_format"] = "bai"
        payload["file_set"] = args.analysis_set_acc
        payload["content_type"] = "index"
        payload["derived_from"] = [args.lab_key + args.analysis_set_acc + "_bam"]
        payload["controlled_access"] = args.controlled_access
        payload["analysis_step_version"] = "/analysis-step-versions/39c0498d-91f6-42de-8896-2fab1403f032/"
        payload[Connection.PROFILE_KEY] = "index_file"
        _schema_property = conn.get_profile_from_payload(payload).properties
        print(payload)    
        
        conn.post(payload, upload_file=True, truncate_long_strings_in_payload_log=True)

    #ATAC Fragment
    if args.atac_fragment:
        print("ATAC Fragments:", args.atac_fragment)
        payload = {}
        payload["submitted_file_name"] = args.atac_fragment
        payload["md5sum"] = gs.get_md5sum(args.atac_fragment.split("/")[2], "/".join(args.atac_fragment.split("/")[3:]))
        payload["aliases"] = [args.lab_key + args.analysis_set_acc + "_fragments"]
        payload["lab"] = args.lab 
        payload["award"] = args.award
        payload["file_format"] = "bed"
        payload["file_format_type"] = "bed3+"
        payload["file_set"] = args.analysis_set_acc
        payload["content_type"] = "fragments"
        payload["derived_from"] = args.atac_r1_acc + args.atac_r2_acc + args.atac_bc_acc + args.atac_seqspec_acc
        payload["controlled_access"] = False
        payload["filtered"] = False
        payload["assembly"] = args.genome
        payload["file_format_specifications"] = ["buenrostro-bernstein:igvf-single-cell-pipeline-fragment-file-specification"]
        payload["analysis_step_version"] = "/analysis-step-versions/39c0498d-91f6-42de-8896-2fab1403f032/"
        payload[Connection.PROFILE_KEY] = "tabular_file"
        _schema_property = conn.get_profile_from_payload(payload).properties
        print(payload)

        conn.post(payload, upload_file=True, truncate_long_strings_in_payload_log=True)

    #ATAC Fragment Index
    if args.atac_fragment_index:
        print("ATAC Fragment Index:", args.atac_fragment_index)
        payload = {}
        payload["submitted_file_name"] = args.atac_fragment_index
        payload["md5sum"] = gs.get_md5sum(args.atac_fragment_index.split("/")[2], "/".join(args.atac_fragment_index.split("/")[3:]))
        payload["aliases"] = [args.lab_key + args.analysis_set_acc + "_fragment_index"]
        payload["lab"] = args.lab 
        payload["award"] = args.award
        payload["file_format"] = "tbi"
        payload["file_set"] = args.analysis_set_acc
        payload["content_type"] = "index"
        payload["derived_from"] = [args.lab_key + args.analysis_set_acc + "_fragments"]
        payload["controlled_access"] = False
        payload["analysis_step_version"] = "/analysis-step-versions/39c0498d-91f6-42de-8896-2fab1403f032/"
        payload[Connection.PROFILE_KEY] = "index_file"
        _schema_property = conn.get_profile_from_payload(payload).properties
        print(payload) 
        
        conn.post(payload, upload_file=True, truncate_long_strings_in_payload_log=True)

    #RNA H5AD
    if args.rna_h5ad:
        args.rna_bc_acc = []
        for i in args.rna_r1_acc:
            args.rna_bc_acc.append(conn.get(conn.get(i)["aliases"][0].rsplit("_",1)[0] + "_barcode" )["accession"])
        print("RNA H5AD:", args.rna_h5ad)
        payload = {}
        payload["submitted_file_name"] = args.rna_h5ad
        payload["md5sum"] = gs.get_md5sum(args.rna_h5ad.split("/")[2], "/".join(args.rna_h5ad.split("/")[3:]))
        payload["aliases"] = [args.lab_key + args.analysis_set_acc + "_matrix_file"]
        payload["lab"] = args.lab 
        payload["award"] = args.award
        payload["file_format"] = "h5ad"
        payload["file_set"] = args.analysis_set_acc
        payload["content_type"] = "sparse gene count matrix"
        payload["principal_dimension"] = "cell"
        payload["secondary_dimensions"] = ["gene"]
        payload["filtered"] = False
        payload["derived_from"] = args.rna_r1_acc + args.rna_r2_acc + args.rna_seqspec_acc + args.rna_bc_acc 
        payload["reference_files"] = rna_reference_files[args.genome]
        payload["analysis_step_version"] = "/analysis-step-versions/9c457b9f-fc6d-4cf1-b249-218827e9b449/"
        payload["file_format_specifications"] = ["buenrostro-bernstein:igvf-sc-pipeline-matrix-h5-specification"]
        payload[Connection.PROFILE_KEY] = "matrix_file"
        _schema_property = conn.get_profile_from_payload(payload).properties
        print(payload)    

        conn.post(payload, upload_file=True, truncate_long_strings_in_payload_log=True)

    #RNA KB Tar
    if args.rna_kb_tar:
        print("RNA KB Tar:", args.rna_kb_tar)
        payload = {}
        payload["submitted_file_name"] = args.rna_kb_tar
        payload["md5sum"] = gs.get_md5sum(args.rna_kb_tar.split("/")[2], "/".join(args.rna_kb_tar.split("/")[3:]))
        payload["aliases"] = [args.lab_key + args.analysis_set_acc + "_rna_kb_output_folder_tar_gz"]
        payload["lab"] = args.lab 
        payload["award"] = args.award
        payload["file_format"] = "tar"
        payload["file_set"] = args.analysis_set_acc
        payload["content_type"] = "kallisto single cell RNAseq output"
        payload["derived_from"] = args.rna_r1_acc + args.rna_r2_acc + args.rna_seqspec_acc + args.rna_bc_acc 
        payload["filtered"] = False
        payload["reference_files"] = rna_reference_files[args.genome]
        payload["principal_dimension"] = "cell"
        payload["secondary_dimensions"] = ["gene"]
        payload["file_format_specifications"] = ["buenrostro-bernstein:igvf-sc-pipeline-matrix-tar-specification", "igvf:igvf-sc-pipeline-rna-tar-mtx-per-file-specification"]
        payload["analysis_step_version"] = "/analysis-step-versions/9c457b9f-fc6d-4cf1-b249-218827e9b449/"
        payload[Connection.PROFILE_KEY] = "matrix_file"
        _schema_property = conn.get_profile_from_payload(payload).properties
        print(payload)    
        
        conn.post(payload, upload_file=True, truncate_long_strings_in_payload_log=True)

    #ATAC QC BAM
    if args.atac_bam_summary_stats:
        payload = {}
        payload[Connection.PROFILE_KEY] = "single_cell_atac_seq_quality_metric"
        payload["quality_metric_of"] = [args.lab_key + args.analysis_set_acc + "_bam"]
        payload["lab"] = args.lab 
        payload["award"] = args.award
        payload["description"] =  "ATACseq chromap QC metric"
        payload["analysis_step_version"] = "/analysis-step-versions/39c0498d-91f6-42de-8896-2fab1403f032/"
        payload["atac_bam_summary_stats"] = {"path": args.atac_bam_summary_stats}
        payload["aliases"] = [args.lab_key + args.analysis_set_acc +"_single_cell_atac_seq_quality_metric_alignment_uniform-pipeline"]
        _schema_property = conn.get_profile_from_payload(payload).properties
        print(payload)    
        
        conn.post(payload, upload_file=False, truncate_long_strings_in_payload_log=True)


    #ATAC QC fragment
    if args.atac_fragment_alignment_stats:
        payload = {}
        payload[Connection.PROFILE_KEY] = "single_cell_atac_seq_quality_metric"
        payload["quality_metric_of"] = [args.lab_key + args.analysis_set_acc + "_fragments"]
        payload["description"] =  "ATACseq chromap QC metric"
        payload["lab"] = args.lab 
        payload["award"] = args.award
        payload["analysis_step_version"] = "/analysis-step-versions/39c0498d-91f6-42de-8896-2fab1403f032/"
        payload["atac_fragment_summary_stats"] = {"path": args.atac_fragment_barcode_summary}
        payload["atac_fragments_alignment_stats"] = {"path": args.atac_fragment_alignment_stats}
        payload['metadata_map'] = {
            'Number_of_reads': 'n_reads',
            'Number_of_mapped_reads': 'n_mapped_reads',
            'Number_of_uniquely_mapped_reads': 'n_uniquely_mapped_reads',
            'Number_of_reads_have_multi-mappings': 'n_reads_with_multi_mappings',
            'Number_of_candidates': 'n_candidates',
            'Number_of_mappings': 'n_mappings',
            'Number_of_uni-mappings': 'n_uni_mappings',
            'Number_of_multi-mappings': 'n_multi_mappings',
            # TODO: check if this map is right
            'Number_of_barcodes_in_whitelist': 'n_barcodes_on_onlist',
            'Number_of_corrected_barcodes': 'n_corrected_barcodes',
            'Number_of_output_mappings_(passed_filters)': 'n_output_mappings',
            'uni-mappings': 'uni_mappings',
            'multi-mappings': 'multi_mappings',
            'total': 'total',
            'percentage_duplicates': 'pct_duplicates'
        }
        # Read the JSON file specified in args.rna_qc_kb_info into memory
        with open(args.atac_fragment_metrics, 'r') as file:
            atac_fragment_metrics_data = json.load(file)
        
        # Update values in payload['metadata_map'] using values in rna_qc_kb_info_data
        for key, value in payload['metadata_map'].items():
            if key in atac_fragment_metrics_data:
                payload[value] = atac_fragment_metrics_data[key]
        
        payload["aliases"] = [args.lab_key + args.analysis_set_acc +"_single_cell_atac_seq_quality_metric_fragment_uniform-pipeline"]
        payload.pop("metadata_map", None)
        _schema_property = conn.get_profile_from_payload(payload).properties
        print(payload)    
        
        conn.post(payload, upload_file=False, truncate_long_strings_in_payload_log=True)


    #RNA QC
    if args.rna_qc_kb_info:
        payload = {}
        payload[Connection.PROFILE_KEY] = "single_cell_rna_seq_quality_metric"
        payload["quality_metric_of"] = [args.lab_key + args.analysis_set_acc + "_matrix_file", 
                                        args.lab_key + args.analysis_set_acc + "_rna_kb_output_folder_tar_gz"]
        payload["lab"] = args.lab 
        payload["award"] = args.award
        payload["analysis_step_version"] = "/analysis-step-versions/9c457b9f-fc6d-4cf1-b249-218827e9b449/"
        payload["rnaseq_kb_info"] = {"path": args.rna_qc_kb_parameters}
        payload["description"] =  "RNAseq Kallisto Bustools QC metric"
        payload['metadata_map'] = {
            'numRecords': 'n_records',   # starting here is inspect.json
            'numReads': 'n_reads',
            'numBarcodes': 'n_barcodes',
            'medianReadsPerBarcode': 'median_reads_per_barcode',
            'meanReadsPerBarcode': 'mean_reads_per_barcode',
            'numUMIs': 'total_umis',
            'numBarcodeUMIs': 'n_barcode_umis',
            'medianUMIsPerBarcode': 'median_umis_per_barcode',
            'meanUMIsPerBarcode': 'mean_umis_per_barcode',
            'gtRecords': 'gt_records',
            'numBarcodesOnOnlist': 'num_barcodes_on_onlist',
            'percentageBarcodesOnOnlist': 'percentage_barcodes_on_onlist',
            'numReadsOnOnlist': 'num_reads_on_onlist',
            'percentageReadsOnOnlist': 'percentage_reads_on_onlist',
            'n_targets': 'n_targets',    # starting here is run_info.json
            'n_bootstraps': 'n_bootstraps',
            'n_processed': 'n_processed',
            'n_pseudoaligned': 'n_pseudoaligned',
            'n_unique': 'n_unique',
            'p_pseudoaligned': 'p_pseudoaligned',
            'p_unique': 'p_unique',
            'index_version': 'index_version',
            'k-mer length': 'kmer_length'
        }
        
        # Read the JSON file specified in args.rna_qc_kb_info into memory
        with open(args.rna_qc_kb_info, 'r') as file:
            rna_qc_kb_info_data = json.load(file)
        
        # Update values in payload['metadata_map'] using values in rna_qc_kb_info_data
        for key, value in payload['metadata_map'].items():
            if key in rna_qc_kb_info_data:
                payload[value] = rna_qc_kb_info_data[key]

        # Read the JSON file specified in args.rna_qc_kb_info into memory
        with open(args.rna_qc_inspect, 'r') as file:
            rna_qc_inspect_data = json.load(file)
        
        # Update values in payload['metadata_map'] using values in rna_qc_kb_info_data
        for key, value in payload['metadata_map'].items():
            if key in rna_qc_inspect_data:
                payload[value] = rna_qc_inspect_data[key]

        payload["aliases"] = [args.lab_key + args.analysis_set_acc + "_single_cell_rna_seq_quality_metric_gene_count_uniform-pipeline"]
        payload.pop("metadata_map", None)
        _schema_property = conn.get_profile_from_payload(payload).properties
        print(payload)    
        
        conn.post(payload, upload_file=False, truncate_long_strings_in_payload_log=True)

if __name__ == "__main__":
    main()