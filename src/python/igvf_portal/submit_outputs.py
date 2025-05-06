from igvf_utils.connection import Connection
import sys
#import gcs_functions as gs
import argparse
import json

def main():
    parser = argparse.ArgumentParser(description="Submit single-cell pipeline outputs to IGVF.")
    
    # Adding arguments
    #parser.add_argument("--atac_bam", required=True, help="Path to the ATAC BAM file.")
    #parser.add_argument("--atac_bam_index", required=True, help="Path to the ATAC BAM index file.")
    #parser.add_argument("--atac_fragment", required=True, help="Path to the ATAC fragment file.")
    #parser.add_argument("--atac_fragment_index", required=True, help="Path to the ATAC fragment index file.")
    parser.add_argument("--atac_bam_summary_stats", required=True, help="Path to the ATAC QC metrics file.")
    parser.add_argument("--atac_fragment_alignment_stats", required=True, help="Path to the ATAC QC metrics file.")
    parser.add_argument("--atac_fragment_barcode_summary", required=True, help="Path to the ATAC QC metrics file.")
    parser.add_argument("--atac_fragment_metrics", required=True, help="Path to the ATAC QC metrics file.")
    #parser.add_argument("--atac_mm_list", required=True, help="List of ATAC metadata accessions.")
    #parser.add_argument("--rna_h5ad", required=True, help="Path to the RNA H5AD file.")
    #parser.add_argument("--rna_kb_tar", required=True, help="Path to the RNA KB tar file.")
    parser.add_argument("--rna_qc_kb_info", required=True, help="Path to the RNA QC metrics file.")
    parser.add_argument("--rna_qc_kb_parameters", required=True, help="Path to the RNA QC metrics file.")
    parser.add_argument("--rna_qc_inspect", required=True, help="Path to the RNA QC metrics file.")
    #parser.add_argument("--rna_mm_list", required=True, help="List of RNA metadata accessions.")
    parser.add_argument("--lab", required=True, help="Lab name")
    parser.add_argument("--lab_key", required=True, help="Lab alias key")
    parser.add_argument("--award", required=True, help="Lab award")
    parser.add_argument("--analysis_set_acc", required=True)
    #parser.add_argument("--atac_r1_acc", required=True)
    #parser.add_argument("--atac_r2_acc", required=True)
    #parser.add_argument("--atac_seqspec_acc", required=True)
    #parser.add_argument("--rna_r1_acc", required=True)
    #parser.add_argument("--rna_r2_acc", required=True)
    #parser.add_argument("--rna_seqspec_acc", required=True)


    # Parse the arguments
    args = parser.parse_args()
    
    # Example usage of parsed arguments
    #print("ATAC BAM:", args.atac_bam)
    #print("ATAC BAM Index:", args.atac_bam_index)
    #print("ATAC Fragment:", args.atac_fragment)
    #print("ATAC Fragment Index:", args.atac_fragment_index)
    print("ATAC QC:", args.atac_bam_summary_stats)
    #print("ATAC Metadata List:", args.atac_mm_list)
    #print("RNA H5AD:", args.rna_h5ad)
    #print("RNA KB Tar:", args.rna_kb_tar)
    print("RNA QC:", args.rna_qc_kb_info)
    print("RNA QC:", args.rna_qc_inspect)
    #print("RNA Metadata List:", args.rna_mm_list)
    print("Lab:", args.lab)
    print("Lab key:", args.lab_key)
    print("Award:", args.award)
    print("Analysis Set Accession:", args.analysis_set_acc)
    #print("ATAC r1 Accession:", args.atac_r1_acc)
    #print("ATAC r2 Accession:", args.atac_r2_acc)
    #print("ATAC seqspec Accession:", args.atac_seqspec_acc)
    #print("RNA r1 Accession:", args.rna_r1_acc)
    #print("RNA r2 Accession:", args.rna_r2_acc)
    #print("RNA seqspec Accession:", args.rna_seqspec_acc)
    
    # Add your logic here to process the inputs and interact with IGVF or GCS
    conn=Connection("prod")

    #ATAC QC BAM
    payload = {}
    payload[Connection.PROFILE_KEY] = "single_cell_atac_seq_quality_metric"
    payload["quality_metric_of"] = [args.lab_key + args.analysis_set_acc + "_bam"]
    payload["lab"] = args.lab 
    payload["award"] = args.award
    payload["analysis_step_version"] = "/analysis-step-versions/39c0498d-91f6-42de-8896-2fab1403f032/"
    payload["atac_bam_summary_stats"] = {"path": args.atac_bam_summary_stats}
    payload["aliases"] = [args.lab_key + args.analysis_set_acc +"_single_cell_atac_seq_quality_metric_alignment_uniform-pipeline"]
    _schema_property = conn.get_profile_from_payload(payload).properties
    print(payload)    
    stdout = conn.post(payload, upload_file=False, truncate_long_strings_in_payload_log=True)


    #ATAC QC fragment
    payload = {}
    payload[Connection.PROFILE_KEY] = "single_cell_atac_seq_quality_metric"
    payload["quality_metric_of"] = [args.lab_key + args.analysis_set_acc + "_fragments"]
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
    stdout = conn.post(payload, upload_file=False, truncate_long_strings_in_payload_log=True)


    #RNA QC
    payload = {}
    payload[Connection.PROFILE_KEY] = "single_cell_rna_seq_quality_metric"
    payload["quality_metric_of"] = [args.lab+":" + args.analysis_set_acc + "_matrix_file", 
                                    args.lab+":" + args.analysis_set_acc + "_rna_kb_output_folder_tar_gz"]
    payload["lab"] = args.lab 
    payload["award"] = args.award
    payload["analysis_step_version"] = "/analysis-step-versions/9c457b9f-fc6d-4cf1-b249-218827e9b449/"
    payload["rnaseq_kb_info"] = {"path": args.rna_qc_kb_parameters}
    payload['description'] =  'RNAseq Kallisto Bustools QC metric',
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
        'numBarcodesOnOnlist': 'numBarcodesOnOnlist',
        'percentageBarcodesOnOnlist': 'percentageBarcodesOnOnlist',
        'numReadsOnOnlist': 'numReadsOnOnlist',
        'percentageReadsOnOnlist': 'percentageReadsOnOnlist',
        'n_targets': 'n_targets',    # starting here is run_info.json
        'n_bootstraps': 'n_bootstraps',
        'n_processed': 'n_processed',
        'n_pseudoaligned': 'n_pseudoaligned',
        'n_unique': 'n_unique',
        'p_pseudoaligned': 'p_pseudoaligned',
        'p_unique': 'p_unique',
        'index_version': 'index_version',
        'k-mer length': 'k-mer length'
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

    payload["aliases"] = [args.lab_key + args.analysis_set_acc +"_single_cell_rna_seq_quality_metric_gene_count_uniform-pipeline"]
    payload.pop("metadata_map", None)
    _schema_property = conn.get_profile_from_payload(payload).properties
    print(payload)    
    stdout = conn.post(payload, upload_file=False, truncate_long_strings_in_payload_log=True)

if __name__ == "__main__":
    main()