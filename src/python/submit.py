from igvf_utils.connection import Connection
import sys
import gcs_functions as gs
import argparse

def main():
    parser = argparse.ArgumentParser(description="Submit single-cell pipeline outputs to IGVF.")
    
    # Adding arguments
    parser.add_argument("--atac_bam", required=True, help="Path to the ATAC BAM file.")
    parser.add_argument("--atac_bam_index", required=True, help="Path to the ATAC BAM index file.")
    parser.add_argument("--atac_fragment", required=True, help="Path to the ATAC fragment file.")
    parser.add_argument("--atac_fragment_index", required=True, help="Path to the ATAC fragment index file.")
    parser.add_argument("--atac_qc", required=True, help="Path to the ATAC QC metrics file.")
    parser.add_argument("--atac_mm_list", required=True, help="List of ATAC metadata accessions.")
    parser.add_argument("--rna_h5ad", required=True, help="Path to the RNA H5AD file.")
    parser.add_argument("--rna_kb_tar", required=True, help="Path to the RNA KB tar file.")
    parser.add_argument("--rna_qc", required=True, help="Path to the RNA QC metrics file.")
    parser.add_argument("--rna_mm_list", required=True, help="List of RNA metadata accessions.")
    parser.add_argument("--lab", required=True, help="Lab name")
    parser.add_argument("--lab_key", required=True, help="Lab alias key")
    parser.add_argument("--award", required=True, help="Lab award")
    parser.add_argument("--analysis_set_acc", required=True)
    parser.add_argument("--atac_r1_acc", required=True)
    parser.add_argument("--atac_r2_acc", required=True)
    parser.add_argument("--atac_seqspec_acc", required=True)
    parser.add_argument("--rna_r1_acc", required=True)
    parser.add_argument("--rna_r2_acc", required=True)
    parser.add_argument("--rna_seqspec_acc", required=True)


    # Parse the arguments
    args = parser.parse_args()
    
    # Example usage of parsed arguments
    #print("ATAC BAM:", args.atac_bam)
    #print("ATAC BAM Index:", args.atac_bam_index)
    #print("ATAC Fragment:", args.atac_fragment)
    #print("ATAC Fragment Index:", args.atac_fragment_index)
    print("ATAC QC:", args.atac_qc)
    #print("ATAC Metadata List:", args.atac_mm_list)
    #print("RNA H5AD:", args.rna_h5ad)
    #print("RNA KB Tar:", args.rna_kb_tar)
    print("RNA QC:", args.rna_qc)
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

    #ATAC QC
    payload = {}
    payload[Connection.PROFILE_KEY] = "single_cell_atac_seq_quality_metric"
    payload["quality_metric_of"] = [args.lab_key + args.analysis_set_acc + "_bam", 
                                    args.lab_key + args.analysis_set_acc + "_fragments"]
    payload["lab"] = args.lab 
    payload["award"] = args.award
    payload["analysis_step_version"] = "/analysis-step-versions/39c0498d-91f6-42de-8896-2fab1403f032/"

    conn.post(payload)

    #RNA QC
    payload = {}
    payload[Connection.PROFILE_KEY] = "single_cell_rna_seq_quality_metric"
    payload["quality_metric_of"] = [args.lab+":" + args.analysis_set_acc + "_matrix_file", 
                                    args.lab+":" + args.analysis_set_acc + "_rna_kb_output_folder_tar_gz"]
    payload["lab"] = args.lab 
    payload["award"] = args.award
    payload["analysis_step_version"] = "/analysis-step-versions/9c457b9f-fc6d-4cf1-b249-218827e9b449/"

    conn.post(payload)

if __name__ == "__main__":
    main()