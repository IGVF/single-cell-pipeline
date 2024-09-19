# This code takes in input a table containing the results of the pipeline and create a manifest file to upload the results to Synapse
# The table should contain the following columns:
# - sample_id This might have different names, but it should be the unique identifier for each sample
# - ATAC_barcode
# -	ATAC_fastq_R1
# - ATAC_fastq_R2
# - Barcode_inclusion_list_ATAC
# - Barcode_inclusion_list_RNA
# - Chemistry
# - Genome
# - RNA_fastq_R1
# - RNA_fastq_R2
# - Subpool

# - atac_bam
# - atac_bam_log

# - atac_chromap_barcode_metadata
# - atac_filter_fragments
# - atac_filter_fragments_index

# - atac_snapatac2_barcode_metadata

# - csv_summary
# - html_summary

# - genome_tsv

# - joint_barcode_metadata

# - reference_fasta

# - rna_barcode_metadata

# - rna_aggregated_counts_h5ad
# - rna_kb_output
# - rna_log
# - rna_mtx_tar
# - rna_mtxs_h5ad

# - seqspec
# - seqspec_atac_index
# - seqspec_atac_onlist
# - seqspec_rna_index
# - seqspec_rna_onlist
# The manifest will contain the following columns:
# - path
# - parent
# - name
# - used
# - executed
# - activityName

import argparse
import csv
import math
import os
import pandas as pd
import synapseclient
import synapseutils
from synapseclient import Folder

to_be_uploaded_file = ["seqspec_atac_onlist", "seqspec_rna_onlist", "atac_bam", "atac_bam_log", "atac_chromap_barcode_metadata", "atac_filter_fragments", "atac_filter_fragments_index", "atac_snapatac2_barcode_metadata", "csv_summary", "html_summary", "joint_barcode_metadata",  "rna_barcode_metadata", "rna_aggregated_counts_h5ad", "rna_kb_output", "rna_log", "rna_mtx_tar", "rna_mtxs_h5ad"]


def get_activity_and_script(colname):
    if colname == 'seqspec_atac_onlist':
        return ["\"Generating onlist for ATAC\"", "\"https://github.com/IGVF/single-cell-pipeline/blob/main/tasks/task_seqspec_extract.wdl\""]
    if colname == 'seqspec_rna_onlist':
        return ["\"Generating onlist for RNA\"", "\"https://github.com/IGVF/single-cell-pipeline/blob/main/tasks/task_seqspec_extract.wdl\""]
    if colname == 'atac_bam':
        return ["\"Align ATAC and return bam and log\"", "\"https://github.com/IGVF/single-cell-pipeline/blob/main/tasks/task_chromap_bam.wdl\""]
    if colname == 'atac_bam_log':
        return ["\"Align ATAC and return bam log\"", "\"https://github.com/IGVF/single-cell-pipeline/blob/main/tasks/task_chromap_bam.wdl\""]
    if colname == 'atac_chromap_barcode_metadata':
        return ["\"Barcode metadata statistics from Chromap\"", "\"https://github.com/IGVF/single-cell-pipeline/blob/main/tasks/task_chromap.wdl\""]
    if colname == 'atac_filter_fragments':
        return ["\"Raw fragment file from Chromap\"", "\"https://github.com/IGVF/single-cell-pipeline/blob/main/tasks/task_chromap.wdl\""]
    if colname == 'atac_filter_fragments_index':
        return ["\"Index for raw fragment file from Chromap\"", "\"https://github.com/IGVF/single-cell-pipeline/blob/main/tasks/task_chromap.wdl\""]
    if colname == 'atac_snapatac2_barcode_metadata':
        return ["\"Filtered barcode metadata statistics from SnapATAC2\"", "\"https://github.com/IGVF/single-cell-pipeline/blob/main/tasks/task_qc_atac.wdl\""]
    if colname == 'csv_summary':
        return ["\"Pipeline summary in CSV file format\"", "\"https://github.com/IGVF/single-cell-pipeline/blob/main/tasks/task_html_report.wdl\""]
    if colname == 'html_summary':
        return ["\"Pipeline summary in HTML file format\"", "\"https://github.com/IGVF/single-cell-pipeline/blob/main/tasks/task_html_report.wdl\""]
    if colname == 'joint_barcode_metadata':
        return ["\"Joint barcode metadata\"", "\"https://github.com/IGVF/single-cell-pipeline/blob/main/tasks/task_joint_qc.wdl\""]
    if colname == 'reference_fasta':
        return ["\"Reference fasta file\"", "\"input file\""]
    if colname == 'rna_barcode_metadata':
        return ["\"Barcode metadata statistics from KB\"", "\"https://github.com/IGVF/single-cell-pipeline/blob/main/tasks/task_kb_count.wdl\""]
    if colname == 'rna_aggregated_counts_h5ad':
        return ["\"RNA count matrix in h5ad format\"", "\"https://github.com/IGVF/single-cell-pipeline/blob/main/tasks/task_kb_count.wdl\""]
    if colname == 'rna_kb_output':
        return ["\"Total outptu from KB in TAR format\"", "\"https://github.com/IGVF/single-cell-pipeline/blob/main/tasks/task_kb_count.wdl\""]
    if colname == 'rna_log':
        return ["\"Log file from KB\"", "\"https://github.com/IGVF/single-cell-pipeline/blob/main/tasks/task_log_rna.wdl\""]
    if colname == 'rna_mtx_tar':
        return ["\"RNA count matrices in TAR format\"", "\"https://github.com/IGVF/single-cell-pipeline/blob/main/tasks/task_kb_count.wdl\""]
    if colname == 'rna_mtxs_h5ad':
        return ["\"RNA count matrices in h5ad format\"", "\"https://github.com/IGVF/single-cell-pipeline/blob/main/tasks/task_kb_count.wdl\""]
    return ["", ""]


def get_used(row, colname, local_root, remote_root):
    if colname == 'seqspec_atac_onlist':
        return f"\"{clean_string(row['ATAC_barcode'],local_root,remote_root)};{clean_string(row['ATAC_fastq_R1'],local_root,remote_root)};{clean_string(row['ATAC_fastq_R2'],local_root,remote_root)};{clean_string(row['Barcode_inclusion_list_ATAC'], 'https://storage.googleapis.com/', 'gs://')};{row['seqspec']}\""
    if colname == 'seqspec_rna_onlist':
        return f"\"{clean_string(row['RNA_fastq_R1'],local_root,remote_root)};{clean_string(row['RNA_fastq_R2'],local_root,remote_root)};{clean_string(row['Barcode_inclusion_list_RNA'],'https://storage.googleapis.com/','gs://')};{row['seqspec']}\""
    if colname == 'atac_bam':
        return f"\"{clean_string(row['ATAC_barcode'],local_root,remote_root)};{clean_string(row['ATAC_fastq_R1'],local_root,remote_root)};{clean_string(row['ATAC_fastq_R2'],local_root,remote_root)};{row['seqspec_atac_onlist'].replace(remote_root, local_root)}\""
    if colname == 'atac_bam_log':
        return f"\"{clean_string(row['ATAC_barcode'],local_root,remote_root)};{clean_string(row['ATAC_fastq_R1'],local_root,remote_root)};{clean_string(row['ATAC_fastq_R2'],local_root,remote_root)};{row['seqspec_atac_onlist'].replace(remote_root, local_root)}\""
    if colname == 'atac_chromap_barcode_metadata':
        return f"\"{clean_string(row['ATAC_barcode'],local_root,remote_root)};{clean_string(row['ATAC_fastq_R1'],local_root,remote_root)};{clean_string(row['ATAC_fastq_R2'],local_root,remote_root)};{row['seqspec_atac_onlist'].replace(remote_root, local_root)}\""
    if colname == 'atac_filter_fragments':
        return f"\"{clean_string(row['ATAC_barcode'],local_root,remote_root)};{clean_string(row['ATAC_fastq_R1'],local_root,remote_root)};{clean_string(row['ATAC_fastq_R2'],local_root,remote_root)};{row['seqspec_atac_onlist'].replace(remote_root, local_root)}\""
    if colname == 'atac_filter_fragments_index':
        return f"\"{clean_string(row['ATAC_barcode'],local_root,remote_root)};{clean_string(row['ATAC_fastq_R1'],local_root,remote_root)};{clean_string(row['ATAC_fastq_R2'],local_root,remote_root)};{row['seqspec_atac_onlist'].replace(remote_root, local_root)}\""
    if colname == 'atac_snapatac2_barcode_metadata':
        return f"\"{clean_string(row['atac_filter_fragments'], local_root, remote_root)};{clean_string(row['atac_filter_fragments_index'], local_root, remote_root)};{clean_string(row['atac_chromap_barcode_metadata'], local_root, remote_root)}\""
    if colname == 'joint_barcode_metadata':
        return f"\"{clean_string(row['atac_snapatac2_barcode_metadata'], local_root, remote_root)};{clean_string(row['rna_barcode_metadata'], local_root, remote_root)}\""
    if colname == 'rna_barcode_metadata':
        return f"\"{clean_string(row['RNA_fastq_R1'],local_root,remote_root)};{clean_string(row['RNA_fastq_R2'],local_root,remote_root)};{row['seqspec_rna_onlist'].replace(remote_root, local_root)}\""
    if colname == 'rna_aggregated_counts_h5ad':
        return f"\"{clean_string(row['RNA_fastq_R1'],local_root,remote_root)};{clean_string(row['RNA_fastq_R2'],local_root,remote_root)};{row['seqspec_rna_onlist'].replace(remote_root, local_root)}\""
    if colname == 'rna_kb_output':
        return f"\"{clean_string(row['RNA_fastq_R1'],local_root,remote_root)};{clean_string(row['RNA_fastq_R2'],local_root,remote_root)};{row['seqspec_rna_onlist'].replace(remote_root, local_root)}\""
    if colname == 'rna_log':
        return f"\"{clean_string(row['RNA_fastq_R1'],local_root,remote_root)};{clean_string(row['RNA_fastq_R2'],local_root,remote_root)};{row['seqspec_rna_onlist'].replace(remote_root, local_root)}\""
    if colname == 'rna_mtx_tar':
        return f"\"{clean_string(row['RNA_fastq_R1'],local_root,remote_root)};{clean_string(row['RNA_fastq_R2'],local_root,remote_root)};{row['seqspec_rna_onlist'].replace(remote_root, local_root)}\""
    if colname == 'rna_mtxs_h5ad':
        return f"\"{clean_string(row['RNA_fastq_R1'],local_root,remote_root)};{clean_string(row['RNA_fastq_R2'],local_root,remote_root)};{row['seqspec_rna_onlist'].replace(remote_root, local_root)}\""
    return None


def create_manifest_row(row, col, parent, basename_string, local_root, remote_root):
    activity, script = get_activity_and_script(col)
    used = get_used(row, col, local_root, remote_root)
    name_to_use = os.path.basename(row[col]).replace(".filter", "")
    if col == 'seqspec_rna_onlist':
        name_to_use = f"{basename_string}.rna.onlist.txt.gz"
    if col == 'seqspec_atac_onlist':
        name_to_use = f"{basename_string}.atac.onlist.txt.gz"
    manifest_row = {
        'path': row[col].replace(remote_root, local_root),
        'parent': parent,
        'name': name_to_use,
        'used': used if not used is None else "\"\"",
        'executed': script,
        'activityName': activity
    }
    return manifest_row


def clean_string(string, local_root, remote_root):
    return string.replace("]", "").replace("[", "").replace(",", ";").replace('"', '').replace(remote_root, local_root)


def main():
    parser = argparse.ArgumentParser(description='Upload Y3 results to Synapse')
    parser.add_argument('--table', help='Path to the table containing the results')
    parser.add_argument('--project', help='Synapse project ID')
    parser.add_argument('--output', help='Path to the output manifest file')
    parser.add_argument('--remote-root', help='Path to the remote directory where the files are stored')
    parser.add_argument('--local-root', help='Path to the remote directory where the files are stored')
    parser.add_argument('--dry-run', help='Do not upload the files to Synapse', action='store_true')
    args = parser.parse_args()

    syn = synapseclient.Synapse()
    syn.login(silent=True)

    table = pd.read_csv(args.table, sep='\t')

    manifest = []
    for i, row in table.iterrows():
        basename_string = row['Subpool']
        basename_string = clean_string(basename_string, args.local_root, args.remote_root)
        # Create a folder for the data
        folder = syn.findEntityId(name=basename_string, parent=args.project)
        if not folder:
            folder = Folder(name=basename_string, parent=args.project)
            folder = syn.store(folder)
        for col in table.columns:
            if (isinstance(row[col], float) and math.isnan(row[col])) or col not in to_be_uploaded_file:
                continue
            manifest_row = create_manifest_row(row, col, folder, basename_string, args.local_root, args.remote_root)
            if manifest_row["parent"]:
                manifest.append(manifest_row)

    # Save the manifest to a tsv file
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=manifest[0].keys(), delimiter='\t', quoting=csv.QUOTE_NONE, quotechar=None, escapechar='\\')
        writer.writeheader()
        writer.writerows(manifest)

#    synapseutils.syncToSynapse(
#        syn=syn,
#        manifestFile=args.output,
#        sendMessages=False, dryRun=True
#    )

# Run the main function
if __name__ == '__main__':
    main()

