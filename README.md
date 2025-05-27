# IGVF Single-cell pipeline

Welcome to the IGVF Single Cell Pipeline repository.
The latest version of the pipeline is the [v1.1](https://github.com/IGVF/single-cell-pipeline/tree/v1.1) and it is based on the
[atomic-workflows v1.1](https://github.com/IGVF/atomic-workflows/tree/v1.1).
The documentation for v1.1 can be found [here](https://docs.google.com/document/d/1NgNYDduZsThKTyND8DI1DIMwiG9q-Rt462377_NZXis/edit?tab=t.0#heading=h.9ecc41kilcvq).

## Table of Contents

- [Introduction](#introduction)
- [Documentation](#documentation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Introduction

This pipeline has been designed by the IGVF single-cell focus group for automated pre-processing of single-cell sequencing data. 
It is a WDL wrapper for the python scripts present in [atomic-workflows](https://github.com/IGVF/atomic-workflows/tree/v1.1).It supports native execution on the Terra/AnVIL platform but can also be executed on both compute clusters with job submission engines as well as standalone machines, with built-in parallelization and distributed computing capabilities. 

The pipeline can be run with raw FASTQ files and a seqspec file describing them as inputs. The pipeline can handle ATAC and RNA modalities. Taking a seqspec in input makes the pipeline virtually suitable to process every library independently from the assay, protocol, library preparation or sequencing strategy used.
The ATAC and RNA are processed and QC independently. The main outputs of the pipeline are a fragment file for the ATAC and a count matrix for RNA.

![Pipeline Overview](docs/images/pipeline_flowchart.png)


## Documentation

For more detailed documentation of the inputs and outputs of each task, please refer to our the IGVF [atomic-workflows](https://github.com/IGVF/atomic-workflows/tree/v1.1).


## Usage

To run the pipeline on Terra/Anvil go to the [dockstore page](https://dockstore.org/workflows/github.com/IGVF/single-cell-pipeline/IGVF-single_cell_pipeline:main?tab=info) and select export to Terra/AnVIL


To run the pipeline locally we recommend using the IGVF [atomic-workflows](https://github.com/IGVF/atomic-workflows/tree/v1.1).

## Input Files

Accepted input formats are: `local path`, `https` url, Google Storage (`gs`) url


### Common Inputs

- **create_onlist_mapping** *(Type: Boolean)*: Whether to create an on-list mapping. This is useful for 10x multiomic assays. Default is `false`.
- **prefix** *(Type: String)*: Prefix for the analysis set.
- **igvf_credentials** *(Type: File?)*: Optional IGVF credentials file.
- **subpool** *(Type: String?)*: Subpool to address. Default is `"none"`.
- **genome_tsv** *(Type: File)*: TSV formatted file containing genome references and annotations.

### ATAC-specific Inputs

- **atac_read1** *(Type: Array[File])*: First read file(s) for ATAC sequencing data. If you are only processing RNA use `[]`
- **atac_read2** *(Type: Array[File])*: Second read file(s) for ATAC sequencing data. If you are only processing RNA use `[]`
- **fastq_barcode** *(Type: Array[File])*: Barcode file(s) for ATAC sequencing data. If you are only processing RNA use `[]`
- **atac_barcode_inclusion_list** *(Type: File)*: Inclusion list for ATAC barcodes.
- **chromap_genome_index_tar_gz** *(Type: File?)*: Optional Chromap genome index tarball.
- **genome_fasta** *(Type: File?)*: Optional genome FASTA file.
- **atac_read_format** *(Type: String)*: Format of the ATAC read files.

### RNA-specific Inputs

- **rna_read1** *(Type: Array[File])*: First read file(s) for RNA sequencing data. If you are only processing ATAC use `[]`
- **rna_read2** *(Type: Array[File])*: Second read file(s) for RNA sequencing data. If you are only processing ATAC use `[]`
- **fastq_barcode_rna** *(Type: Array[File])*: Barcode file(s) for RNA sequencing data. Default is `[]`.
- **rna_barcode_inclusion_list** *(Type: File)*: Inclusion list for RNA barcodes.
- **kb_mode** *(Type: String)*: Mode for RNA processing. Default is `"nac"`.
- **rna_read_format** *(Type: String)*: Format of the RNA read files.
- **kb_genome_index_tar_gz** *(Type: File?)*: Optional KB genome index tarball.
- **rna_replacement_list** *(Type: File?)*: Optional replacement list for RNA. This is specific for kallisto-bustool.


## Output Files

### RNA

| Name | Suffix | Description |
|------|--------|-------------|
| `rna_kb_h5ad` | .rna_kb_h5ad | [Raw] h5ad file containing RNA count matrices |
| `rna_kb_output_folder_tar_gz` | .rna_kb_output_folder.tar.gz | [Raw] Tarball containing all output files from the RNA processing |
| `rna_kb_run_info_json` | .rna_kb_run_info.json | [Raw] JSON file with run information for RNA processing |
| `rna_kb_library_qc_metrics_json` | .rna_kb_library_qc_metrics.json | [Raw] JSON file with library QC metrics for RNA |
| `rna_kb_parameters_json` | .rna_kb_parameters.json | [Raw] JSON file with parameters used for RNA processing |

### ATAC

| Name | Suffix | Description |
|------|--------|-------------|
| `atac_bam` | .atac_chromap_bam | [Raw] Aligned BAM file from Chromap for ATAC |
| `atac_bam_index` | .atac_chromap_bam_index | [Raw] BAM index file for ATAC |
| `atac_bam_summary_stats` | .atac_chromap_bam_summary | [Raw] Summary statistics generated by Chromap for the aligned BAM file |
| `atac_fragments` | .atac_fragments | [Raw] Fragment file generated from Chromap for ATAC |
| `atac_fragments_index` | .atac_fragments_index | [Raw] Index file for the ATAC fragment file |
| `atac_fragments_log` | .atac_fragments_alignment_stats | [Raw] Log file with alignment statistics from Chromap for ATAC fragments |
| `atac_fragments_metrics` | .atac_fragments_qc_metrics | [Raw] QC metrics from chromap for ATAC fragments |
| `atac_fragments_alignment_stats` | .atac_fragments_alignment_stats | [Raw] Alignment statistics from chromap for ATAC fragments |
| `atac_fragments_barcode_summary` | .atac_fragments_barcode_summary | [Raw] Summary of barcode statistics from chromap for ATAC fragments |


## Contributing

We welcome contributions to improve this pipeline. Please fork the repository and submit a pull request with your changes. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
