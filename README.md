# IGVF Single-cell pipeline

Welcome to the IGVF Single Cell Pipeline repository. 

## Table of Contents

- [Introduction](#introduction)
- [Documentation](#documentation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Introduction

This pipeline had been designed by the IGVF single-cell focus group for automated pre-processing and quality control of single-cell sequencing data. It supports native execution on the Terra/AnViL platform but can also be executed on both compute clusters with job submission engines as well as standalone machines, with built-in parallelization and distributed computing capabilities. 

The pipeline can be run with raw FASTQ files and a seqspec file describing them as inputs. The pipeline can handle ATAC and RNA modalities. Taking a seqspec in input makes the pipeline virtually suitable to process every library independently from the assay, protocol, library preparation or sequencing strategy used.
The RNA and ATAC are processed and QC independently. If the dataset is multi-modal, a final step of joint QC is performed. The main outputs of the pipeline are a fragment file for the ATAC and a count matrix for RNA along side QC metrics and a HTML report.

![Pipeline Overview](docs/images/pipeline_flowchart.png)


## Documentation

For more detailed documentation, please refer to our [Google Doc](https://docs.google.com/document/d/1NgNYDduZsThKTyND8DI1DIMwiG9q-Rt462377_NZXis/edit).

## Usage

To run the pipeline on Terra/Anvil go to the [dockstore page](https://dockstore.org/workflows/github.com/IGVF/single-cell-pipeline/IGVF-single_cell_pipeline:main?tab=info) and select export to Terra/Anvil


To run the pipeline locally you can refere to these [instructions](docs/install-pipeline-locally.org).

## Contributing

We welcome contributions to improve this pipeline. Please fork the repository and submit a pull request with your changes. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
