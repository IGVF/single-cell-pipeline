# This script takes in a multireport from the IGVF portal and generates a table that can be used as input to the uniform processing pipeline. 

## WORK IN PROGRESS

## The output table will have 10 columns (2 columns describing seqspec not implemented yet)

## [OUTPUT TABLE DESCRIPTION] Columns: Row ID (biosam_x_subpool accession string), ATAC_MM (list of measurementSet accessions), RNA_MM (list of measurementSet accessions), ATAC_R1_files (array of SequenceFile accessions), ATAC_R2_files (array of SequenceFile accessions), RNA_R1_files (array of SequenceFile accessions), RNA_R2_files (array of SequenceFile accessions), ATAC_seqspecs (NOT IMPLEMENTED YET, array of accessions), RNA_seqspecs (NOT IMPLEMENTED YET, array of accessions).

import igvf_utils as iu
from igvf_utils.connection import Connection 
import pandas as pd
import numpy as np #don't really need this, replace np code with base code
import sys

pd.options.mode.chained_assignment = None


#Input table is report from IGVF portal. For example, a user looking to analyze the Share-seq bone marrow dataset from the Buenrostro-Bernstein group can generate the table from the below link

#https://api.data.igvf.org/multireport/?type=MeasurementSet&lab.title=Jason+Buenrostro%2C+Broad&donors.taxa=Homo+sapiens&samples.sample_terms.term_name=bone+marrow&field=%40id&field=accession&field=aliases&field=alternate_accessions&field=award&field=donors&field=lab&field=samples&field=summary&field=donors.taxa&field=files&field=multiome_size&field=assay_term&limit=all


df = pd.read_csv(sys.argv[1], sep = "\t", skiprows=1)

#connect to portal instance
if (sys.argv[2] == "prod"):
    conn = Connection("api.data.igvf.org")
else:
    conn = Connection("api.sandbox.igvf.org")

#first create unique list of biosample x subpool (Samples of measurementSet)
biosam_x_subpool_list = list()

for index, mm_row in df.iterrows():
    print(mm_row["Accession"])  
    
    biosam_x_subpool_list.append( conn.get(mm_row["Accession"])["samples"][0]["accession"] )

biosam_x_subpool_list = list(set(biosam_x_subpool_list))

#next, using each biosam_x_subpool, create the processing table. 
## [FUTURE] In theory, you could just query for seqspecs and that is all that's needed for the pipeline (a table of biosam_x_subpool with array of RNA seqspec and array of ATAC seqspec) 

#create empty output table [TODO: add seqspec columns]

processing_table = pd.DataFrame({"row_ID": biosam_x_subpool_list})

processing_table["ATAC_MM"] = np.empty((len(processing_table), 0)).tolist()

processing_table["RNA_MM"] = np.empty((len(processing_table), 0)).tolist()

processing_table["ATAC_R1"] = np.empty((len(processing_table), 0)).tolist()
processing_table["ATAC_R2"] = np.empty((len(processing_table), 0)).tolist()

processing_table["RNA_R1"] = np.empty((len(processing_table), 0)).tolist()
processing_table["RNA_R2"] = np.empty((len(processing_table), 0)).tolist()

processing_table["taxa"] = "NA"

print(processing_table)

for index, sam_row in processing_table.iterrows():
    print(sam_row["row_ID"])  
    
    #query portal for biosam_x_subpool
    biosam_x_subpool = conn.get(sam_row["row_ID"])
    
    tdf = pd.DataFrame(biosam_x_subpool["file_sets"])[['summary','accession']].groupby("summary").agg(list).T
    
    #this had to be harcoded to search for ATAC/RNA. DACC might have better solution
    processing_table["ATAC_MM"][index] = tdf["single-cell ATAC-seq (SHARE-seq)"]["accession"]
    
    processing_table["RNA_MM"][index] = tdf["single-cell RNA sequencing assay (SHARE-seq)"]["accession"]
    
    processing_table["taxa"][index] = biosam_x_subpool['taxa']
    
    atac_r1_file = []
    atac_r2_file = []
    for acc in tdf["single-cell ATAC-seq (SHARE-seq)"]["accession"]:
        #get mm_set
        mm = conn.get(acc)
        files_df = pd.DataFrame(mm["files"])
        
        #keep fastqs only
        files_df = files_df[files_df['file_format'] == 'fastq']
        files_df["sequencing_run"] = 0
        files_df["illumina_read_type"] = 0
        
        for file_index, file_row in files_df.iterrows():
            
            f = conn.get(file_row['accession'])
        
            files_df["sequencing_run"][file_index] = f["sequencing_run"]
            
            files_df["illumina_read_type"][file_index] = f["illumina_read_type"]
        
        files_df = files_df.sort_values(by=['sequencing_run','illumina_read_type'])
        
        atac_r1_file = list(files_df['accession'][::2])
        atac_r2_file = list(files_df['accession'][1::2])
        
    processing_table["ATAC_R1"][index] = atac_r1_file
    processing_table["ATAC_R2"][index] = atac_r2_file
        
    rna_r1_file = []
    rna_r2_file = []
    for acc in tdf["single-cell RNA sequencing assay (SHARE-seq)"]["accession"]:
        #get mm_set
        mm = conn.get(acc)
        files_df = pd.DataFrame(mm["files"])
        
        #keep fastqs only
        files_df = files_df[files_df['file_format'] == 'fastq']
        files_df["sequencing_run"] = 0
        files_df["illumina_read_type"] = 0
        
        for file_index, file_row in files_df.iterrows():
            
            f = conn.get(file_row['accession'])
            
            files_df["sequencing_run"][file_index] = f["sequencing_run"]
            
            files_df["illumina_read_type"][file_index] = f["illumina_read_type"]
        
        files_df = files_df.sort_values(by=['sequencing_run','illumina_read_type'])
        
        rna_r1_file = list(files_df['accession'][::2])
        rna_r2_file = list(files_df['accession'][1::2])
    
    processing_table["RNA_R1"][index] = rna_r1_file
    processing_table["RNA_R2"][index] = rna_r2_file
    
processing_table.to_csv("processing_table.tsv", sep = "\t")
    
    