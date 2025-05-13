version 1.0

# TASK
# Create a file with QC metrics
# Gather information from log files


task log_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard IGVF pipeline: log atac task'
    }

    input {
        # This function takes as input the necessary log files and extracts
        # the quality metrics
        File? alignment_log
        File? barcode_log
        String? prefix = "sample"
    }

    command <<<
        # Formatting the output of chromap and extracting statistics
        grep "Number of" ~{alignment_log} | grep -v threads| tr -d '.' | sed 's/ /_/g' | sed 's/:_/,/g'> qc_metrics.csv
        grep "#" ~{alignment_log}  | sed 's/, /\n/g' | tr -d '# ' | sed 's/:/,/g' | tr -d '.' >> qc_metrics.csv
        # Compute the percentage of duplicates from the barcode log file.
        awk -v FS="," 'NR>1{total+=$2; dups+=$3; unmapped+=$4; lowmapq+=$5}END{printf "percentage_duplicates,%5.1f", 100*dups/(total-unmapped-lowmapq)}' ~{barcode_log} >> qc_metrics.csv
        # Convert the csv to JSON
        python -c "import csv, json; f = open('qc_metrics.csv', 'r'); reader = csv.reader(f); data = {row[0]: int(row[1].strip()) if row[1].isdigit() else float(row[1].strip()) for row in reader}; print(json.dumps(data, indent=4))" > ~{prefix}_qc_metrics.json

    >>>
    output {
        File atac_statistics_json = "~{prefix}_qc_metrics.json"
    }

    runtime {
        docker: 'docker.io/igvf/chromap:v1'
    }
    parameter_meta {
        alignment_log: {
            description: 'ATAC alignment log file',
            help: 'Log file from ATAC alignment step.'
        }
        barcode_log: {
            description: 'ATAC dups log file',
            help: 'Barcode log file from ATAC alignment step.'
        }
    }
}