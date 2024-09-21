# Description: This script is used to upload files to Synapse in batch. It reads a manifest file and uploads the files to Synapse.

import synapseclient
import synapseutils

syn = synapseclient.Synapse()
syn.login(silent=True)

failed_ids = []

with open('manifest_team7_resume.csv', 'r') as f:
    header = f.readline().strip()
    temp = None
    current_id = None
    for line in f:
        id = line.split('\t')[1]
        if current_id is None:
            current_id = id
            temp = open('temp.tsv', 'w')
            temp.write(header + '\n')
        elif id != current_id:
            temp.close()
            try:
                print("Uploading file with id: " + current_id)
                synapseutils.syncToSynapse(syn=syn, manifestFile="temp.tsv", sendMessages=False)
                print("Successfully uploaded file with id: " + current_id)
            except Exception as e:
                failed_ids.append(current_id)
            current_id = id
            temp = open('temp.tsv', 'w')
            temp.write(header + '\n')
        temp.write(line)
    if temp is not None:
        temp.close()
        try:
            print("Uploading file with id: " + current_id)
            synapseutils.syncToSynapse(syn=syn, manifestFile="temp.tsv", sendMessages=False)
            print("Successfully uploaded file with id: " + current_id)
        except Exception as e:
            failed_ids.append(current_id)

# Save the list of failed ids to a file
with open('failed_ids.txt', 'w') as f:
    for id in failed_ids:
        f.write(id + '\n')
