import hashlib
from google.cloud import storage
import subprocess
import re
import base64 

def get_md5sum(bucket_name, object_name):
    
    #running gsutil stat from command line
    command = f"gsutil stat gs://{bucket_name}/{object_name}"
    
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    
    output = result.stdout
    md5_match = re.search(r'Hash \(md5\):\s+(\S+)', output)
    base64_md5 = md5_match.group(1)
    hex_md5 = base64.b64decode(base64_md5).hex()
    return hex_md5

def calculate_md5(bucket_name, object_name):
    # Create a storage client
    client = storage.Client()

    # Get the bucket and object
    bucket = client.get_bucket(bucket_name)
    blob = bucket.blob(object_name)

    # Download the object as a byte stream
    byte_stream = blob.download_as_string()

    # Calculate the MD5 hash of the byte stream
    md5_hash = hashlib.md5(byte_stream).hexdigest()

    return md5_hash

def copy_to_gcs(bucket_name, source_file_path, destination_blob_name):
    """
    Uploads a file to a GCS bucket.

    Args:
        bucket_name (str): The name of the GCS bucket.
        source_file_path (str): The path to the file to upload.
        destination_blob_name (str): The name of the blob in the GCS bucket.
    """
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(destination_blob_name)
    blob.upload_from_filename(source_file_path)

# Provide your bucket and object name
#bucket_name = "your-bucket-name"
#object_name = "your-object-name"

# Call the function to retrieve the MD5 hash
#md5sum = get_md5sum(bucket_name, object_name)
#print("MD5 sum:", md5sum)