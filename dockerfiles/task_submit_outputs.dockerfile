#FROM python:3.11-slim
FROM google/cloud-sdk:slim

# Install Python if needed
RUN apt-get update && apt-get install -y python3-pip

ENV DEBIAN_FRONTEND=noninteractive

RUN pip install --break-system-packages requests awscli inflection jsonschema urllib3 google-api-python-client
RUN pip install --break-system-packages https://github.com/IGVF-DACC/igvf_utils/archive/master.zip

# Copy the download_with_credentials.py script into the image
COPY /src/python/igvf_portal/submit_outputs.py /usr/local/bin/submit_outputs.py
COPY /src/python/igvf_portal/gcs_functions.py /usr/local/bin/gcs_functions.py

# Make the script executable
RUN chmod +x /usr/local/bin/submit_outputs.py

# Set the entrypoint to bash
ENTRYPOINT ["/bin/bash"]
