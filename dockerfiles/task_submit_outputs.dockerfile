FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN pip install requests awscli inflection jsonschema urllib3 google-api-python-client
RUN pip install https://github.com/IGVF-DACC/igvf_utils/archive/master.zip

# Copy the download_with_credentials.py script into the image
COPY /src/python/igvf_portal/submit_outputs.py /usr/local/bin/submit_outputs.py

# Make the script executable
RUN chmod +x /usr/local/bin/submit_outputs.py

# Set the entrypoint to bash
ENTRYPOINT ["/bin/bash"]
