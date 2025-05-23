FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN pip install requests

# Copy the download_with_credentials.py script into the image
COPY download_with_credentials.py /usr/local/bin/download_with_credentials.py

# Make the script executable
RUN chmod +x /usr/local/bin/download_with_credentials.py

# Set the entrypoint to bash
ENTRYPOINT ["/bin/bash"]
