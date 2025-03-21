FROM debian:latest

ENV DEBIAN_FRONTEND=noninteractive

# Install necessary packages
RUN apt-get update && apt-get install -y \
    wget \
    python3 \
    python3-pip \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Copy the download_with_credentials.py script into the image
COPY /src/python/igvf_portal/download_with_credentials.py /usr/local/bin/download_with_credentials.py

# Make the script executable
RUN chmod +x /usr/local/bin/download_with_credentials.py

# Set the entrypoint to bash
ENTRYPOINT ["/bin/bash"]
