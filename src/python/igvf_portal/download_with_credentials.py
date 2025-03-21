from pathlib import Path

import os
import requests
import sys


def download_file_via_https(href_url: str, output_dir: str = '.') -> None:
    """Download file via href url on the portal

    Args:
        href_url (str): URL of the file to download
        output_dir (str, optional): output directory to save the file.

    Returns:
        None: None
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    username = os.getenv('IGVF_API_KEY')
    password = os.getenv('IGVF_SECRET_KEY')
    local_filename = Path(output_dir) / Path(href_url).name

    with requests.Session() as session:
        if username and password:
            session.auth = (username, password)
        response = session.get(href_url, stream=True)
        response.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)


if __name__ == "__main__":
    url = sys.argv[1]
    download_file_via_https(url)
