
# python download_online_text_file.py

import requests
from bs4 import BeautifulSoup
import pandas as pd
import os

# URL of the directory page containing the text file links
url = "https://www.ncei.noaa.gov/pub/data/uscrn/products/subhourly01/2023/" #just change the final number here to pull from other years
data_FOLDER = 'Raw Data'

# Create a directory to save the files
os.makedirs(data_FOLDER, exist_ok=True)

# Fetch the directory page
response = requests.get(url)
soup = BeautifulSoup(response.content, 'html.parser')

# Find all links to text files
links = soup.find_all('a', href=True)
text_links = [link['href'] for link in links if link['href'].endswith('.txt')]

total_files = len(text_links) # find the total number of text files to download

for idx, link in enumerate(text_links):
    # Full URL of the text file
    file_url = url + link

    # Fetch the text file
    file_response = requests.get(file_url)

    # Extract file name from the URL
    file_name = link

    # Save as text file
    with open(f'{data_FOLDER}/{file_name}', 'wb') as file:
        file.write(file_response.content)

    # Save as parsed Excel file with the same file name
    df = pd.read_csv(file_url, delimiter=r'\s+', header=None)
    excel_file_name = file_name.replace('.txt', '.xlsx')
    df.to_excel(f'{data_FOLDER}/{excel_file_name}', index=False, header=False)

    # Display the progress with the name and number
    print(f"Finished downloading and converting {file_name} ({idx + 1}/{total_files}). {total_files - idx - 1} files left.")

print("All files downloaded and saved.")
