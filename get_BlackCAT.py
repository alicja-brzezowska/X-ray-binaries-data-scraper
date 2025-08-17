from bs4 import BeautifulSoup
import requests
import csv 
from pprint import pprint

url = 'https://www.astro.puc.cl/BlackCAT/transients.php'
page = requests.get(url)
soup = BeautifulSoup(page.text, 'html.parser')

# Extracting confirmed black holes
black_holes = []
data = soup.find_all('tr', style="font-weight:bold")


for row in data:
    row_data = row.find_all('td')
    black_holes.append({
        'name': " ".join(row_data[1].text.replace("\n", " ").split()),
        'ra': row_data[2].text.strip(),
        'dec': row_data[3].text.strip()
    })

pprint(black_holes)

with open('black_holes.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['name', 'ra', 'dec'])
    writer.writeheader()
    writer.writerows(black_holes)