from bs4 import BeautifulSoup
import requests
import csv
from pprint import pprint
from astropy.coordinates import SkyCoord
import astropy.units as u

url = 'https://www.astro.puc.cl/BlackCAT/transients.php'
page = requests.get(url)
soup = BeautifulSoup(page.text, 'html.parser')

# Extracting confirmed black holes
black_holes = []
data = soup.find_all('tr', style="font-weight:bold")

for row in data:
    row_data = row.find_all('td')
    name = " ".join(row_data[1].text.replace("\n", " ").split())
    ra = row_data[2].text.strip()
    dec = row_data[3].text.strip()

    link_tag = row_data[1].find('a')
    relative_link = link_tag['href']
    link = f"https://www.astro.puc.cl/BlackCAT/{relative_link}"

    # Debugging: Seeing if hyperlink present 
    # print(f"Black Hole: {name}, Link: {link}")

    ks_quiescent = "no data"
    
    detail_page = requests.get(link)
    detail_soup = BeautifulSoup(detail_page.text, 'html.parser')

    table = detail_soup.find('table', {'id': 'magnitudes'})
    if table:
        td_elements = table.find_all('td')

        for td in td_elements:
            text = td.text.strip()
            if "K" in text:  # Look for both "K" and "Ks" magnitudes
                ks_quiescent = text
                break

# Debugging 
    #bh_coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
    #ra = bh_coord.ra.deg
    #dec = bh_coord.dec.deg
    black_holes.append({
        'name': name,
        'ra': ra,
        'dec': dec,
        'brightness_quiescent': ks_quiescent
    })

pprint(black_holes)

with open('black_holes.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['name', 'ra', 'dec', 'brightness_quiescent'])
    writer.writeheader()
    writer.writerows(black_holes)