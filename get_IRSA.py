from astroquery.ipac.irsa import Irsa
from astropy.coordinates import SkyCoord
import astropy.units as u
import csv
import numpy as np


#  extracting the list of catalogs with 2mass data
# "2mass" is not in catalog identifiers, but in catalog descriptions
# print({k for k, v in Irsa.list_catalogs().items() if "2mass" in k.lower() or "2mass" in v.lower()})

with open('black_holes.csv', 'r') as f:
    csv_reader = csv.reader(f)
    header = next(csv_reader)
    black_holes = list(csv_reader)

results = []
for bh in black_holes:
    name = bh[0]
    ra_bh = bh[1]
    dec_bh = bh[2]

    # Convert black hole coordinates to decimal degrees
    bh_coord = SkyCoord(ra_bh, dec_bh, unit=(u.hourangle, u.deg), frame='icrs')
    ra = bh_coord.ra.deg
    dec = bh_coord.dec.deg


    table = Irsa.query_region(
        bh_coord,
        radius=30 * u.arcsec,
        catalog="fp_psc",   # selecting relevant catalog --> point source catalog
        spatial='Cone'
    )

    # saving only the data connected to the K band; 
    # ph_qual accesses quality of data: goes from A (best)to E (worst), U is upper limit
    # cc_flg is the contamination flag: 0 is unaffected
    # K band is the 3rd character in these tuples
    cols = ["ra", "dec", "k_m", "k_msigcom", "ph_qual", "cc_flg"]

    if (table is None) or (len(table) == 0):
        results.append({
            'bh_name': name,
            'bh_ra_deg': ra_bh,
            'bh_dec_deg': dec_bh,
            'k_object_ra_deg': None,
            'k_object_dec_deg': None,
            'sep_arcsec': None,  
            'k_mag': None,
            'k_mag_err': None,
            'ph_qual': None,
            'cc_flg': None
        })
        continue

    tab = table[cols]

    tab = tab.filled(np.nan)

# only sorting by magnitude, could be done including quality as well
    tab.sort('k_m')

    best = tab[0]
    k_object_coord = SkyCoord(best['ra'], best['dec'], unit='deg', frame='icrs')
    sep_arcsec = bh_coord.separation(k_object_coord).arcsec 

    results.append({
        'bh_name': name,
        'bh_ra_deg': ra,
        'bh_dec_deg': dec,
        'k_object_ra_deg': float(best['ra']) if best['ra'] else None,
        'k_object_dec_deg': float(best['dec']) if best['dec'] else None,
        'sep_arcsec': sep_arcsec if sep_arcsec else None,
        'k_mag': float(best['k_m']) if best['k_m'] else None,
        'k_mag_err': float(best['k_msigcom']) if best['k_msigcom'] else None,
        'ph_qual': best['ph_qual'][2],
        'cc_flg': best['cc_flg']
    })

out_header = [
    'bh_name', 'bh_ra_deg', 'bh_dec_deg',
    'k_object_ra_deg', 'k_object_dec_deg', 'sep_arcsec',
    'k_mag', 'k_mag_err', 'ph_qual', 'cc_flg'
]

with open('black_holes_with_2mass_k.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=out_header)
    writer.writeheader()
    for r in results:
        writer.writerow(r)



