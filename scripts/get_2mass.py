from astroquery.ipac.irsa import Irsa
from astropy.coordinates import SkyCoord
import astropy.units as u
import csv
import numpy as np


from functions import read_csv, save_results_to_csv

def query_IRSA_catalog(bh_coord, radius_arcsec=30.0, catalog="fp_psc"):
    "Query the IRSA catalog for data around a given RA and DEC."
    radius = (radius_arcsec / 3600.0) * u.deg
    return Irsa.query_region(bh_coord, radius=radius, catalog=catalog, spatial='Cone')


def find_brightest_IRSA(table, bh_coord):
    "Given query data finds the brightest star in the K band."
    "Returns columns of the star's RA, DEC, its brigthness and error, quality factors "
    cols = ["ra", "dec", "k_m", "k_msigcom", "ph_qual", "cc_flg"]
    if (table is None) or (len(table) == 0):
        return None, None
    tab = table[cols].filled(np.nan)
    tab.sort('k_m')
    best = tab[0]
    star_coord = SkyCoord(best['ra'], best['dec'], unit='deg', frame='icrs')
    sep_arcsec = bh_coord.separation(star_coord).arcsec
    return best, sep_arcsec


def process_input_file(input_file, radius_arcsec=30.0):
    "Perform query_IRSA_catalog and find_brightest_IRSA for each row in the input file."
    "File is expected to have format: name of the object, ra, dec"

    results = []
    for row in input_file:
        name, bh_ra, bh_dec = row[0], row[1], row[2]
        # Convert black hole coordinates to degrees for queries 
        bh_coord = SkyCoord(bh_ra, bh_dec, unit=(u.hourangle, u.deg), frame='icrs')
        #bh_ra_deg = bh_coord.ra.deg
        #h_dec_deg = bh_coord.dec.deg

        # Query the 2mass catalog
        table = query_IRSA_catalog(bh_coord, radius_arcsec=radius_arcsec, catalog="fp_psc")
        if (table is None) or (len(table) == 0):
            results.append({
                'name': name,
                'ra': bh_ra,
                'dec': bh_dec,
                'star_ra': None,
                'star_dec': None,
                'sep_arcsec': None,
                'k_mag': None,
                'k_mag_err': None,
                'ph_qual': None,
                'cc_flg': None
            })
            continue

        # Find the brightest star for each black hole 
        best, sep_arcsec = find_brightest_IRSA(table, bh_coord)
        if best is None:
            results.append({
                'name': name,
                'ra': bh_ra,
                'dec': bh_dec,
                'star_ra': None,
                'star_dec': None,
                'sep_arcsec': None,
                'k_mag': None,
                'k_mag_err': None,
                'ph_qual': None,
                'cc_flg': None
            })
            continue

        star_ra_deg = float(best['ra']) if best['ra'] else None
        star_dec_deg = float(best['dec']) if best['dec'] else None

        # Converting the star's coordinates (RA,DEC) to sexagesimal for saving 
        if (star_ra_deg is not None) and (star_dec_deg is not None):
            star_coord = SkyCoord(star_ra_deg, star_dec_deg, unit=('deg','deg'), frame='icrs')
            star_ra = star_coord.ra.to_string(unit=u.hour, sep=":", precision=2)   
            star_dec = star_coord.dec.to_string(unit=u.deg, sep=":", precision=2, alwayssign=True)  
        else:
            star_ra, star_dec = None, None

        results.append({
            'bh_name': name,
            'bh_ra': bh_ra,
            'bh_dec': bh_dec,
            'star_ra': star_ra,
            'star_dec': star_dec,
            'sep_arcsec': float(sep_arcsec) if sep_arcsec else None,
            'k_mag': float(best['k_m']) if best['k_m'] else None,
            'k_mag_err': float(best['k_msigcom']) if best['k_msigcom'] else None,
            'ph_qual': best['ph_qual'][2],
            'cc_flg': best['cc_flg']
        })
    return results



def main():
    input_file = 'results/black_holes'
    output_file = 'results/black_holes_with_2mass.csv'

    black_holes = read_csv(input_file)
    results = process_input_file(black_holes, radius_arcsec=30.0)

    # Specific for black_holes
    out_header = [
        'bh_name',
        'bh_ra','bh_dec',
        'star_ra','star_dec',
        'sep_arcsec','k_mag','k_mag_err','ph_qual','cc_flg'
    ]

    save_results_to_csv(results, output_file, out_header)
    print(f"Saved results to {output_file}")

if __name__ == "__main__":
    main()
