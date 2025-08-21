import numpy as np
import csv
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u

from functions import read_csv, save_results_to_csv


def query_Gaia(ra, dec, radius_arcsec = 30.0):
    "Query the Gaia catalog for objects around a given RA and DEC."
    radius_deg = radius_arcsec / 3600.0
    query = """
        SELECT gaia.ra, gaia.dec, gaia.phot_rp_mean_mag
        FROM gaiaedr3.gaia_source AS gaia
        WHERE 1 = CONTAINS(
            POINT('ICRS', {:.6f}, {:.6f}),
            CIRCLE('ICRS', gaia.ra, gaia.dec, {:.6f}))
    """.format(ra, dec, radius_deg)     
    job = Gaia.launch_job(query)
    r = job.get_results()
    return r

def find_brightest_gaia(data, bh_coord):
    "Given the queried data, find the brightest star (in the rp filter) for each input object."
    "Return the coordinates and brightness in the Rp filter"

    if not data:
        return None 
    
    star_ra_deg = np.array([row[0] for row in data], dtype = float)
    star_dec_deg = np.array([row[1] for row in data], dtype = float)
    rp_mean_mag = np.array([row[2] for row in data], dtype = float)

    # Calculate separations 
    src_coords = SkyCoord(star_ra_deg * u.deg, star_dec_deg * u.deg, frame='icrs')
    seps = src_coords.separation(bh_coord).arcsec

    # Find the brightest source
    idx = int(np.nanargmin(rp_mean_mag)) if np.any(~np.isnan(rp_mean_mag)) else None

    if idx is not None:
        return {
            "star_ra_deg": star_ra_deg[idx],
            "star_dec_deg": star_dec_deg[idx],
            "sep_arcsec": seps[idx],
            "rp_mean_mag": rp_mean_mag[idx],
        }
    return None


def main():
    input_file = "results/black_holes.csv"
    output_file = "results/black_holes_with_gaia.csv"
    radius_arcsec = 30.0 

    black_holes = read_csv(input_file)

    results = []

    for bh in black_holes:
        name = bh[0]
        bh_ra = bh[1]  
        bh_dec = bh[2]  

        # Convert RA and Dec from sexagesimal to decimal degrees to pass into Gaia query
        bh_coord = SkyCoord(bh_ra, bh_dec, unit=(u.hourangle, u.deg), frame='icrs')
        ra = bh_coord.ra.deg
        dec = bh_coord.dec.deg

        # Perform the search for old and new VVV catalogs
        data = query_Gaia(ra, dec, radius_arcsec)
        if not data:
            print(f"No VVV data found for {name} at RA: {bh_ra}, Dec: {bh_dec}")
            continue

        brightest = find_brightest_gaia(data, bh_coord)
        if brightest is None:
            results.append([
                name,
                bh_ra,  
                bh_dec,  
                  '', '', '', ''
            ])
        else:
            k = SkyCoord(brightest["star_ra_deg"] * u.deg, brightest["star_dec_deg"] * u.deg, frame='icrs')
            star_ra = k.ra.to_string(unit=u.hour, sep=':', precision=2, pad=True)
            star_dec = k.dec.to_string(unit=u.deg,  sep=':', precision=2, pad=True, alwayssign=True) 
            results.append([
                name,
                bh_ra,  
                bh_dec, 
                star_ra,  
                star_dec, 
                brightest['sep_arcsec'],
                brightest['rp_mean_mag']
            ])
 
    output_header = [
        'bh_name', 'bh_ra', 'bh_dec',
        'star_ra', 'star_dec', 'sep_arcsec',
        'rp_mag'
    ]

    save_results_to_csv(results, output_file, output_header)


if __name__ == "__main__":
    main()
