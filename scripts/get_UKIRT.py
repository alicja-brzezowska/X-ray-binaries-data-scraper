import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import astropy.units as u
from astroquery.ukidss import Ukidss



from functions import read_csv, save_results_to_csv, rect_to_radec


def query_UKIRT(ra_deg, dec_deg, radius_arcsec= 30.0):
    """Query the UKIRT catalog for objects in 30'' around a specified coordinate"""
    center = SkyCoord(ra_deg, dec_deg, unit=(u.deg, u.deg), frame='icrs')
    radius = radius_arcsec * u.arcsec

    table = Ukidss.query_region(
        center,
        radius=radius,
        programme_id='GPS',
        database='UKIDSSDR11PLUS',  # selected database
    )
    return table


def find_brightest_UKIRT(table, bh_coord):
    """
    Given an astroquery.ukidss table and black_hole coordinates 
    return a table of brightest objects around each with 'ra','dec','kAperMag3','kAperMag3Err', 
    and separation columns
    """
    cols = ["ra", "dec", "kAperMag3", "kAperMag3Err"]
    if (table is None) or (len(table) == 0):
        return None, None

    tab = table[cols].filled(np.nan)

    # Recording only brightness within a reasonable range
    bad_mag = (~np.isfinite(tab['kAperMag3']) |
               (tab['kAperMag3'] < -50) | (tab['kAperMag3'] > 50))
    bad_err = (~np.isfinite(tab['kAperMag3Err']) |
               (tab['kAperMag3Err'] <= 0) | (tab['kAperMag3Err'] > 1.5))
    good = ~(bad_mag | bad_err)
    tab = tab[good]

    if len(tab) == 0:
        return None, None

    tab.sort('kAperMag3')     
    best = tab[0]

    k_coord = SkyCoord(float(best['ra']), float(best['dec']),
                       unit='deg', frame='icrs')
    sep_arcsec = bh_coord.separation(k_coord).arcsec
    return best, sep_arcsec


def visualize_results(bh_ra, bh_dec):
    "Visualize the range of the UKIDSS GPS survey alongside the black holes."

    # Range of the UKIDSS Galactic Plane Survey
    ra_1, dec_1 = rect_to_radec(141, 230, -5, 5)
    ra_2, dec_2 = rect_to_radec(15, 107, -5, 5)
    ra_3, dec_3 = rect_to_radec(-2,15, -2,2)

    plt.figure(figsize=(10, 7))
    plt.scatter(bh_ra, bh_dec, color="red", label="Black Holes", zorder=3)
    plt.plot(ra_1, dec_1, "b-", label="Survey range", zorder=2)
    plt.plot(ra_2, dec_2, "b-", zorder=2)
    plt.plot(ra_3, dec_3, "b-", zorder=2)

    plt.gca().invert_xaxis()  
    plt.xlabel("RA (deg)")
    plt.ylabel("Dec (deg)")
    plt.title("VVV Survey Range and X-ray Binary Locations")
    plt.legend()
    return plt.show()


def main():
    input_file = 'results/black_holes.csv'
    output_file = 'results/black_holes_with_UKIRT.csv'
    radius_arcsec = 30.0

    black_holes = read_csv(input_file)

    results = []
    bh_ra_list, bh_dec_list = [], [] #for plotting purposes

    # Extract RA and DEC data from the black_holes file
    for bh in black_holes:
        name, bh_ra, bh_dec = bh[0], bh[1], bh[2]

        # Transform black hole RA and DEC into degrees for query purposes
        bh_coord = SkyCoord(bh_ra, bh_dec, unit=(u.hourangle, u.deg), frame='icrs')
        bh_ra_list.append(bh_coord.ra.deg)
        bh_dec_list.append(bh_coord.dec.deg)


        # Query data from the UKIRT database 
        try:
            table = query_UKIRT(bh_coord.ra.deg, bh_coord.dec.deg, radius_arcsec=radius_arcsec)
        except Exception:
            table = None

        if (table is None) or (len(table) == 0):
            print(f"[INFO] No UKIDSS GPS sources found within {radius_arcsec}\" of {name}")
            results.append({
                'bh_name': name,
                'bh_ra': bh_ra,
                'bh_dec': bh_dec,
                'star_ra': "-",
                'star_dec': "-",
                'sep_arcsec': "-",
                'k_mag': "-",
                'k_mag_err': "-"
            })
            continue

        best, sep_arcsec = find_brightest_UKIRT(table, bh_coord)
        if not best:
            continue

        # Convert star coordinates into sexagesimal for output
        try:
            star_coord = SkyCoord(float(best['ra']), float(best['dec']), unit=('deg','deg'), frame='icrs')
            star_ra = star_coord.ra.to_string(unit=u.hour, sep=":", precision=2)   
            star_dec = star_coord.dec.to_string(unit=u.deg, sep=":", precision=2, alwayssign=True)  
            
        except Exception: # If nothing is found
            star_ra, star_dec = None, None 

        # Append results
        results.append({
            'bh_name': name,
            'bh_ra': bh_ra,
            'bh_dec': bh_dec,
            'star_ra':star_ra,
            'star_dec': star_dec,
            'sep_arcsec': float(sep_arcsec),
            'k_mag': float(best['kAperMag3']) if best['kAperMag3'] is not None else None,
            'k_mag_err': float(best['kAperMag3Err']) if best['kAperMag3Err'] is not None else None,
        })


    header = [
        'bh_name',
        'bh_ra','bh_dec',
        'star_ra','star_dec',
        'sep_arcsec','k_mag','k_mag_err'
    ]
    # Save results 
    save_results_to_csv(results, output_file, header)
    visualize_results(bh_ra_list, bh_dec_list)

if __name__ == "__main__":
    main()








 