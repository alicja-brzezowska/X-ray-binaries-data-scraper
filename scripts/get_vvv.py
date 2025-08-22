import requests
import json
import numpy as np
import csv
import matplotlib.pyplot as plt
from pprint import pprint
from astropy.coordinates import SkyCoord
import astropy.units as u

from functions import read_csv, save_results_to_csv, rect_to_radec


def query_VVV(ra, dec, radius_arcsec = 30.0):
    """Query the VVV catalog for K band data around a given RA and DEC."""
    radius_deg = radius_arcsec/ 3600 
    query = (
        "SELECT RA2000, DEC2000, KS_1APERMAG3, KS_1APERMAG3ERR, "  # Old catalog has two K band brightness measurementse
        "KS_2APERMAG3, KS_2APERMAG3ERR "
        "FROM VVV_bandMergedSourceCat_V3 "
        f"WHERE CONTAINS(POINT('ICRS', RA2000, DEC2000), "
        f"CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, {radius_deg:.6f}))=1"
    )
    url = f"https://archive.eso.org/tap_cat/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=json&QUERY={query}"
    response = requests.get(url)
    result = json.loads(response.content)
    return result.get('data', [])


def query_VVV_new(ra, dec, radius_arcsec = 30.0):
    """Query the new VVV catalog for K band data around a given RA and Dec."""
    radius_deg = radius_arcsec/ 3600
    query = (
        "SELECT ra, de, phot_ks_mean_mag, phot_ks_std_mag "    # New catalog has one K band brightness measurement
        "FROM VVVX_VIRAC_V2_SOURCES "
        f"WHERE CONTAINS(POINT('ICRS', ra, de), "
        f"CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, {radius_deg:.6f}))=1"
    )
    url = f"https://archive.eso.org/tap_cat/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=json&QUERY={query}"
    response = requests.get(url)
    result = json.loads(response.content)
    return result.get('data', [])


def find_brightest_vvv(data, bh_coord):
    """Find the brightest K band source in the queried VVV data."""
    if not data:
        return None

    ra_vals = np.array([row[0] for row in data], dtype=float)
    dec_vals = np.array([row[1] for row in data], dtype=float)
    ks1_vals = np.array([row[2] for row in data], dtype=float)
    ks1_errs = np.array([row[3] for row in data], dtype=float)
    ks2_vals = np.array([row[4] for row in data], dtype=float)
    ks2_errs = np.array([row[5] for row in data], dtype=float)

    # Calculate separations 
    src_coords = SkyCoord(ra_vals * u.deg, dec_vals * u.deg, frame='icrs')
    seps = src_coords.separation(bh_coord).arcsec

    # Find the brightest sources
    idx1 = int(np.nanargmin(ks1_vals)) if np.any(~np.isnan(ks1_vals)) else None
    idx2 = int(np.nanargmin(ks2_vals)) if np.any(~np.isnan(ks2_vals)) else None

    chosen_idx = None
    best_mag = np.nan

    """Determine the index of the brightest star from both sets of measurements.
    If both indices show the same source, take the average of their magnitudes.
    If they are different, choose the one with the lower magnitude (brighter)"""

    if idx1 is not None and idx2 is not None:
        if idx1 == idx2:
            chosen_idx = idx1
            best_mag = (ks1_vals[idx1] + ks2_vals[idx2]) / 2.0
        else:
            if ks1_vals[idx1] < ks2_vals[idx2]:
                chosen_idx = idx1
                best_mag = ks1_vals[idx1]
            else:
                chosen_idx = idx2
                best_mag = ks2_vals[idx2]
    elif idx1 is not None:
        chosen_idx = idx1
        best_mag = ks1_vals[idx1]
    elif idx2 is not None:
        chosen_idx = idx2
        best_mag = ks2_vals[idx2]

    if chosen_idx is not None:
        return {
            "star_ra_deg": ra_vals[chosen_idx],
            "star_dec_deg": dec_vals[chosen_idx],
            "sep_arcsec": seps[chosen_idx],
            "ks1": ks1_vals[chosen_idx],
            "ks1err": ks1_errs[chosen_idx],
            "ks2": ks2_vals[chosen_idx],
            "ks2err": ks2_errs[chosen_idx],
            "best_mag": best_mag,
        }
    return None

def find_brightest_vvv_new(data, bh_coord):
    """Find the brightest K band source in the queried VVV data (new catalog)."""
    if not data:
        return None

    ra_vals = np.array([row[0] for row in data], dtype=float)
    dec_vals = np.array([row[1] for row in data], dtype=float)
    ks_vals = np.array([row[2] for row in data], dtype=float)
    ks_errs = np.array([row[3] for row in data], dtype=float)

    # Calculate separations 
    src_coords = SkyCoord(ra_vals * u.deg, dec_vals * u.deg, frame='icrs')
    seps = src_coords.separation(bh_coord).arcsec

    # Find the brightest source
    idx = int(np.nanargmin(ks_vals)) if np.any(~np.isnan(ks_vals)) else None

    if idx is not None:
        return {
            "star_ra_deg": ra_vals[idx],
            "star_dec_deg": dec_vals[idx],
            "sep_arcsec": seps[idx],
            "ks": ks_vals[idx],
            "kserr": ks_errs[idx]
        }
    return None


def visualize_results(bh_ra, bh_dec):
    """Plot the VVV survey range and black hole locations to double-check the results."""
    ra_bulge, dec_bulge = rect_to_radec(-10, +10, -10, +5)
    ra_disk, dec_disk = rect_to_radec(-65.0, -10.0, -2, +2)

    plt.figure(figsize=(10, 7))
    plt.scatter(bh_ra, bh_dec, color="red", label="Black Holes", zorder=3)
    plt.plot(ra_bulge, dec_bulge, "b-", label="VVV Bulge", zorder=2)
    plt.plot(ra_disk, dec_disk, "g-", label="VVV Disk", zorder=2)

    plt.gca().invert_xaxis()  
    plt.xlabel("RA (deg)")
    plt.ylabel("Dec (deg)")
    plt.title("VVV Survey Range and X-ray Binary Locations")
    plt.legend()
    return plt.show()


def main():
    input_file = 'results/black_holes.csv'
    output_file = 'results/black_holes_with_vvv.csv'
    output_file_new = 'results/black_holes_with_vvv_new.csv'
    radius_arcsec = 30.0


    # Extract RA and DEC data from the black_holes file 
    black_holes = read_csv(input_file)
    

    results = []  
    results_new = []
    bh_ra_list, bh_dec_list = [], [] #for plotting purposes

    for bh in black_holes:
        name = bh[0]
        bh_ra = bh[1]
        bh_dec = bh[2]

        # Convert the coordinates to degrees for queries
        bh_coord = SkyCoord(bh_ra, bh_dec, unit=(u.hourangle, u.deg), frame='icrs')
        bh_ra_deg = bh_coord.ra.deg
        bh_dec_deg = bh_coord.dec.deg

        bh_ra_list.append(bh_ra_deg)
        bh_dec_list.append(bh_dec_deg)

        # Perform the search for old and new VVV catalogs
        data = query_VVV(bh_ra_deg, bh_dec_deg, radius_arcsec)
        if not data:
            print(f"No VVV data found for {name} at RA: {bh_ra}, Dec: {bh_dec}")
            results.append([
                name, bh_ra, bh_dec,
                "", "", "", "", "", "", "", ""
        ])
        else:
            brightest = find_brightest_vvv(data, bh_coord)
            if brightest:
                k = SkyCoord(brightest["star_ra_deg"] * u.deg, brightest["star_dec_deg"] * u.deg, frame='icrs')
                star_ra = k.ra.to_string(unit=u.hour, sep=':', precision=2, pad=True)
                star_dec = k.dec.to_string(unit=u.deg,  sep=':', precision=2, pad=True, alwayssign=True)
                results.append([
                    name, bh_ra, bh_dec,
                    star_ra, star_dec,
                    brightest["sep_arcsec"],
                    brightest["ks1"], brightest["ks1err"], brightest["ks2"], brightest["ks2err"],
                    brightest["best_mag"]
                ])
            else:
                results.append([
                    name, bh_ra, bh_dec,
                    "", "", "", "", "", "", "", ""
                ])

        brightest = find_brightest_vvv(data, bh_coord)

        if brightest:
            k = SkyCoord(brightest["star_ra_deg"] * u.deg, brightest["star_dec_deg"] * u.deg, frame='icrs')
            star_ra = k.ra.to_string(unit=u.hour, sep=':', precision=2, pad=True)
            star_dec = k.dec.to_string(unit=u.deg,  sep=':', precision=2, pad=True, alwayssign=True)
            results.append([
                name, bh_ra, bh_dec,
                star_ra, star_dec,
                brightest["sep_arcsec"],
                brightest["ks1"], brightest["ks1err"], brightest["ks2"], brightest["ks2err"],
                brightest["best_mag"]
            ])

                
        data_new = query_VVV_new(bh_ra_deg, bh_dec_deg, radius_arcsec)

        if not data_new:
            print(f"No new VVV data found for {name} at RA: {bh_ra}, Dec: {bh_dec}")
            results_new.append([
                name, bh_ra, bh_dec,
                "", "", "", "", ""
         ])
        else:
            brightest_new = find_brightest_vvv_new(data_new, bh_coord)
            if brightest_new:
                k2 = SkyCoord(brightest_new["star_ra_deg"] * u.deg, brightest_new["star_dec_deg"] * u.deg, frame='icrs')
                vvv2_ra_hms = k2.ra.to_string(unit=u.hour, sep=':', precision=2, pad=True)
                vvv2_dec_dms = k2.dec.to_string(unit=u.deg,  sep=':', precision=2, pad=True, alwayssign=True)
                results_new.append([
                    name, bh_ra, bh_dec,
                    vvv2_ra_hms, vvv2_dec_dms,
                    brightest_new["sep_arcsec"],
                    brightest_new["ks"], brightest_new["kserr"]
                ])
            else:
                results_new.append([
                    name, bh_ra, bh_dec,
                    "", "", "", "", ""
                ])


    
    # Header names for both output files
    header_old = [
        'name', 'bh_ra', 'bh_dec',
        'star_ra', 'star_dec', 'sep_arcsec',
        'ks_1_mag', 'ks_1_mag_err', 'ks_2_mag', 'ks_2_mag_err',
        'best_mag_for_ranking'
    ]

    header_new = [
        'name', 'bh_ra', 'bh_dec',
        'star_ra', 'star_dec', 'sep_arcsec',
        'ks_mag', 'ks_mag_err'
    ]

    save_results_to_csv(results, output_file, header=header_old)
    save_results_to_csv(results_new, output_file_new, header=header_new) 

    visualize_results(bh_ra_list, bh_dec_list)

    
if __name__ == "__main__":
    main()