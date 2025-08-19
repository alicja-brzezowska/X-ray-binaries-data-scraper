import requests
import json
import numpy as np
import csv
import matplotlib.pyplot as plt
from pprint import pprint
from astropy.coordinates import SkyCoord
import astropy.units as u
from get_GRAVITY_stars import degrees_to_sexagesimal_dec, degrees_to_sexagesimal_ra


def read_csv(file_name):
    with open(file_name, 'r') as f:
        csv_reader = csv.reader(f)
        header = next(csv_reader)
        black_holes = list(csv_reader)
    return black_holes


def query_vvv_catalog(ra, dec, radius_arcsec = 30.0):
    """Query the VVV catalog for K band data around a given RA and Dec."""
    radius_deg = radius_arcsec/ 3600
    query = (
        "SELECT RA2000, DEC2000, KS_1APERMAG3, KS_1APERMAG3ERR, "
        "KS_2APERMAG3, KS_2APERMAG3ERR "
        "FROM VVV_bandMergedSourceCat_V3 "
        f"WHERE CONTAINS(POINT('ICRS', RA2000, DEC2000), "
        f"CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, {radius_deg:.6f}))=1"
    )
    url = f"https://archive.eso.org/tap_cat/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=json&QUERY={query}"
    response = requests.get(url)
    result = json.loads(response.content)
    return result.get('data', [])

def query_vvv_new_catalog(ra, dec, radius_arcsec = 30.0):
    """Query the new VVV catalog for K band data around a given RA and Dec."""
    radius_deg = radius_arcsec/ 3600
    query = (
        "SELECT ra, de, phot_ks_mean_mag, phot_ks_std_mag "
        "FROM VVVX_VIRAC_V2_SOURCES "
        f"WHERE CONTAINS(POINT('ICRS', ra, de), "
        f"CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, {radius_deg:.6f}))=1"
    )
    url = f"https://archive.eso.org/tap_cat/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=json&QUERY={query}"
    response = requests.get(url)
    result = json.loads(response.content)
    return result.get('data', [])


def find_brightest_vvv(data, bh_coord):
    if not data:
        return None

    ra_vals = np.array([row[0] for row in data], dtype=float)
    dec_vals = np.array([row[1] for row in data], dtype=float)
    ks1_vals = np.array([row[2] for row in data], dtype=float)
    ks1_errs = np.array([row[3] for row in data], dtype=float)
    ks2_vals = np.array([row[4] for row in data], dtype=float)
    ks2_errs = np.array([row[5] for row in data], dtype=float)

    # Calculate separations using SkyCoord
    src_coords = SkyCoord(ra_vals * u.deg, dec_vals * u.deg, frame='icrs')
    seps = src_coords.separation(bh_coord).arcsec

    # Find the brightest source
    idx1 = int(np.nanargmin(ks1_vals)) if np.any(~np.isnan(ks1_vals)) else None
    idx2 = int(np.nanargmin(ks2_vals)) if np.any(~np.isnan(ks2_vals)) else None

    chosen_idx = None
    best_mag = np.nan

    """Determine the index of the brightest star from both sets of measurements
    if both indices show the same source, take the average of their magnitudes
    if they are different, choose the one with the lower magnitude (brighter)"""

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
            "vvv_ra": ra_vals[chosen_idx],
            "vvv_dec": dec_vals[chosen_idx],
            "sep_arcsec": seps[chosen_idx],
            "ks1": ks1_vals[chosen_idx],
            "ks1err": ks1_errs[chosen_idx],
            "ks2": ks2_vals[chosen_idx],
            "ks2err": ks2_errs[chosen_idx],
            "best_mag": best_mag,
        }
    return None

def find_brightest_vvv_new(data, bh_coord):
    if not data:
        return None

    ra_vals = np.array([row[0] for row in data], dtype=float)
    dec_vals = np.array([row[1] for row in data], dtype=float)
    ks_vals = np.array([row[2] for row in data], dtype=float)
    ks_errs = np.array([row[3] for row in data], dtype=float)

    # Calculate separations using SkyCoord
    src_coords = SkyCoord(ra_vals * u.deg, dec_vals * u.deg, frame='icrs')
    seps = src_coords.separation(bh_coord).arcsec

    # In the new catalog, only one K band magnitude 
    idx = int(np.nanargmin(ks_vals)) if np.any(~np.isnan(ks_vals)) else None

    if idx is not None:
        return {
            "vvv_ra": ra_vals[idx],
            "vvv_dec": dec_vals[idx],
            "sep_arcsec": seps[idx],
            "ks": ks_vals[idx],
            "kserr": ks_errs[idx],
            "best_mag": ks_vals[idx]
        }
    return None


def rect_to_radec(lmin, lmax, bmin, bmax, n=50):
    """Convert rectangular coordinates (l,b) to RA, Dec in degrees."""
    l_vals = [lmin] * n + list(np.linspace(lmin, lmax, n)) + [lmax] * n + list(np.linspace(lmax, lmin, n))
    b_vals = list(np.linspace(bmin, bmax, n)) + [bmax] * n + list(np.linspace(bmax, bmin, n)) + [bmin] * n
    gal = SkyCoord(l=l_vals*u.deg, b=b_vals*u.deg, frame="galactic")
    return gal.icrs.ra.deg, gal.icrs.dec.deg


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


def save_results_to_csv(results, output_file):
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'name', 'bh_ra', 'bh_dec',
            'vvv_ra', 'vvv_dec', 'sep_arcsec',
            'ks1_mag', 'ks1_mag_err', 'ks2_mag', 'ks2_mag_err',
            'best_mag_for_ranking'
        ])
        writer.writerows(results)
    print(f"Saved results to {output_file}")


def main():
    input_file = 'black_holes.csv'
    output_file = '../results/black_holes_with_vvv.csv'
    output_file_new = '../results/black_holes_with_vvv_new.csv'
    radius_arcsec = 30.0

    black_holes = read_csv(input_file)
    
    results = []  
    results_new = []
    bh_ra_list, bh_dec_list = [], []

    for bh in black_holes:
        name = bh[0]
        ra_str = bh[1]
        dec_str = bh[2]

        bh_coord = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg), frame='icrs')
        ra = bh_coord.ra.deg
        dec = bh_coord.dec.deg

        bh_ra_list.append(ra)
        bh_dec_list.append(dec)

        data = query_vvv_catalog(ra, dec, radius_arcsec)
        if not data:
            print(f"No VVV data found for {name} at RA: {ra_str}, Dec: {dec_str}")
            continue

        brightest = find_brightest_vvv(data, bh_coord)
        if brightest:
            results.append([
                name, ra_str, dec_str,
                degrees_to_sexagesimal_ra(brightest["vvv_ra"]), degrees_to_sexagesimal_dec(brightest["vvv_dec"]), brightest["sep_arcsec"],
                brightest["ks1"], brightest["ks1err"], brightest["ks2"], brightest["ks2err"],
                brightest["best_mag"]
            ])

        data_new = query_vvv_new_catalog(ra, dec, radius_arcsec)
        if not data_new:
            print(f"No new VVV data found for {name} at RA: {ra_str}, Dec: {dec_str}")
            continue

        brightest_new = find_brightest_vvv_new(data_new, bh_coord)
        if brightest_new:
            results_new.append([
                name, ra_str, dec_str,
                degrees_to_sexagesimal_ra(brightest_new["vvv_ra"]), degrees_to_sexagesimal_dec(brightest_new["vvv_dec"]), brightest_new["sep_arcsec"],
                brightest_new["ks"], brightest_new["kserr"], "", "",
                brightest_new["best_mag"]
            ])


    save_results_to_csv(results, output_file)
    save_results_to_csv(results_new, output_file_new)

    visualize_results(bh_ra_list, bh_dec_list)

    

if __name__ == "__main__":
    main()