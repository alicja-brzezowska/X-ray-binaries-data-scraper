import requests
import json
import numpy as np
import csv
from pprint import pprint
from astropy.coordinates import SkyCoord
import astropy.units as u
from get_GRAVITY_stars import degrees_to_sexagesimal_dec, degrees_to_sexagesimal_ra


def read_csv(file_name):
    """Reads a CSV file and returns a list of black hole data."""
    with open(file_name, 'r') as f:
        csv_reader = csv.reader(f)
        header = next(csv_reader)
        black_holes = list(csv_reader)
    return black_holes


def query_vvv_catalog(ra, dec, radius_arcsec = 30.0):
    radius_deg = 30/ 3600
    """Queries the VVV catalog for sources within a given radius."""
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


def find_brightest_vvv(data, bh_coord):
    """Finds the brightest source in the VVV data."""
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

    # Determine the index of the brightest star from both sets of measurements
    # if both indices show the same source, take the average of their magnitudes
    # if they are different, choose the one with the lower magnitude (brighter)

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


def save_results_to_csv(results, output_file):
    """Saves the results to a CSV file."""
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
    output_file = 'black_holes_with_vvv.csv'
    radius_arcsec = 30.0

    black_holes = read_csv(input_file)
    results = []

    for bh in black_holes:
        name = bh[0]
        ra_str = bh[1]
        dec_str = bh[2]

        bh_coord = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg), frame='icrs')
        ra = bh_coord.ra.deg
        dec = bh_coord.dec.deg

        data = query_vvv_catalog(ra, dec, radius_arcsec)
        if not data:
            print(f"No VVV data found for {name} at RA: {ra_str}, Dec: {dec_str}")
            continue

        brightest = find_brightest_vvv(data, bh_coord)
        if brightest:
           # vvv_coord = SkyCoord(ra * u.deg , dec * u.deg, frame='icrs')
            results.append([
                name, ra_str, dec_str,
                degrees_to_sexagesimal_ra(brightest["vvv_ra"]), degrees_to_sexagesimal_dec(brightest["vvv_dec"]), brightest["sep_arcsec"],
                brightest["ks1"], brightest["ks1err"], brightest["ks2"], brightest["ks2err"],
                brightest["best_mag"]
            ])
    save_results_to_csv(results, output_file)


if __name__ == "__main__":
    main()