import requests
import json
import numpy as np
import csv
from pprint import pprint
from astropy.coordinates import SkyCoord
import astropy.units as u

radius = 30.0 / 3600.0  # 30 arcsec in deg

results = []

with open('black_holes.csv', 'r') as f:
    csv_reader = csv.reader(f)
    header = next(csv_reader)
    black_holes = list(csv_reader)


for bh in black_holes:
    name = bh[0]
    ra_str = bh[1]
    dec_str = bh[2]

    # Convert black hole coordinates to decimal degrees
    bh_coord = SkyCoord(ra_str, dec_str, unit=(
        u.hourangle, u.deg), frame='icrs')
    ra = bh_coord.ra.deg
    dec = bh_coord.dec.deg

    query = (
        "SELECT RA2000, DEC2000, KS_1APERMAG3, KS_1APERMAG3ERR, "
        "KS_2APERMAG3, KS_2APERMAG3ERR "
        "FROM VVV_bandMergedSourceCat_V3 "
        f"WHERE CONTAINS(POINT('ICRS', RA2000, DEC2000), "
        f"CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, {radius:.6f}))=1"
    )

    url = f"https://archive.eso.org/tap_cat/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=json&QUERY={query}"

    r = requests.get(url)
    result = json.loads(r.content)
    print(result)
    data = result.get('data', [])

    vvv_ra = vvv_dec = sep_arcsec = ks1 = ks1err = ks2 = ks2err = best_mag = ""

    if not data:
        print(f"No VVV data found for {name} at RA: {ra_str}, Dec: {dec_str}")
        continue

    pprint(data)

    ra_vals = np.array([(row[0]) for row in data], dtype=float)
    dec_vals = np.array([(row[1]) for row in data], dtype=float)
    ks1_vals = np.array([(row[2]) for row in data], dtype=float)
    ks1_errs = np.array([(row[3]) for row in data], dtype=float)
    ks2_vals = np.array([(row[4]) for row in data], dtype=float)
    ks2_errs = np.array([(row[5]) for row in data], dtype=float)

    # Using SkyCoord to calculate separations
    src_coords = SkyCoord(ra_vals * u.deg, dec_vals * u.deg, frame='icrs')
    seps = src_coords.separation(bh_coord).arcsec

    # Calculating the brightest source, taking the average if both measurements match
    idx1 = None
    if np.any(~np.isnan(ks1_vals)):
        idx1 = int(np.nanargmin(ks1_vals))

    idx2 = None
    if np.any(~np.isnan(ks2_vals)):
        idx2 = int(np.nanargmin(ks2_vals))

    chosen_idx = None
    best_mag = np.nan

    if (idx1 is not None and idx2 is not None):
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
        vvv_ra = ra_vals[chosen_idx]
        vvv_dec = dec_vals[chosen_idx]
        sep_arcsec = seps[chosen_idx]
        ks1 = ks1_vals[chosen_idx]
        ks1err = ks1_errs[chosen_idx]
        ks2 = ks2_vals[chosen_idx]
        ks2err = ks2_errs[chosen_idx]
        best_mag = best_mag

    results.append([
        name, ra, dec, vvv_ra, vvv_dec, sep_arcsec,
        ks1, ks1err, ks2, ks2err, best_mag
    ])

with open('black_holes_with_vvv.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow([
        'name', 'bh_ra', 'bh_dec',
        'vvv_ra', 'vvv_dec', 'sep_arcsec',
        'ks1_mag', 'ks1_mag_err', 'ks2_mag', 'ks2_mag_err',
        'best_mag_for_ranking'
    ])
    writer.writerows(results)

print("Saved black_holes_with_vvv.csv")
