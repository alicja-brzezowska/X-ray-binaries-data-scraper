import numpy as np
import csv
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u

def get_adist(ra1, dec1, ra2, dec2):
    deg = np.pi/180.0
    ra1 *= deg
    dec1 *= deg
    ra2 *= deg
    dec2 *= deg
    adist = np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra2-ra1)
    adist = np.arccos(np.clip(adist, -1, 1)) * 3600 / deg
    return adist

def degrees_to_sexagesimal_ra(x):
    tmp = x / 15.0
    ra_h = int(tmp)
    tmp = (tmp - ra_h)*60.0
    ra_m = int(tmp)
    ra_s = (tmp-ra_m)*60.0
    return '{:02d}:{:02d}:{:05.2f}'.format(ra_h, ra_m, ra_s)

def degrees_to_sexagesimal_dec(x):
    sign = '+' if x >= 0 else '-'
    tmp = abs(x)
    dec_d = int(tmp)
    tmp = (tmp - dec_d)*60.0
    dec_m = int(tmp)
    dec_s = (tmp-dec_m)*60.0
    return '{}{:02d}:{:02d}:{:05.2f}'.format(sign, dec_d, dec_m, dec_s)


def find_brightest_gaia_rp(ra, dec, radius_arcsec=30.0):
    radius_deg = radius_arcsec / 3600.0
    query = """
        SELECT gaia.source_id, gaia.ra, gaia.dec, gaia.phot_rp_mean_mag
        FROM gaiaedr3.gaia_source AS gaia
        WHERE 1 = CONTAINS(
            POINT('ICRS', {:.6f}, {:.6f}),
            CIRCLE('ICRS', gaia.ra, gaia.dec, {:.6f}))
    """.format(ra, dec, radius_deg)

    job = Gaia.launch_job(query)
    r = job.get_results()
    if len(r) == 0:
        return None

    best_i = None
    best_rp = float('inf')
    best_adist = None

    for i in range(len(r)):
        rp = r['phot_rp_mean_mag'][i]
        if rp is None or (isinstance(rp, float) and np.isnan(rp)):
            continue
        adist = get_adist(ra, dec, r['ra'][i], r['dec'][i])
        if rp < best_rp:
            best_rp = rp
            best_i = i
            best_adist = adist

    if best_i is None:
        return None

    ra_s  = degrees_to_sexagesimal_ra(r['ra'][best_i])
    dec_s = degrees_to_sexagesimal_dec(r['dec'][best_i])

    return {
        'source_id': r['source_id'][best_i],
        'ra': float(r['ra'][best_i]),
        'dec': float(r['dec'][best_i]),
        'ra_s': ra_s,
        'dec_s': dec_s,
        'rp': float(r['phot_rp_mean_mag'][best_i]),
        'sep_arcsec': float(best_adist)
    }

results = []
with open('black_holes.csv', 'r') as f:
    csv_reader = csv.reader(f)
    header = next(csv_reader)
    black_holes = list(csv_reader)

for bh in black_holes:
    name = bh[0]
    ra_str = bh[1]
    dec_str = bh[2]

    bh_coord = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg), frame='icrs')
    ra = bh_coord.ra.deg
    dec = bh_coord.dec.deg

    rec = find_brightest_gaia_rp(ra, dec, radius_arcsec=30.0)
    if rec is None:
        results.append([name, ra, dec, 'NONE', '', '', '', ''])
    else:
        results.append([
            name,
            ra,
            dec,
            rec['source_id'],
            rec['ra'],
            rec['dec'],
            rec['rp'],
            rec['sep_arcsec']
        ])

with open('black_holes_with_gaia.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['bh_name', 'ra_bh', 'dec_bh', 'star_id', 'ra_star', 'dec_star', 'rp_mag_star', 'sep_arcsec'])
    writer.writerows(results)

print("Saved results to black_holes_with_gaia.csv")
