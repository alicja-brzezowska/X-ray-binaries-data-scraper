from astroquery.ipac.irsa import Irsa
from astropy.coordinates import SkyCoord
import astropy.units as u
import csv
import numpy as np

# Preliminary: extracting the list of catalogs with 2mass data 
# "2mass" is not in catalog identifiers, but in catalog descriptions 
# print({k for k, v in Irsa.list_catalogs().items() if "2mass" in k.lower() or "2mass" in v.lower()})

def read_csv(file_name):
    with open(file_name, 'r') as f:
        csv_reader = csv.reader(f)
        header = next(csv_reader)
        return list(csv_reader)

def query_IRSA_catalog(bh_coord, radius_arcsec=30.0, catalog="fp_psc"):
    radius = (radius_arcsec / 3600.0) * u.deg
    return Irsa.query_region(bh_coord, radius=radius, catalog=catalog, spatial='Cone')

def find_brightest_IRSA(table, bh_coord):
    cols = ["ra", "dec", "k_m", "k_msigcom", "ph_qual", "cc_flg"]
    if (table is None) or (len(table) == 0):
        return None, None
    tab = table[cols].filled(np.nan)
    tab.sort('k_m')
    best = tab[0]
    k_object_coord = SkyCoord(best['ra'], best['dec'], unit='deg', frame='icrs')
    sep_arcsec = bh_coord.separation(k_object_coord).arcsec
    return best, sep_arcsec

def sexagesimal_strings(coord):
    s = coord.to_string('hmsdms', sep=':', precision=2, alwayssign=True)
    ra_hms, dec_dms = s.split()
    return ra_hms, dec_dms

def process_black_holes(black_holes, radius_arcsec=30.0):
    results = []
    for bh in black_holes:
        name, ra_in, dec_in = bh[0], bh[1], bh[2]
        bh_coord = SkyCoord(ra_in, dec_in, unit=(u.hourangle, u.deg), frame='icrs')
        ra_deg = bh_coord.ra.deg
        dec_deg = bh_coord.dec.deg
        bh_ra_hms, bh_dec_dms = sexagesimal_strings(bh_coord)

        table = query_IRSA_catalog(bh_coord, radius_arcsec=radius_arcsec, catalog="fp_psc")
        if (table is None) or (len(table) == 0):
            results.append({
                'bh_name': name,
                'bh_ra_hms': bh_ra_hms,
                'bh_dec_dms': bh_dec_dms,
                'k_object_ra_hms': None,
                'k_object_dec_dms': None,
                'sep_arcsec': None,
                'k_mag': None,
                'k_mag_err': None,
                'ph_qual': None,
                'cc_flg': None
            })
            continue

        best, sep_arcsec = find_brightest_IRSA(table, bh_coord)
        if best is None:
            results.append({
                'bh_name': name,
                'bh_ra_hms': bh_ra_hms,
                'bh_dec_dms': bh_dec_dms,
                'k_object_ra_hms': None,
                'k_object_dec_dms': None,
                'sep_arcsec': None,
                'k_mag': None,
                'k_mag_err': None,
                'ph_qual': None,
                'cc_flg': None
            })
            continue

        k_ra_deg = float(best['ra']) if best['ra'] else None
        k_dec_deg = float(best['dec']) if best['dec'] else None
        if (k_ra_deg is not None) and (k_dec_deg is not None):
            k_coord = SkyCoord(k_ra_deg, k_dec_deg, unit=('deg','deg'), frame='icrs')
            k_ra_hms, k_dec_dms = sexagesimal_strings(k_coord)
        else:
            k_ra_hms, k_dec_dms = None, None

        results.append({
            'bh_name': name,
            'bh_ra_hms': bh_ra_hms,
            'bh_dec_dms': bh_dec_dms,
            'k_object_ra_hms': k_ra_hms,
            'k_object_dec_dms': k_dec_dms,
            'sep_arcsec': float(sep_arcsec) if sep_arcsec else None,
            'k_mag': float(best['k_m']) if best['k_m'] else None,
            'k_mag_err': float(best['k_msigcom']) if best['k_msigcom'] else None,
            'ph_qual': best['ph_qual'][2],
            'cc_flg': best['cc_flg']
        })
    return results

def write_results_csv(results, out_path):
    out_header = [
        'bh_name',
        'bh_ra_hms','bh_dec_dms',
        'k_object_ra_hms','k_object_dec_dms',
        'sep_arcsec','k_mag','k_mag_err','ph_qual','cc_flg'
    ]
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=out_header)
        writer.writeheader()
        for r in results:
            writer.writerow(r)

def main():
    black_holes = read_csv('black_holes.csv')
    results = process_black_holes(black_holes, radius_arcsec=30.0)
    write_results_csv(results, 'black_holes_with_2mass.csv')
    print("Saved results to black_holes_with_2mass.csv")

if __name__ == "__main__":
    main()
