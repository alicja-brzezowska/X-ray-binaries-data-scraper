import csv 
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import os


def read_csv(file_name):
    """Read a CSV file"""
    with open(file_name, 'r') as f:
        csv_reader = csv.reader(f)
        header = next(csv_reader)
        return list(csv_reader)
    

def save_results_to_csv(rows, output_file, header):
    """Save rows to CSV. Works for both lists and dictonaries."""
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    rows = list(rows or [])

    with open(output_file, "w", newline="", encoding="utf-8") as f:
        # If first row is a dict, use DictWriter
        if rows and isinstance(rows[0], dict):
            w = csv.DictWriter(f, fieldnames=header)
            w.writeheader()
            for r in rows:
                w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in header})
        else:
            # Accomodate for lists
            w = csv.writer(f)
            w.writerow(header)
            w.writerows(rows)


def rect_to_radec(lmin, lmax, bmin, bmax, n=50):
    """Convert rectangular coordinates (l,b) to RA, Dec in degrees."""
    l_vals = [lmin] * n + list(np.linspace(lmin, lmax, n)) + [lmax] * n + list(np.linspace(lmax, lmin, n))
    b_vals = list(np.linspace(bmin, bmax, n)) + [bmax] * n + list(np.linspace(bmax, bmin, n)) + [bmin] * n
    gal = SkyCoord(l=l_vals*u.deg, b=b_vals*u.deg, frame="galactic")
    return gal.icrs.ra.deg, gal.icrs.dec.deg