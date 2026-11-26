"""Read UDC metrics for the static benchmark reports."""

import csv
import numpy as np


DIRECTIONS = ('Horizontal', 'Vertical', 'Diagonal45', 'Diagonal135')
PSF_FLOOR_DB = -70
PSF_COLORMAP = 'magma'


def read_csv(folder, name):
    with (folder / name).open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


def green_mtf(profiles):
    paths = [[r for r in profiles if r['Channel'] == 'G' and r['Direction'] == direction]
             for direction in DIRECTIONS]
    frequencies = np.array([[float(r['CyclesPerMm']) for r in path] for path in paths])
    for frequency in frequencies[1:]:
        np.testing.assert_array_equal(frequencies[0], frequency)
    mtf = np.array([[float(r['MTF']) for r in path] for path in paths])
    clear = np.array([[float(r['ClearMTF']) for r in path] for path in paths])
    return frequencies[0], mtf.min(axis=0), clear.min(axis=0)
