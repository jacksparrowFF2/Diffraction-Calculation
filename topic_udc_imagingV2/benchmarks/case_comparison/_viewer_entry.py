"""Forward legacy comparison commands to the standalone viewer topic."""

from pathlib import Path
from runpy import run_path
import sys


def run_viewer_entry(filename):
    source = Path(__file__).resolve().parents[3] / 'topic_transmission_diffraction_viewer' / 'src'
    if str(source) not in sys.path:
        sys.path.insert(0, str(source))
    run_path(str(source / filename), run_name='__main__')
