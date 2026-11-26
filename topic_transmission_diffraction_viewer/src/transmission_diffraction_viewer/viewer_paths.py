"""Locate viewer-owned files and optional result sources."""

from pathlib import Path


TOPIC = Path(__file__).resolve().parents[2]
INPUT_DIR = TOPIC / 'Input'
OUTPUT_DIR = TOPIC / 'Output'
RESULT_ROOTS = (INPUT_DIR, TOPIC.parent / 'topic_udc_imagingV2' / 'Output')
