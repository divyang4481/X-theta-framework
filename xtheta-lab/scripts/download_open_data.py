#!/usr/bin/env python3
"""
Helper script to download and extract open Bell/CHSH datasets.
"""
from __future__ import annotations

import argparse
import os
import urllib.request
import zipfile
from pathlib import Path

HENSEN_ZIP_URL = "https://data.4tu.nl/file/e8cf2991-3153-48ad-b67d-dfd7d7d97fd3/289c4850-6ed5-45a3-8b19-de8671f873a8"

def download_file(url: str, dest: Path):
    print(f"Downloading {url} to {dest}...")
    urllib.request.urlretrieve(url, dest)
    print("Download complete.")

def extract_zip(zip_path: Path, extract_to: Path):
    print(f"Extracting {zip_path} to {extract_to}...")
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_to)
    print("Extraction complete.")

def main():
    parser = argparse.ArgumentParser(description="Download open Bell/CHSH datasets.")
    parser.add_argument("--dataset", choices=["hensen"], required=True, help="Dataset to download.")
    parser.add_argument("--force", action="store_true", help="Force re-download and re-extraction.")

    args = parser.parse_args()

    if args.dataset == "hensen":
        base_dir = Path("data/open_bell/hensen")
        raw_dir = base_dir / "raw"
        zip_path = base_dir / "data.zip"
        target_file = raw_dir / "bell_open_data.txt"

        base_dir.mkdir(parents=True, exist_ok=True)

        # Download
        if not zip_path.exists() or args.force:
            download_file(HENSEN_ZIP_URL, zip_path)
        else:
            print(f"Skip download: {zip_path} already exists. Use --force to re-download.")

        # Extract
        if not target_file.exists() or args.force:
            raw_dir.mkdir(parents=True, exist_ok=True)
            extract_zip(zip_path, raw_dir)
        else:
            print(f"Skip extraction: {target_file} already exists. Use --force to re-extract.")

        # Final check
        if target_file.exists():
            print(f"\nSuccess! Hensen data ready at: {target_file}")
            print("\nValidation command:")
            print(f"python scripts/run_open_data_validation.py --dataset hensen --data {target_file} --output outputs/open_data/hensen --bootstrap-samples 1000")
        else:
            print(f"\nError: Expected file not found after extraction: {target_file}")

if __name__ == "__main__":
    main()
