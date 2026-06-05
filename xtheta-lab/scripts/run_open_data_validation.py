#!/usr/bin/env python3
"""
Command-line runner for Open Bell/CHSH Data Validation.
"""
import argparse
import sys
import os
from pathlib import Path

# Ensure xtheta-lab is in path
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from xtheta.experiments.open_data_validation import run_open_data_chsh_validation
from xtheta.data.adapters.weihs import load_weihs_dataset
from xtheta.data.adapters.hensen import load_hensen_dataset
from xtheta.data.adapters.big_bell_test import load_big_bell_test_dataset
from xtheta.data.loaders import get_loader

def main():
    parser = argparse.ArgumentParser(description="Run Open Bell/CHSH Data Validation for X-Theta.")
    parser.add_argument("--dataset", choices=["generic", "weihs", "hensen", "big_bell_test"],
                        default="generic", help="Dataset type to process.")
    parser.add_argument("--data", required=True, help="Path to data file or directory.")
    parser.add_argument("--output", default="outputs", help="Output directory.")
    parser.add_argument("--chunksize", type=int, default=200_000, help="Chunk size for streaming.")
    parser.add_argument("--bootstrap-samples", type=int, default=1000, help="Number of bootstrap samples.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for bootstrap.")
    parser.add_argument("--geometry-fit", choices=["smax-envelope", "xy", "xz"],
                        default="smax-envelope", help="Fit geometry to use.")

    args = parser.parse_args()

    data_path = Path(args.data)
    if not data_path.exists():
        print(f"Error: Data path does not exist: {data_path}")
        sys.exit(1)

    # Select adapter
    if args.dataset == "weihs":
        data_iterator = load_weihs_dataset(str(data_path), chunksize=args.chunksize)
    elif args.dataset == "hensen":
        data_iterator = load_hensen_dataset(str(data_path), chunksize=args.chunksize)
    elif args.dataset == "big_bell_test":
        data_iterator = load_big_bell_test_dataset(str(data_path), chunksize=args.chunksize)
    else:
        # Generic loader
        loader_func = get_loader(data_path)
        data_iterator = loader_func(data_path, chunksize=args.chunksize)

    # Run validation
    run_open_data_chsh_validation(
        data_iterator,
        dataset_name=args.dataset if args.dataset != "generic" else data_path.stem,
        output_dir=args.output,
        bootstrap_samples=args.bootstrap_samples,
        seed=args.seed
    )

if __name__ == "__main__":
    main()
