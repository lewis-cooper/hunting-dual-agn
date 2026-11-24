"""Find imaging pairs from local master, Tractor, Gaia and CCD catalogues."""

import argparse
from pathlib import Path

from imaging_search.catalogues import input_paths, resolve_config
from imaging_search.pipeline import run_search
from settings import WORK_ROOT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True,help='Local catalogue paths and coverage metadata')
    parser.add_argument('--stage',choices=['all','discover'],default='all',
                        help='discover finds pairs and writes Gaia/CCD query positions without requiring those data')
    parser.add_argument('--output-root',type=Path,default=WORK_ROOT/'imaging_search')
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--check-inputs',action='store_true')
    args = parser.parse_args()
    if not args.config.is_file():
        parser.error(f'Missing {args.config}; copy imaging_search.example.json and set local data paths. No download will be attempted.')
    config = resolve_config(args.config)
    if args.check_inputs:
        paths = input_paths(config,args.stage=='discover')
        print(f'Found {len(paths)} local input files. This stage has no network/data-download functionality.')
        return
    run_search(config,args.output_root,args.stage,args.resume)


if __name__ == '__main__':
    main()
