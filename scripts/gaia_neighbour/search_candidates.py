"""Find the historical Gaia cnt1>1 quasar sample from separately supplied local data."""

import argparse
from pathlib import Path

from dual_agn.one_arcsecond_search import input_paths, resolve_config, run_search
from dual_agn.settings import ONE_ARCSECOND_SEARCH_ROOT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, default=ONE_ARCSECOND_SEARCH_ROOT)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--check-inputs', action='store_true')
    args = parser.parse_args()
    if not args.config.is_file():
        parser.error(f'Missing {args.config}; run scripts.shared.prepare_quasar_catalogue to generate the local discovery configuration, or use the main workflow.')
    config = resolve_config(args.config)
    if args.check_inputs:
        print(f'Found {len(input_paths(config))} local discovery inputs; mode={config["input_mode"]}. No network queries are made.')
        return
    run_search(config, args.output_root, args.resume)


if __name__ == '__main__':
    main()
