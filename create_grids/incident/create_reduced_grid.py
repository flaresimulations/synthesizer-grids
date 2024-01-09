"""
Code to make a reduced grid, for example limiting the number of age bins to a
specific set or a max age.
"""

import numpy as np
import argparse
from pathlib import Path
import yaml
import h5py

# synthesiser modules
from synthesizer.abundances import Abundances
from synthesizer.grid import Grid
from synthesizer.photoionisation import cloudy17 as cloudy


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reduce a grid")

    parser.add_argument(
        "-synthesizer_data_dir",
        type=str,
        required=True,
        help="path to synthesizer_data_dir",
    )

    parser.add_argument(
        "-original_grid",
        type=str,
        required=True,
        help="the name of the original_grid",
    )

    parser.add_argument(
        "-max_age",
        type=float,
        required=False,
        default=None,
        help="log10 of the max age",
    )

    parser.add_argument(
        "-ages",
        type=str,
        required=False,
        default=None,
        help="log10 of the specific set of ages, e.g. -ages=6.,7.,8.",
    )

    args = parser.parse_args()

    # open the parent incident grid
    original_grid = Grid(
        args.original_grid,
        grid_dir=f"{args.synthesizer_data_dir}",
        read_lines=False,
    )

    print(args.ages)

    ages = np.array(list(map(float, args.ages.split(","))))
    print(ages)

    # get name of new grid
    if args.ages:
        new_grid_name = f"{args.original_grid}-ages:{args.ages}"
    elif args.max_age:
        new_grid_name = f"{args.original_grid}-max_age:{args.max_age}"

    # create an indices array
    all_indices = np.indices(original_grid.log10age.shape)[0]

    # maximum index for the age

    if args.ages:
        indices = []
        for age in ages:
            print(age)
            indices.append(np.where(original_grid.log10age == age)[0][0])
        indices = np.array(indices)
    elif args.max_age:
        indices = all_indices[original_grid.log10age <= args.max_age]

    print(indices)

    # open the new grid
    with h5py.File(
        f"{args.synthesizer_data_dir}/{new_grid_name}.hdf5", "w"
    ) as hf:
        # open the original incident model grid
        with h5py.File(
            f"{args.synthesizer_data_dir}/{args.original_grid}.hdf5", "r"
        ) as hf_original:
            print("-" * 50)
            print(f"ORIGINAL GRID - {args.original_grid}")
            print("-" * 50)
            print("---- attributes")
            for k, v in hf_original.attrs.items():
                print(k, v)
            print("---- groups and datasets")
            hf_original.visititems(print)

            # copy top-level attributes
            for k, v in hf_original.attrs.items():
                hf.attrs[k] = v

            # copy axes, modifying the age axis
            hf["axes/metallicity"] = hf_original["axes/metallicity"][:]
            hf["axes/log10age"] = hf_original["axes/log10age"][indices]

            # copy log10Q
            for ion in hf_original["log10Q"].keys():
                hf[f"log10Q/{ion}"] = hf_original["log10Q"][ion][indices, :]

            # copy wavelength grid
            hf[f"spectra/wavelength"] = hf_original["spectra/wavelength"][:]

            # copy spectra
            hf[f"spectra/incident"] = hf_original["spectra/incident"][
                indices, :
            ]

        # print attributes and datasets of new grid
        print("-" * 50)
        print(f"NEW GRID - {new_grid_name}")
        print("-" * 50)
        print("---- attributes")
        for k, v in hf.attrs.items():
            print(k, v)
        print("---- groups and datasets")
        hf.visititems(print)
        print(hf["axes/log10age"][:])
