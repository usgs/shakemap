#!/usr/bin/env python

# stdlib imports
import argparse
import pathlib
import sys

# third party imports
from netCDF4 import Dataset


def is_axis_required(filename):
    root = Dataset(filename, "r", format="NETCDF4")
    axis_required = False
    if "lon" in root.variables and "lat" in root.variables:
        return False
    if "x" not in root.variables or "y" not in root.variables:
        return False
    xattrs = root.variables["x"].ncattrs()
    yattrs = root.variables["x"].ncattrs()
    if "axis" not in xattrs or "axis" not in yattrs:
        axis_required = True
    root.close()
    return axis_required


def fix_file(filename):
    if not is_axis_required(filename):
        return (False, "File {filename} does not need fixing.")
    root = Dataset(filename, "a", format="NETCDF4")
    # set the attribute on both variables
    xvar = root.variables["x"]
    yvar = root.variables["y"]
    xvar.axis = "X"
    yvar.axis = "Y"
    root.close()
    msg = f"File {filename} now has 'axis' attributes on 'x' and 'y' variables"
    return (True, msg)


def main(args):
    if args.check:
        msg = f"File {args.filename} does not need fixing."
        result = False
        if is_axis_required(args.filename):
            msg = f"File {args.filename} needs fixing."
            result = True
        print(msg)
        return not result

    infile = pathlib.Path(args.filename)
    if not infile.exists():
        print(f"Input file {infile} does not exist.")
        sys.exit(1)
    result, msg = fix_file(infile)
    print(msg)
    return result


if __name__ == "__main__":
    desc = "Fix NetCDF files to be compliant with GDAL 'axis' attribute on X and Y."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filename", help="Netcdf file to fix")
    parser.add_argument(
        "-c",
        "--check",
        help="Check to see if netcdf file needs to be fixed",
        action="store_true",
        default=False,
    )
    pargs = parser.parse_args()
    result = main(pargs)
    sys.exit(int(not result))
