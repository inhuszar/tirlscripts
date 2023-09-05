#!/usr/bin/env python
# -*- coding: utf-8 -*-

#   _______ _____ _____  _                    _       _
#  |__   __|_   _|  __ \| |                  (_)     | |
#     | |    | | | |__) | |     ___  ___ _ __ _ _ __ | |_ ___
#     | |    | | |  _  /| |    / __|/ __| '__| | '_ \| __/ __|
#     | |   _| |_| | \ \| |____\__ \ (__| |  | | |_) | |_\__ \
#     |_|  |_____|_|  \_\______|___/\___|_|  |_| .__/ \__|___/
#                                              | |
#                                              |_|
#
# Copyright (C) 2018-2023 University of Oxford
# Part of the FMRIB Software Library (FSL)
# Author: Istvan N. Huszar


# SHBASECOPYRIGHT


"""
Applies a physical transformation to a NIfTI volume.

"""

__tirlscript__ = "stanford.applywarp"


# IMPORTS

import sys
import argparse
import numpy as np
import nibabel as nib

# TIRL IMPORTS

import tirl
from tirl.timage import TImage


# IMPLEMENTATION

def main(args):
    source = TImage(args.source, external="affine")
    target = TImage(args.target, external="affine")
    chain = tirl.load(args.warp)
    res = warp(source, target, chain)
    if args.output.endswith(".timg"):
        res.save(args.output, overwrite=False)
    elif args.output.endswith((".nii", ".nii.gz")):
        _snapshot(res, args.output)


def warp(source, target, chain=None):
    if chain is not None:
        target.domain.external += chain
    result = source.evaluate(target.domain)
    result.domain = result.domain[:1]
    return result


def _snapshot(img, fname):
    affine = np.eye(4)
    affine[:3] = img.domain.external.reduce()[0].matrix
    nifti = nib.Nifti1Image(img.data, affine)
    nib.save(nifti, fname)


def create_cli():
    parser = argparse.ArgumentParser(prog="applywarp")
    parser.add_argument(
        "-i", "--source", type=str,
        help="Source (input) image that will be warped.")
    parser.add_argument(
        "-r", "--target", type=str,
        help="Target image, used as a field-of-view reference.")
    parser.add_argument(
        "-w", "--warp", type=str,
        help="TIRL image transformation chain (*.img.chain).")
    parser.add_argument(
        "-o", "--output", type=str,
        help="Output file name (*.timg or *.nii or *.nii.gz).")
    return parser


if __name__ == "__main__":
    parser = create_cli()
    if len(sys.argv[1:]) > 0:
        args = parser.parse_args()
        main(args)
    else:
        parser.print_help()
