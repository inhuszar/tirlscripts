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


import os
import setuptools


TS = "tirlscripts"
OX = "oxford"
MAN = "manchester"


# Package definitions

OXFORD_PACKAGES = [
    "bigmac", "cliutils", "mnd",
    "mouse", "scriptutils", "segmentation", "slicealign", "tracer"
]
MANCHESTER_PACKAGES = [
    "nsurg"
]


def generate_package_name(package_group, package_name):
    return f"{TS}.{package_group}.{package_name}"


def generate_package_path(package_group, package_name):
    return os.path.join(
        f"{TS}-{package_group}-{package_name}",
        f"{TS}/{package_group}",
        f"{package_name}"
    )


# Find packages

def oxford(package_name):
    return (
        generate_package_name(OX, package_name),
        generate_package_path(OX, package_name)
    )


def manchester(package_name):
    return (
        generate_package_name(MAN, package_name),
        generate_package_path(MAN, package_name)
    )


# Install tirlscripts package

ALL_PACKAGES = {
    **dict([oxford(package_name) for package_name in OXFORD_PACKAGES]),
    **dict([manchester(package_name) for package_name in MANCHESTER_PACKAGES])
}

setuptools.setup(
    name=TS,
    packages=[*ALL_PACKAGES.keys()],
    package_dir=ALL_PACKAGES
)
