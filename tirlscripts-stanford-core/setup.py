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


from setuptools import setup
from tirlscripts.oxford.mnd import __version__


#---------------------------- TIRLscripts Installer ---------------------------#

setup(
    name="tirlscripts-oxford-mnd",
    version=__version__,
    description="Four-stage histology-to-MRI registration pipeline for TIRL",
    author="Istvan N. Huszar",
    author_email="istvan.huszar@ndcn.ox.ac.uk",
    license="Apache v2 Software License",
    packages=["tirlscripts.oxford.mnd"]
)
