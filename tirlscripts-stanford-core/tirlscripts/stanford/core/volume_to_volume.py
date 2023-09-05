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
Registers two 3D volumes.

"""

__tirlscript__ = "stanford.v2v"


# DEPENDENCIES

import os
import sys
import logging
import argparse
import numpy as np
from attrdict import AttrMap


# TIRL IMPORTS

import tirl.settings as ts
from tirl.chain import Chain
from tirl.timage import TField
from tirl.timage import TImage
from tirl.cost.mi import CostMI
from tirl.cost.msd import CostMSD
from tirl.cost.mind import CostMIND
from tirl.regularisation.diffusion import DiffusionRegulariser
from tirl.transformations.scale import TxScale, TxIsoScale
from tirl.transformations.affine import TxAffine
from tirl.transformations.euler import TxEulerAngles
from tirl.transformations.translation import TxTranslation
from tirl.transformations.displacement import TxDisplacementField
from tirl.optimisation.optgroup import OptimisationGroup
from tirl.optimisation.optnl import OptNL
from tirl.optimisation.gnoptdiff import GNOptimiserDiffusion


# TIRLSCRIPT IMPORTS

from tirl import srchash
from tirlscripts.oxford.mnd import __version__
from tirlscripts.oxford.scriptutils import general, image, inout
from tirlscripts.stanford.core import applywarp


# DEFINITIONS

from tirl.constants import *


# IMPLEMENTATION


def run(cnf=None, **options):
    """
    Runs TIRL volume-to-volume registration.

    :param cnf:
        Configuration file. Instead of a file, a dictionary with
        suitable content may also be specified.
    :type cnf: Union[str, dict, None]
    :param options:
        Overriding configuration parameters.
    :type options: Any

    """
    # Load script configuration
    if cnf is not None:
        if isinstance(cnf, dict):
            cnf = dict(cnf)
        elif isinstance(cnf, str):
            cnf = general.load_configurations(cnf)
        else:
            raise TypeError(
                f"Unrecognised configuration format: "
                f"{cnf.__.class__.__name__}")
    cnf.update(options)
    p, logger = general.initialise_script(**cnf)
    p.logger = logger.name  # avoid globals
    logger.debug(f"tirlscripts module version: {__version__}")
    logger.debug(f"volume_to_volume.py SHA1: {srchash(__file__)}")

    # Load and configure input images
    if p.target.export is True:
        ext = ts.EXTENSIONS["TImage"]
        p.target.export = os.path.join(p.general.outputdir, f"target.{ext}")
    target = inout.load_volume(**p.target)
    target.rule = None
    if p.source.export is True:
        ext = ts.EXTENSIONS["TImage"]
        p.source.export = os.path.join(p.general.outputdir, f"source.{ext}")
    source = inout.load_volume(**p.source)
    source.rule = None

    # Perform pre-registration actions on the images,
    # unless they were loaded from an alternative (TImage) source.
    # Actions are user-defined functions within the TIRL namespace. Actions
    # can be chain-loaded to perform preparatory analysis steps on the images
    # before registration begins.
    isalternative = p.target.file.lower().endswith(
        (ts.EXTENSIONS["TImage"], ts.EXTENSIONS["TIRLObject"]))
    if not isalternative:
        logger.info("Performing target image operations...")
        target = image.perform_image_operations(
            target, *p.preprocessing.target, scope=globals(), other=source)

    # Perform actions on the brain slice prior to registration, unless it was
    # loaded from a TImage file.
    isalternative = p.source.file.lower().endswith(
        (ts.EXTENSIONS["TImage"], ts.EXTENSIONS["TIRLObject"]))
    if not isalternative:
        logger.info("Performing source image operations...")
        source = image.perform_image_operations(
            source, *p.preprocessing.source, scope=globals(), other=target)

    # Run the registration routine
    try:
        register(target, source, p)
    except Exception as exc:
        logger.error(exc.args)
        logger.fatal(f"The registration terminated with an exception.")
        raise exc
    else:
        logger.fatal("The registration was completed successfully.")


def initialise_transformations(fixed, q):
    """
    Creates transformation chain that will be optimised.

    :returns: linear + nonlinear transformation chain with initial parameters
    :rtype: Chain

    """
    # Scale
    lb = np.asarray(q.init.scale.lb)
    ub = np.asarray(q.init.scale.ub)
    tx_scale = TxIsoScale(
        *q.init.scale.x0, dim=3, bounds=(lb, ub), name="scale")
    lb = np.asarray(q.init["anisotropic-scale"]["lb"])
    ub = np.asarray(q.init["anisotropic-scale"]["ub"])
    tx_anisotropic_scale = TxScale(
        1., 1., 1., bounds=(lb, ub), name="anisotropic_scale")

    # Rotation
    tx_initial_rotation = TxEulerAngles(
        *q.init.rotation.x0,
        order=str(q.init.rotation.order),
        mode=q.init.rotation.mode,
        bounds=(-np.inf, np.inf), name="initial_rotation")
    if str(q.init.rotation.mode).lower() == "deg":
        lb = np.radians(q.init.rotation.lb)
        ub = np.radians(q.init.rotation.ub)
    else:
        lb = np.asarray(q.init.rotation.lb)
        ub = np.asarray(q.init.rotation.ub)
    tx_rotation = TxEulerAngles(
        0, 0, 0,
        order=str(q.init.rotation.order),
        mode=q.init.rotation.mode,
        bounds=(lb, ub), name="rotation")

    # Translation
    tx_initial_translation = TxTranslation(
        *q.init.translation.x0,
        bounds=(-np.inf, np.inf),
        name="initial_translation"
    )
    lb = np.asarray(q.init.translation.lb)
    ub = np.asarray(q.init.translation.ub)
    tx_translation = TxTranslation(
        0, 0, 0,
        bounds=(lb, ub),
        name="translation"
    )

    # Affine
    x0 = np.asarray(q.init.affine.x0).reshape((3, 4))
    lb = np.asarray(q.init.affine.lb)
    ub = np.asarray(q.init.affine.ub)
    tx_affine = TxAffine(x0, bounds=(lb, ub), name="affine")

    # Append linear transformations to the domain equivalent to the fixed image
    linear_chain = Chain(
        tx_initial_rotation,
        tx_rotation,
        tx_scale,
        tx_anisotropic_scale,
        tx_initial_translation,
        tx_translation,
        tx_affine
    )
    domain = fixed.domain[:]
    domain.chain.extend(linear_chain)

    # Nonlinear
    x0 = float(q.init.nonlinear.x0) * np.ones((3, *fixed.vshape))
    if q.init.nonlinear.lb is None:
        lb = None
    else:
        lb = float(q.init.nonlinear.lb) * np.ones_like(x0)
    if q.init.nonlinear.ub is None:
        ub = None
    else:
        ub = float(q.init.nonlinear.ub) * np.ones_like(x0)
    field = TField(
        x0, copy=False, domain=domain[:], taxes=(0,), order=TENSOR_MAJOR)
    tx_nonlinear = TxDisplacementField(
        field, bounds=(lb, ub), name="nonlinear", mode=NL_REL)

    # Return the full transformation chain
    return Chain(*linear_chain, tx_nonlinear)


def register(fixed, moving, cnf):
    """
    Runs the four registration stages: rigid-coarse, rigid-fine, affine, and
    non-linear. The function has no return value. The transformation chain is
    attached to the Domain of the fixed TImage object and is optimised in situ.

    :param fixed:
        fixed image (to which the chain is attached)
    :type fixed: TImage
    :param moving:
        moving image (that defines the coordinate space that the mapping
        is into)
    :type moving: TImage
    :param cnf: all configuration options
    :type cnf: AttrMap or dict

    """
    p = AttrMap(cnf)
    q = p.regparams
    logger = logging.getLogger(p.logger)

    # Create transformation chain (does not change the fixed image)
    # rotation -> scale -> translation -> affine -> nonlinear
    logger.info("Initialising transformation chain...")
    chain = initialise_transformations(fixed, q)
    logger.info("Transformation chain has been initialised.")

    # Generate output: (best) initial alignment
    fixed.save(os.path.join(
        p.general.outputdir, "fixed.timg"), overwrite=True)
    applywarp._snapshot(
        fixed,
        os.path.join(p.general.outputdir, "fixed.nii.gz"))
    if fixed.mask is not None:
        applywarp._snapshot(
            fixed.tmask(),
            os.path.join(p.general.outputdir, "fixed_mask.nii.gz"))
    moving.save(os.path.join(
        p.general.outputdir, "moving.timg"), overwrite=True)
    applywarp._snapshot(
        moving,
        os.path.join(p.general.outputdir, "moving.nii.gz"))
    if moving.mask is not None:
        applywarp._snapshot(
            moving.tmask(),
            os.path.join(p.general.outputdir, "moving_mask.nii.gz"))

    # Coarse rigid registration
    fixed.domain.chain.extend(chain[:-2])
    if "rigid-coarse" in p.general.stages:
        logger.info("Starting coarse rigid registration...")
        best_initialisations = rigid3d_coarse(fixed, moving, p)
        logger.info("Completed coarse rigid registration.")
    else:
        logger.info("Coarse rigid registration was skipped.")
        best_initialisations = ((np.zeros(3), np.zeros(3)),)

    # Report on best initialisation(s)
    logger.info(f"Best initialisation(s):")
    logger.info(best_initialisations)
    initsfile = os.path.join(p.general.outputdir, "best_inits.txt")
    with open(initsfile, "w") as fp:
        fp.write("Rotation[XYZ:degrees],Translation[xyz:mm]\n")
        tx_rot = fixed.domain.external["initial_rotation"]
        tx_trans = fixed.domain.external["initial_translation"]
        for ix, init in enumerate(best_initialisations):
            # Save NIfTI image
            rotations, translations = init
            tx_rot.parameters[:] = np.radians(rotations)
            tx_trans.parameters[:] = translations
            applywarp._snapshot(
                applywarp.warp(moving, fixed),
                os.path.join(p.general.outputdir,
                             f"moving1_init{ix}.nii.gz")
            )
            # Save initialisation parameters
            rotations = str(rotations.tolist())
            translations = str(translations.tolist())
            fp.write(f"{rotations},{translations}\n")

    # Fine rigid registration (7 DOF)
    if "rigid-fine" in p.general.stages:
        logger.info("Starting fine rigid registration...")
        rigid3d_fine(fixed, moving, best_initialisations, p)
        logger.info("Completed fine rigid registration.")
        # Generate output
        fixed.save(os.path.join(
            p.general.outputdir, "fixed2_rigid.timg"), overwrite=True)
        applywarp._snapshot(
            applywarp.warp(moving, fixed),
            os.path.join(p.general.outputdir, "moving2_rigid.nii.gz"))
    else:
        logger.info("Fine rigid registration was skipped.")

    # Anisotropic scaling (9 DOF)
    if "anisotropic-scaling" in p.general.stages:
        logger.info("Starting anisotropic scaling optimisation (9 DOF)...")
        anisotropic_scaling(fixed, moving, p)
        logger.info("Completed anisotropic scaling optimisation.")
        # Generate output
        fixed.save(os.path.join(
            p.general.outputdir, "fixed3_scaling.timg"), overwrite=True)
        applywarp._snapshot(
            applywarp.warp(moving, fixed),
            os.path.join(p.general.outputdir, "moving3_scaling.nii.gz"))
    else:
        logger.info("Anisotropic scaling optimisation was skipped.")

    # Affine registration
    fixed.domain.chain.append(chain[-2])
    if "affine" in p.general.stages:
        logger.info("Starting affine registration...")
        affine3d(fixed, moving, p)
        logger.info("Completed affine registration.")
        # Generate output
        fixed.save(os.path.join(
            p.general.outputdir, "fixed4_affine.timg"), overwrite=True)
        applywarp._snapshot(
            applywarp.warp(moving, fixed),
            os.path.join(p.general.outputdir, "moving4_affine.nii.gz"))
    else:
        logger.info("Affine registration was skipped.")

    # Non-linear registration
    tx_nonlinear = chain[-1]
    tx_nonlinear.domain.chain = fixed.domain.chain[:]
    fixed.domain.chain.append(tx_nonlinear)
    if "nonlinear" in p.general.stages:
        logger.info("Starting non-linear registration...")
        diffreg3d(fixed, moving, p)
        logger.info("Completed non-linear registration.")
        # Generate output
        fixed.save(os.path.join(
            p.general.outputdir, "fixed5_nonlinear.timg"), overwrite=True)
        applywarp._snapshot(
            applywarp.warp(moving, fixed),
            os.path.join(p.general.outputdir, "moving5_nonlinear.nii.gz"))
    else:
        logger.info("Non-linear registration was skipped.")

    # Save transformations in an intuitive format
    _save_transformations(fixed, moving, p.general.outputdir, logger)


def _save_transformations(fixed, moving, outputdir, logger):
    """
    Export complete transformation chains that map between the physical spaces
    of the NIfTI volumes. Note that the filenames represent the direction of
    the image registration rather than the direction of mapping coordinates.

    """

    # Forward transformation (maps target coordinates to source)
    ftx_file = os.path.join(outputdir, "source_to_target.img.chain")
    f2m = _physical(fixed.domain.external) \
          + _physical(moving.domain.external).inverse()
    try:
        f2m.save(ftx_file)
    except Exception:
        logger.error(f"Could not save forward transformation: {ftx_file}")
    else:
        logger.info(f"Saved forward transformation: {ftx_file}")

    # Backward transformation (maps source coordinates to target coordinates)
    logger.info(
        "Inverting the transformation (this may take several minutes)...")
    m2f = f2m.inverse()
    btx_file = os.path.join(outputdir, "target_to_source.img.chain")
    try:
        m2f.save(btx_file)
    except Exception:
        logger.error(f"Could not save backward transformation: {btx_file}")
    else:
        logger.info(f"Saved backward transformation: {btx_file}")


def _physical(chain):
    """
    If the first transformation in a chain is a NIfTI affine (sform or qform),
    this function returns the subchain that follows the affine. The rationale
    is that the resulting chain provides a mapping from the physical space of
    the NIfTI volume, rather than the voxel space.

    """
    if len(chain) == 0:
        return chain
    affine = chain[0]
    niftiaffine = ("affine", "sform", "qform")
    if isinstance(affine, TxAffine) and (affine.name in niftiaffine):
        return chain[1:]
    else:
        return chain


def rigid3d_coarse(fixed, moving, cnf):
    """
    Performs a coarse testing of 3D rotations and translations (6 DOF)
    to find the most suitable initialisation.

    Args:
        fixed:
            fixed TImage
        moving:
            moving TImage
        cnf:
            program configurations

    Returns:
        N-Tuple of initialisation parameters in the form of
        ((rot_x, rot_y, rot_z). (trans_x, trans_y, trans_z)),
        where rotations are in degrees and translations are in arbitrary
        physical units (usually mm).

    """
    p = AttrMap(cnf)
    q = AttrMap(p.regparams["rigid-coarse"])
    logger = logging.getLogger(p.logger)

    # Build sequence of initial rotations
    if q.sequence.combinations:
        initial_rotations = np.stack(
            np.meshgrid(
                q.sequence.x, q.sequence.y, q.sequence.z,
                indexing="ij"),
            axis=-1).reshape(-1, 3)
    else:
        initial_rotations = np.stack(
            (q.sequence.x, q.sequence.y, q.sequence.z), axis=-1)

    # Build sequence of initial translations
    if q.xrange and q.xstops:
        initial_translations = np.stack(
            tuple(
                np.linspace(-half, half, stops, endpoint=True)
                for (half, stops) in zip(np.divide(q.xrange, 2), q.xstops)
            ), axis=-1)
    else:
        initial_translations = ((0, 0, 0),)

    # Sample the images to the target resolution
    _fixed = fixed.rescale(float(q.scale), copy=True)
    _moving = moving.rescale(float(q.scale), copy=True)

    # Prepare rigid transformations
    tx_rot = fixed.domain.external["initial_rotation"]
    tx_trans = fixed.domain.external["initial_translation"]

    # Create cost object
    cost = _cost_object(_fixed, _moving, cost_type=p.general.cost)

    # Measure grid costs
    gridcosts = []
    for rotation in initial_rotations:
        tx_rot.parameters[:] = np.radians(rotation)  # bounds are -inf..inf
        for translation in initial_translations:
            tx_trans.parameters[:] = translation  # bounds are -inf..inf
            cost_value = cost()
            gridcosts.append(cost_value)
            logger.debug(
                f"Rotations: {rotation} [deg], "
                f"Translations: {translation} [mm]: "
                f"Cost: {cost_value}"
            )

    # Return N best initialisations
    n_best = int(q.n_best)
    indices = np.argsort(gridcosts)[:n_best]
    best_initialisations = []
    n_trans = len(initial_translations)
    for ix in indices:
        best_initialisations.append(
            (initial_rotations[ix // n_trans],
             initial_translations[ix % n_trans])
        )

    return best_initialisations


def _cost_object(fixed, moving, cost_type=None):
    """
    Creates a specific type of cost object for the specified pair
    of fixed and moving TImages.

    """
    cost_type = str(cost_type).upper()
    if cost_type == "MSD":
        cost = CostMSD(moving, fixed, normalise=True)
    elif cost_type == "MI":
        cost = CostMI(moving, fixed, normalise=True)
    elif cost_type == "MIND":
        cost = CostMIND(moving, fixed, normalise=True,
                        kernel=MK_FULL, sigma=1, truncate=1.5)
    else:
        raise ValueError("Unsupported cost function.")

    return cost


def rigid3d_fine(fixed, moving, initialisations, cnf):
    """
    Performs fine-tuned optimisation over the rigid-body parameters and
    isotropic scaling (7 DOF).

    """

    p = AttrMap(cnf)
    logger = logging.getLogger(p.logger)

    costs_and_chains = []
    for ix, init in enumerate(initialisations):
        logger.info(f"Fine tuning 7 DOF at Initialisation #{ix}...")
        rotations, translations = init
        tx_rot = fixed.domain.external["initial_rotation"]
        tx_trans = fixed.domain.external["initial_translation"]
        tx_rot.parameters[:] = np.radians(rotations)
        tx_trans.parameters[:] = np.asarray(translations)
        opt_cost = _rigid3d_fine_inner_loop(fixed, moving, cnf)
        costs_and_chains.append((opt_cost, fixed.domain.external.copy()))
    else:
        # Reassign best transformation chain
        chain = min(costs_and_chains, key=lambda t: t[0])[1]
        fixed.domain.external = chain


def anisotropic_scaling(fixed, moving, cnf):

    p = AttrMap(cnf)
    q = AttrMap(p.regparams["anisotropic-scaling"])
    logger = logging.getLogger(p.logger)

    # Scaling-smoothing iteration
    for i, (sc, sm) in enumerate(zip(q.scaling, q.smoothing)):
        logger.debug(f"Scale: {sc}, smoothing: {sm} px...")
        # Prepare images for the current iteration
        fixed.rescale(1. / sc, copy=False)
        moving.rescale(1. / sc, copy=False)
        fixed_smooth = fixed.smooth(sm, copy=True)
        moving_smooth = moving.smooth(sm, copy=True)
        # Prepare co-optimised transformations
        tx_rot = fixed_smooth.domain.external["rotation"]
        tx_scale = fixed_smooth.domain.external["anisotropic_scale"]
        tx_trans = fixed_smooth.domain.external["translation"]
        og = OptimisationGroup(tx_rot, tx_scale, tx_trans)

        # Set cost function
        cost = _cost_object(fixed_smooth, moving_smooth, p.general.cost)

        # Start optimisation
        OptNL(og, cost, method="LN_BOBYQA", visualise=q.visualise,
              xtol_abs=q.xtol_abs, xtol_rel=q.xtol_rel, step=q.opt_step,
              logger=logger, normalised=True)()
        # Transfer optimised transformations to the non-smoothed images
        fixed.domain = fixed_smooth.domain
        fixed.resmgr.sync()
        moving.domain = moving_smooth.domain
        moving.resmgr.sync()
    else:
        # Restore full resolution of the images
        fixed.rescale(1, copy=False)
        moving.rescale(1, copy=False)


def _update_bounds(tx, lb=None, ub=None, autotolerance=0.1):
    # Create parameter range automatically
    params = tx.parameters[:]
    if lb is None:
        lb = params - autotolerance * np.abs(params)
    if ub is None:
        ub = params + autotolerance * np.abs(params)
    # Set bounds
    tx.parameters.set_bounds(lb, ub)


def _rigid3d_fine_inner_loop(fixed, moving, cnf):

    p = AttrMap(cnf)
    q = AttrMap(p.regparams["rigid-fine"])
    logger = logging.getLogger(p.logger)
    final_cost = np.inf

    # Scaling-smoothing iteration
    for i, (sc, sm) in enumerate(zip(q.scaling, q.smoothing)):
        logger.debug(f"Scale: {sc}, smoothing: {sm} px...")
        # Prepare images for the current iteration
        fixed.rescale(1. / sc, copy=False)
        moving.rescale(1. / sc, copy=False)
        fixed_smooth = fixed.smooth(sm, copy=True)
        moving_smooth = moving.smooth(sm, copy=True)
        # Prepare co-optimised transformations
        tx_rot = fixed_smooth.domain.external["rotation"]
        tx_scale = fixed_smooth.domain.external["scale"]
        tx_trans = fixed_smooth.domain.external["translation"]
        og = OptimisationGroup(tx_rot, tx_scale, tx_trans)

        # Set cost function
        cost = _cost_object(fixed_smooth, moving_smooth, p.general.cost)

        # Start optimisation
        OptNL(og, cost, method="LN_BOBYQA", visualise=q.visualise,
              xtol_abs=q.xtol_abs, xtol_rel=q.xtol_rel, step=q.opt_step,
              logger=logger, normalised=True)()
        # Transfer optimised transformations to the non-smoothed images
        fixed.domain = fixed_smooth.domain
        fixed.resmgr.sync()
        moving.domain = moving_smooth.domain
        moving.resmgr.sync()
        final_cost = cost()
    else:
        # Restore full resolution of the images
        fixed.rescale(1, copy=False)
        moving.rescale(1, copy=False)

    return final_cost


def affine3d(fixed, moving, cnf):
    """
    Optimises a 12-DOF affine transformation assuming that the images
    have been rigidly aligned.

    """
    p = AttrMap(cnf)
    q = p.regparams.affine
    logger = logging.getLogger(p.logger)

    # Scaling-smoothing iterations
    for i, (sc, sm) in enumerate(zip(q.scaling, q.smoothing)):
        logger.debug(f"Scale: {sc}, smoothing: {sm} px...")
        # Prepare images for the current iteration
        fixed.rescale(1. / sc, copy=False)
        moving.rescale(1. / sc, copy=False)
        fixed_smooth = fixed.smooth(sm, copy=True)
        moving_smooth = moving.smooth(sm, copy=True)
        # Prepare transformation to optimise
        tx_affine = fixed_smooth.domain.external["affine"]

        # Set cost function
        cost = _cost_object(fixed_smooth, moving_smooth, p.general.cost)

        # Start optimisation
        OptNL(tx_affine, cost, method="LN_BOBYQA",
              xtol_rel=q.xtol_rel, xtol_abs=q.xtol_abs,
              visualise=q.visualise, step=q.opt_step, logger=logger,
              normalised=True)()
        # Transfer optimised transformations to the non-smoothed images
        fixed.domain = fixed_smooth.domain
        fixed.resmgr.sync()
        moving.domain = moving_smooth.domain
        moving.resmgr.sync()
    else:
        # Restore full resolution of the images
        fixed.rescale(1, copy=False)
        moving.rescale(1, copy=False)


def diffreg3d(fixed, moving, cnf):
    """
    Performs a 3D non-linear registration. The transformation is parameterised
    as a dense displacement field. The cost is MIND, and diffusion
    regularisation is used to create smoothness in the deformation field.

    """
    p = AttrMap(cnf)
    q = p.regparams.nonlinear
    logger = logging.getLogger(p.logger)

    # Scaling-smoothing iteration
    for i, (sc, sm) in enumerate(zip(q.scaling, q.smoothing)):
        logger.debug(f"Scale: {sc}, smoothing: {sm} px...")
        # Prepare images for the current iteration
        fixed.rescale(1 / sc, copy=False)
        moving.rescale(1 / sc, copy=False)
        fixed_smooth = fixed.smooth(sm, copy=True)
        fixed_smooth.storage = MEM
        moving_smooth = moving.smooth(sm, copy=True)
        moving_smooth.storage = MEM
        # Prepare transformation to optimise
        tx_nonlinear = fixed_smooth.domain.external[-1]
        # Set cost and regulariser
        cost = CostMIND(moving_smooth, fixed_smooth, sigma=float(q.sigma),
                        truncate=float(q.truncate), kernel=MK_FULL)
        regularisation = DiffusionRegulariser(
            tx_nonlinear, weight=float(q.regweight))
        # Optimise the non-linear transformation
        GNOptimiserDiffusion(
            tx_nonlinear, cost, regularisation, maxiter=int(q.maxiter[i]),
            xtol_rel=q.xtol_rel, xtol_abs=q.xtol_abs, visualise=q.visualise,
            logger=logger)()
        # Transfer optimised transformations to the non-smoothed images
        fixed.domain = fixed_smooth.domain
        fixed.resmgr.sync()
        moving.domain = moving_smooth.domain
        moving.resmgr.sync()
    else:
        # Restore the original resolution of the images
        fixed.rescale(1, copy=False)
        moving.rescale(1, copy=False)


def centralise(img, history=None, other=None):
    """
    Aligns a TImage with its geometrical centre to the origin.

    """
    img.centralise(weighted=False)
    return img


def weighted_centralise(img, history=None, other=None):
    """
    Aligns a TImage with its centre of gravity to the origin. To calculate the
    centre of gravity, coordinates are weighted by the intensity values of
    the image. Masks are discarded from the calculation.

    """
    img.centralise(weighted=True)
    return img


def match_resolution(img, history=None, other=None):
    """
    Changes the resolution of the input image to the resolution of the
    other image, assuming that both images are defined on a domain with an
    affine transformation.

    """
    # Resolution = diagonal of the upper triangular matrix
    # in the QR decomposition of the affine matrix
    img_res = np.diag(np.linalg.qr(_get_affine(img).matrix[:3, :3])[1])
    other_res = np.diag(np.linalg.qr(_get_affine(other).matrix[:3, :3])[1])
    factors = np.divide(other_res, img_res)
    img.rescale(*factors, copy=False)


def dilate_mask(img, history=None, other=None):
    """
    Dilates the binary mask of the input image by a fixed amount.

    """
    from scipy.ndimage import binary_dilation

    mask = img.mask
    if mask is None:
        return img

    for i in range(16):
        mask = binary_dilation(mask).astype(np.uint8)
    else:
        img.mask = mask.astype(np.float32)

    return img


def _get_affine(img):
    """
    Returns the first affine transformation matrix in the external part of
    the TImage Domain's transformation chain.

    """
    for tx in img.domain.external:
        if isinstance(tx, TxAffine):
            return tx


def volume_to_volume(args):
    """
    Main program code. Controls program flow, handles command-line arguments.

    """
    cnf = args.config
    if os.path.isfile(cnf):
        cnf = general.load_configurations(cnf)
    else:
        raise FileNotFoundError(f"The provided configuration file "
                                f"does not exist: {args.config}")

    # Override tissue block and brain slice file paths in the configuration
    # file with those provided on the command line.
    if args.block and os.path.isfile(args.block):
        cnf["block"]["file"] = args.block
    if args.slice and os.path.isfile(args.slice):
        cnf["slice"]["file"] = args.slice
    if args.init:
        if len(args.init) > 1:
            y, x = args.init
            cnf["regparams"]["sites"] = np.atleast_2d([float(y), float(x)]).tolist()
        elif os.path.isfile(args.init[0]):
            sites = np.loadtxt(args.init[0])
            cnf["regparams"]["sites"] = np.atleast_2d(sites).tolist()
        else:
            raise ValueError("Invalid insertion site specification.")
    if args.out:
        cnf["general"]["outputdir"] = args.out

    # Override verbose option in configurations
    cnf["general"]["verbose"] = args.verbose

    # Run registration script
    run(cnf)


def create_cli(parser):
    """
    Sets up the CLI argument parser instance.

    """
    parser.add_argument("--block", metavar="image", default=None, type=str,
                        help="Tissue block photo", required=False)
    parser.add_argument("--slice", metavar="image", default=None, type=str,
                        help="Brain slice photo", required=False)
    parser.add_argument("--init", metavar=("y", "x"), nargs="+",
                        help="Insertion site", default=None, required=False)
    parser.add_argument("--out", metavar="dir", default=None, type=str,
                        help="Output directory", required=False)
    parser.add_argument("--config", metavar="cnf_file.yml", default=None,
                        type=str, help="Configuration file", required=True)
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="Print status messages to the command line",
                        required=False)

    return parser


def main(*args):
    """ Main program code. """

    parser = argparse.ArgumentParser(
        prog="volume_to_volume",
        description="Registers two 3D volumes.")
    parser = create_cli(parser)

    if args:
        volume_to_volume(parser.parse_args(args))
    else:
        parser.print_help()


# Program execution starts here
if __name__ == "__main__":
    main(*sys.argv[1:])
