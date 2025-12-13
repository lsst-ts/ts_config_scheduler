# This file is part of ts_config_scheduler.
#
# Developed for Vera C. Rubin Observatory Telescope and Site Systems.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__all__ = ("get_scheduler",)


import lsst.ts.fbs.utils.maintel.lsst_surveys as lsst_surveys
import numpy as np
import rubin_scheduler.scheduler.detailers as detailers
from lsst.ts.fbs.utils.maintel.lsst_surveys import safety_masks
from rubin_scheduler.scheduler import basis_functions
from rubin_scheduler.scheduler.schedulers import CoreScheduler
from rubin_scheduler.scheduler.surveys import GreedySurvey
from rubin_scheduler.scheduler.utils import (
    CurrentAreaMap,
    make_rolling_footprints,
)
from rubin_scheduler.site_models import Almanac
from rubin_scheduler.utils import SURVEY_START_MJD


def get_scheduler() -> tuple[int, CoreScheduler]:
    """Construct a limited LSST survey scheduler for chaos monkey mode.

    The parameters are not accessible when calling as 'config'.

    Returns
    -------
    nside : `int`
        Healpix map resolution.
    scheduler : `rubin_scheduler.scheduler.scheduler.CoreScheduler`
        Feature based scheduler.
    """
    nside = 32
    science_program = "BLOCK-T648"
    # Configured to only use one bandpass.
    desired_band = "z"

    band_to_filter = {
        "u": "u_24",
        "g": "g_6",
        "r": "r_57",
        "i": "i_39",
        "z": "z_20",
        "y": "y_10",
    }
    exptime = 30
    nexp = 1
    u_exptime = 38
    u_nexp = 1

    survey_start_mjd = SURVEY_START_MJD

    # Safety mask parameters - constraints on all survey pointings
    # Generally shadow_minutes value is set by the survey, but can
    # be set here as well (will be overwritten if too short for survey).
    safety_mask_params = {
        "nside": nside,
        "wind_speed_maximum": 40,
        "apply_time_limited_shadow": True,
        "time_to_sunrise": 3.0,
        "min_az_sunrise": 150,
        "max_az_sunrise": 250,
        "min_alt": 40,
    }

    # General parameters for standard pairs (-80/80 default)
    camera_rot_limits = (-5.0, 5.0)

    # Parameters for rolling cadence footprint definition
    nslice = 2  # N slices for rolling
    rolling_scale = 0.9  # Strength of rolling
    rolling_uniform = True  # Should we use the uniform rolling flag

    # Generate footprint over the sky
    sky = CurrentAreaMap(nside=nside)
    footprints_hp_array, labels = sky.return_maps()
    # Identify pixels for rolling
    roll_indx = np.where((labels == "lowdust") | (labels == "virgo"))[0]
    roll_footprint = footprints_hp_array["r"] * 0
    roll_footprint[roll_indx] = 1

    footprints_hp = {}
    for key in footprints_hp_array.dtype.names:
        footprints_hp[key] = footprints_hp_array[key]

    # Set up a mask to contain some surveys within this region
    footprint_mask = footprints_hp["r"] * 0
    footprint_mask[np.where(footprints_hp["r"] > 0)] = 1

    # And now remove all except desired band
    # This restricted to one band for AOS
    for band in footprints_hp:
        if band != desired_band:
            footprints_hp[band] *= 0

    # Use the Almanac to find the position of the sun at the start of survey
    almanac = Almanac(mjd_start=survey_start_mjd)
    sun_moon_info = almanac.get_sun_moon_positions(survey_start_mjd)
    sun_ra_start = sun_moon_info["sun_RA"].copy()

    # Define the rolling footprint
    footprints = make_rolling_footprints(
        fp_hp=footprints_hp,
        mjd_start=survey_start_mjd,
        sun_ra_start=sun_ra_start,
        nslice=nslice,
        scale=rolling_scale,
        nside=nside,
        wfd_indx=roll_indx,
        order_roll=1,
        n_cycles=3,
        uniform=rolling_uniform,
    )

    # Define the greedy surveys (single-visit per call)
    greedy = lsst_surveys.gen_greedy_surveys(
        nside=nside,
        camera_rot_limits=camera_rot_limits,
        exptime=exptime,
        nexp=nexp,
        u_exptime=u_exptime,
        u_nexp=u_nexp,
        footprints=footprints,
        science_program=science_program,
        safety_mask_params=safety_mask_params,
    )

    # Add chaos monkey perturbation blocks.
    perturbation_time_gap = 15.0  # Gap between perturbation applications.
    perturbation_base_block = "BLOCK-T648"
    perturbation_survey_name_base = "chaos_block"

    perturbation_surveys = []

    for i in range(6):
        perturbation_block = f"{perturbation_base_block}_{i+1}"
        perturbation_survey_name = f"{perturbation_survey_name_base}_{i+1}"
        perturbation_basis_functions = safety_masks(
            **safety_mask_params, shadow_minutes=0
        ) + [
            basis_functions.VisitGap(
                note=perturbation_survey_name,
                gap_min=(i + 1) * perturbation_time_gap,
            ),
            basis_functions.SlewtimeBasisFunction(bandname=None, nside=nside),
        ]

        perturbation_basis_weights = np.ones(len(perturbation_basis_functions))
        # Make the FBS not run away from last pointing
        perturbation_basis_weights[-1] = 100
        perturbation_surveys.append(
            GreedySurvey(
                perturbation_basis_functions,
                perturbation_basis_weights,
                nside=nside,
                survey_name=perturbation_survey_name,
                observation_reason=perturbation_survey_name,
                scheduler_note=perturbation_survey_name,
                science_program=perturbation_block,
                # We will try to tell the FBS how  long it takes to actually
                # acquire these visits here, for more accurate sims.
                nexp=1,
                exptime=60,
                bandname=desired_band,
                detailers=[
                    detailers.LabelRegionsAndDDFs(),
                ],
            )
        )

    # Arrange the surveys in tiers.
    surveys = [
        perturbation_surveys,
        greedy,
    ]

    # Combine into CoreScheduler
    scheduler = CoreScheduler(
        surveys,
        nside=nside,
        survey_start_mjd=survey_start_mjd,
        band_to_filter=band_to_filter,
    )

    return nside, scheduler


if __name__ == "config":
    nside, scheduler = get_scheduler()
