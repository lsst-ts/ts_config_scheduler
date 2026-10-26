#!/usr/bin/env python3
# This file is part of ts_config_scheduler.
#
# Developed for the Vera C. Rubin Observatory Telescope and Site Systems.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

__all__ = ("get_scheduler",)

import copy

import lsst.ts.fbs.utils.maintel.lsst_footprints as lsst_footprints
import lsst.ts.fbs.utils.maintel.lsst_surveys as lsst_surveys
import lsst.ts.fbs.utils.maintel.roman_surveys as roman_surveys
import lsst.ts.fbs.utils.maintel.too_surveys as too_surveys
import numpy as np
import rubin_scheduler.scheduler.basis_functions as bf
import rubin_scheduler.scheduler.detailers as detailers
from lsst.ts.fbs.utils.maintel.lsst_ddf_presched import read_ddf_obs_array
from rubin_scheduler.scheduler.schedulers import BaseQueueManager, CoreScheduler
from rubin_scheduler.scheduler.surveys import ScriptedSurvey
from rubin_scheduler.utils import DEFAULT_NSIDE

CAMERA_ROT_LIMITS = (-80.0, 80.0)
SURVEY_START_MJD = lsst_footprints.SURVEY_START_MJD


def generate_qm(
    rot_tel_pos_limits: tuple[float, float] = CAMERA_ROT_LIMITS,
    nside: int = DEFAULT_NSIDE,
    cloud_limit: float = 1.5,
) -> BaseQueueManager:
    """Generate a QueueManager object."""

    detailer_list = []
    # This detailer updates rotSkyPos if rotTelPos became
    # out of bounds or if rotSkyPos was not (yet) calculated.
    detailer_list.append(detailers.RotspUpdateDetailer(rot_limits=rot_tel_pos_limits))
    bf_list = []
    # This should get zenith masked without having to recalculate alt/az.
    bf_list.append(bf.SlewtimeBasisFunction())
    # Do not observe into any clouds currently at the pointing.
    bf_list.append(
        bf.MaskCloudMapBasisFunction(nside=nside, extinction_limit=cloud_limit)
    )

    qm = BaseQueueManager(
        detailers=detailer_list, basis_functions=bf_list, check_clouds=True
    )
    return qm


def get_scheduler() -> tuple[int, CoreScheduler]:
    """Construct the LSST survey scheduler.

    Returns
    -------
    nside : `int`
        Healpix map resolution.
    scheduler : `rubin_scheduler.scheduler.scheduler.CoreScheduler`
        Feature based scheduler.
    """
    nside = 32
    science_program = "BLOCK-430"

    band_to_filter = {
        "u": "u_24",
        "g": "g_6",
        "r": "r_57",
        "i": "i_39",
        "z": "z_20",
        "y": "y_10",
    }

    exptime = 30
    u_exptime = 38

    # Standard mask parameters - constraints on all survey pointings
    # Generally shadow_minutes value is set by the survey, but can
    # be set here as well (will be overwritten if too short for survey).
    standard_mask_params = {
        "nside": nside,
        "wind_speed_maximum": 40,
        "min_alt": 20,
        "max_alt": 86.5,
        "shadow_minutes": 2,
        "apply_time_limited_shadow": True,
        "time_to_sunrise": 3.0,
        "min_az_sunrise": 150,
        "max_az_sunrise": 250,
    }

    # Camera rot limits for the surveys
    # The CAMERA_ROT_LIMITS are sent to the queue manager and will
    # be enforced at last minute.
    # We can also set a smaller range here, which could sometimes be exceeded.
    camera_rot_limits = CAMERA_ROT_LIMITS

    # General parameters for standard pairs
    pair_time = 33

    # Adjust these as the expected timing updates.
    # -- sets the expected time and number of pointings in a 'blob'.
    blob_survey_params = {
        "slew_approx": 8,
        "band_change_approx": 200.0,
        "read_approx": 3.07,
        "flush_time": 30.0,
        "smoothing_kernel": None,
        "seed": 42,
        "twilight_scale": True,
        "check_scheduled": True,
    }

    # Parameters for near-sun twilight microsurvey
    ei_night_pattern = (True, False, False, False)  # once every 4 nights
    ei_bands = "riz"  # Bands to use for earth interior observations.
    ei_repeat = 4  # Number of times to repeat earth interior observations
    ei_am = 2.5  # Earth interior airmass limit
    ei_elong_req = 45.0  # Solar elongation required for inner solar system
    ei_area_req = 0.0  # Sky area required before attempting inner solar system

    # Generate footprint over the sky
    footprints, template_fp, footprint_mask = lsst_footprints.get_footprints(
        nside=nside,
        bandpasses=("u", "g", "r", "i", "z", "y"),
    )

    # Set up the ToO Surveys
    too_detailers = []
    too_detailers.append(
        detailers.CameraRotDetailer(
            min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
        )
    )
    too_detailers.append(detailers.BandSortDetailer())
    too_detailers.append(detailers.LabelRegionsAndDDFs())

    toos = too_surveys.gen_too_surveys(
        nside=nside,
        detailer_list=too_detailers,
        too_footprint=footprint_mask,
        masks=lsst_surveys.standard_masks(**standard_mask_params),
        science_program=science_program,
    )

    # Set up DDF survey:
    # Parameters for DDF and Roman  dithers
    camera_ddf_rot_limit = 55  # Rotator limit for DDF (degrees) .. 75
    camera_ddf_rot_per_visit = 3.0  # small rotation per visit (degrees) .. 3
    max_dither = 0.2  # Max radial dither for DDF (degrees)
    per_night = False  # Dither DDF per night (True) or per visit (False)

    detailer_list = [
        detailers.CameraSmallRotPerObservationListDetailer(
            min_rot=-camera_ddf_rot_limit,
            max_rot=camera_ddf_rot_limit,
            per_visit_rot=camera_ddf_rot_per_visit,
        ),
        detailers.SplitDetailer(
            detailers.DitherDetailer(per_night=per_night, max_dither=max_dither),
            detailers.EuclidDitherDetailer(per_night=per_night),
        ),
        detailers.BandNexp(bandname="u", nexp=1, exptime=u_exptime),
        detailers.BandSortDetailer(),
        detailers.LabelRegionsAndDDFs(),
        detailers.TruncatePreTwiDetailer(),
    ]

    # For the DDFs, keep the standard max alt (no slew),
    # but modify the shadow minutes.
    ddf_mask_params = copy.deepcopy(standard_mask_params)
    ddf_mask_params["shadow_minutes"] = 30
    ddf_masks = lsst_surveys.standard_masks(**ddf_mask_params)

    ddf_ignore = [
        "blob",
        "pair",
        "long",
        "greedy",
        "templates",
        "twilight",
        "ToO",
        "DD:RGES",
    ]

    ddfs = [
        ScriptedSurvey(
            ddf_masks,
            nside=nside,
            detailers=detailer_list,
            survey_name="deep drilling",
            before_twi_check=False,
            ignore_obs=ddf_ignore,
        )
    ]
    # This hash is provided by the script that
    # generates the pre-computed obs_array.
    # Execute ts_config_scheduler/Scheduler/ddf_gen/lsst_ddf_gen_block_*.py
    # and paste the provided value here.
    expected_hex_digest = "862b760"
    obs_array = read_ddf_obs_array(expected_hex_digest)
    ddfs[0].set_script(obs_array)

    # Define Roman scripted surveys
    roman_ignore = [
        "blob",
        "pair",
        "long",
        "greedy",
        "templates",
        "twilight",
        "ToO",
        "DD:COSMOS",
        "DD:ECDFS",
        "DD:EDFS",
        "DD:ELIASS",
        "DD:XMM",
    ]
    roman_micro = [
        roman_surveys.gen_roman_on_season(
            nside=nside,
            max_dither=max_dither,
            per_night=per_night,
            camera_ddf_rot_limit=camera_ddf_rot_limit,
            camera_ddf_rot_per_visit=camera_ddf_rot_per_visit,
            exptimes=exptime,
            science_program=science_program,
            standard_mask_params=standard_mask_params,
            ignore_obs=roman_ignore,
        ),
        roman_surveys.gen_roman_off_season(
            nside=nside,
            max_dither=max_dither,
            per_night=per_night,
            camera_ddf_rot_limit=camera_ddf_rot_limit,
            camera_ddf_rot_per_visit=camera_ddf_rot_per_visit,
            exptimes=exptime,
            science_program=science_program,
            standard_mask_params=standard_mask_params,
            ignore_obs=roman_ignore,
        ),
    ]

    # Set up long gaps (triplets) survey.
    gaps_night_pattern = (False, True, False, False)
    long_gaps = lsst_surveys.gen_long_gaps_survey(
        footprints=footprints,
        nside=nside,
        camera_rot_limits=camera_rot_limits,
        exptime=exptime,
        u_exptime=u_exptime,
        pair_time=pair_time,
        night_pattern=gaps_night_pattern,
        science_program=science_program,
        blob_survey_params=blob_survey_params,
        standard_mask_params=standard_mask_params,
    )

    # Define the standard pairs during the night survey.
    # The shadow minutes limit will be modified in the survey by the pair time.
    blobs = lsst_surveys.generate_blobs(
        footprints=footprints,
        nside=nside,
        camera_rot_limits=camera_rot_limits,
        exptime=exptime,
        u_exptime=u_exptime,
        pair_time=pair_time,
        science_program=science_program,
        blob_survey_params=blob_survey_params,
        standard_mask_params=standard_mask_params,
    )

    # Define the near-sun twilight microsurvey.
    neo_micro = lsst_surveys.generate_twilight_near_sun(
        nside=nside,
        night_pattern=ei_night_pattern,
        max_airmass=ei_am,
        camera_rot_limits=camera_rot_limits,
        footprint_mask=footprint_mask,
        min_area=ei_area_req,
        bands=ei_bands,
        n_repeat=ei_repeat,
        max_elong=ei_elong_req,
        science_program=science_program,
        standard_mask_params=standard_mask_params,
    )

    # Define the greedy surveys (single-visit per call).
    greedy = lsst_surveys.gen_greedy_surveys(
        nside=nside,
        camera_rot_limits=camera_rot_limits,
        exptime=exptime,
        u_exptime=u_exptime,
        footprints=footprints,
        science_program=science_program,
        standard_mask_params=standard_mask_params,
    )

    # In all of the pairs-related surveys with u visits,
    # ensure the u-band exptime detailer comes last, after TakeAsPairs.
    # Remove this after fixed in ts_fbs_utils.
    for tier in [blobs, short_blobs]:
        for survey in tier:
            d = survey.detailers
            dm = [di for di in d if isinstance(di, detailers.BandNexp)]
            if len(dm) > 0:
                d.remove(dm[0])
                d.insert(-1, dm[0])
                survey.detailers = d

    # Arrange the surveys in tiers.
    surveys = [
        toos,
        roman_micro,
        ddfs,
        long_gaps,
        blobs,
        neo_micro,
        greedy,
    ]

    # Combine into CoreScheduler
    scheduler = CoreScheduler(
        surveys,
        nside=nside,
        survey_start_mjd=SURVEY_START_MJD,
        band_to_filter=band_to_filter,
    )

    return nside, scheduler


if __name__ == "config":
    nside, scheduler = get_scheduler()
