import logging

import lsst.ts.fbs.utils.maintel.lsst_surveys as lsst_surveys
import numpy as np
from rubin_scheduler.scheduler.schedulers import CoreScheduler
from rubin_scheduler.scheduler.utils import (
    CurrentAreaMap,
    Footprint,
    make_rolling_footprints,
)
from rubin_scheduler.site_models import Almanac
from rubin_scheduler.utils import SURVEY_START_MJD

__all__ = ("get_scheduler",)

logger = logging.getLogger(__name__)


def get_scheduler() -> tuple[int, CoreScheduler]:
    """Construct a partial LSST survey scheduler.

    Tries to only retain surveys which would incorporate the slewtime
    map into their basis functions, to better restrict visits to
    few or no slews greater than 9 degrees in elevation.

    Returns
    -------
    nside : `int`
        Healpix map resolution.
    scheduler : `rubin_scheduler.scheduler.scheduler.CoreScheduler`
        Feature based scheduler.
    """
    nside = 32
    science_program = "BLOCK-407"
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
        "time_to_sunrise": 3.0,
        "min_az_sunrise": 120,
        "max_az_sunrise": 290,
    }

    # General parameters for standard pairs (-80/80 default)
    camera_rot_limits = (-60.0, 60.0)
    pair_time = 33
    # Adjust these as the expected timing updates.
    # -- sets the expected time and number of pointings in a 'blob'.
    blob_survey_params = {
        "slew_approx": 8,
        "band_change_approx": 140.0,
        "read_approx": 3.07,
        "flush_time": 30.0,
        "smoothing_kernel": None,
        "nside": nside,
        "seed": 42,
        "dither": "night",
        "twilight_scale": True,
    }
    # Seeing (FWHM in ") max for template
    fwhm_template_max = 1.3

    # Parameters for rolling cadence footprint definition
    nslice = 2  # N slices for rolling
    rolling_scale = 0.9  # Strength of rolling
    rolling_uniform = True  # Should we use the uniform rolling flag

    # Parameters for near-sun twilight microsurvey
    ei_night_pattern = 4  # see pattern_dict below

    # Mapping for night_pattern for near-sun twilight / twi_blob surveys.
    pattern_dict = {
        1: [True],
        2: [True, False],
        3: [True, False, False],
        4: [True, False, False, False],
        5: [True, True, True, True, False, False, False, False],
        6: [True, True, True, False, False, False, False],
        7: [True, True, False, False, False, False],
    }
    ei_night_pattern = pattern_dict[ei_night_pattern]
    reverse_ei_night_pattern = [not val for val in ei_night_pattern]

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

    # Define the alternate twilight (and other short time period)
    # short 15minute pairs
    twi_blobs = lsst_surveys.generate_twi_blobs(
        footprints=footprints,
        nside=nside,
        exptime=exptime,
        nexp=nexp,
        repeat_weight=0,
        night_pattern=reverse_ei_night_pattern,
        science_program=science_program,
        blob_survey_params=blob_survey_params,
        safety_mask_params=safety_mask_params,
    )

    # Define the standard pairs during the night survey
    blobs = lsst_surveys.generate_blobs(
        footprints=footprints,
        nside=nside,
        camera_rot_limits=camera_rot_limits,
        exptime=exptime,
        nexp=nexp,
        u_exptime=u_exptime,
        u_nexp=u_nexp,
        pair_time=pair_time,
        survey_start=survey_start_mjd,
        science_program=science_program,
        blob_survey_params=blob_survey_params,
        safety_mask_params=safety_mask_params,
    )

    # Create template footprint.
    # Similar to rolling footprint but tracks visits separately
    # (only good seeing visits) and no rolling.
    template_fp = Footprint(survey_start_mjd, sun_ra_start, nside=nside)
    for key in footprints_hp_array.dtype.names:
        tmp_fp = np.where(footprints_hp_array[key] > 0, 1, np.nan)
        template_fp.set_footprint(key, tmp_fp)
    # Define template surveys
    template_surveys = lsst_surveys.gen_template_surveys(
        template_fp,
        nside=nside,
        seeing_fwhm_max=fwhm_template_max,
        camera_rot_limits=camera_rot_limits,
        exptime=exptime,
        nexp=nexp,
        u_exptime=u_exptime,
        u_nexp=u_nexp,
        n_obs_template={"u": 4, "g": 4, "r": 4, "i": 4, "z": 4, "y": 4},
        science_program=science_program,
        blob_survey_params=blob_survey_params,
        safety_mask_params=safety_mask_params,
    )

    # Arrange the surveys in tiers.
    surveys = [
        template_surveys,
        blobs,
        twi_blobs,
        greedy,
    ]

    # Combine into CoreScheduler
    scheduler = CoreScheduler(
        surveys,
        nside=nside,
        survey_start_mjd=survey_start_mjd,
        band_to_filter=band_to_filter,
    )

    logger.info(
        f"Configured {len(surveys)} tiers of surveys in the baseline configuration."
    )

    return (nside, scheduler)


if __name__ == "config":
    (nside, scheduler) = get_scheduler()


if __name__ == "__main__":
    # This doesn't really serve a purpose for this config, because no DDFs
    (nside, scheduler) = get_scheduler()
