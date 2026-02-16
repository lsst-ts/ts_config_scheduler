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

import copy
import hashlib
import os
import pathlib

import lsst.ts.fbs.utils.maintel.lsst_surveys as lsst_surveys
import lsst.ts.fbs.utils.maintel.roman_surveys as roman_surveys
import lsst.ts.fbs.utils.maintel.too_surveys as too_surveys
import numpy as np
import rubin_scheduler.scheduler.detailers as detailers
from lsst.ts.fbs.utils.maintel.lsst_surveys import safety_masks
from rubin_scheduler.data import get_data_dir
from rubin_scheduler.scheduler.schedulers import CoreScheduler
from rubin_scheduler.scheduler.surveys import ScriptedSurvey
from rubin_scheduler.scheduler.utils import (
    CurrentAreaMap,
    Footprint,
    ScheduledObservationArray,
    make_rolling_footprints,
)
from rubin_scheduler.site_models import Almanac
from rubin_scheduler.utils import SURVEY_START_MJD


def first_alert_obs(science_program: str) -> ScheduledObservationArray:
    # First alerts:
    # Observe EDFS a and b first, then COSMOS, then M49
    # Test with ~hour on each field over week prior to event, then
    # full night no earlier than Feb 24 61095.5
    # MJD 61087.5 to 61094.5
    first_alert_day = 61095.5
    alert_fields = {
        "EDFS_a": (58.9, -49.32),  # griz
        "EDFS_b": (63.6, -47.6),  # griz
        "COSMOS": (150.108, 2.2336),  # ugrizy
        "M49": (186.2, 7.0),  # ugri
    }
    alert_field_bands = {
        "EDFS_a": "griz",
        "EDFS_b": "griz",
        "COSMOS": "ugrizy",
        "M49": "ugri",
    }
    # Typical sequence
    nvis = {"u": 10, "g": 15, "r": 15, "i": 15, "z": 15, "y": 15}
    obs_params = {
        "flush_length": 1.0,  # days
        "mjd_tol": 15 / 60 / 24.0,  # minutes
        "dist_tol": np.radians(3.0),
        "alt_min": np.radians(25.0),
        "alt_max": np.radians(85.0),
        "sun_alt_max": np.radians(-12.0),
        "HA_min": 21.0,
        "HA_max": 3.0,
        "moon_min_distance": np.radians(25.0),
    }
    expt = {"u": 38, "g": 30, "r": 30, "i": 30, "z": 30, "y": 30}
    thirty_seconds = 30 / 60 / 60 / 24
    # Order by socket
    band_order = "igzr"

    sched_obs = []
    # visits for the week prior to first alerts, plus 2 days
    mjds = np.arange(61087.5, first_alert_day + 2, 1)
    for mjd in mjds:
        c_mjd = mjd
        for name in ["EDFS_a", "COSMOS", "M49"]:
            ra = alert_fields[name][0]
            dec = alert_fields[name][1]
            # Make the metadata for our real DDFs be different than M49
            if name == "M49":
                target_name = "field_" + name.lower()
            else:
                target_name = "ddf_" + name.lower()
            # Run different number of cycles through each field
            n_cycles = 1
            n_per_cycle = 1
            if name == "EDFS_a":
                n_per_cycle = 1
                n_cycles = 1
            elif name == "M49":
                n_per_cycle = 1
                n_cycles = 2
            elif name == "COSMOS":
                n_per_cycle = 1
                n_cycles = 2
            if mjd >= first_alert_day:
                if name == "EDFS_a":
                    n_per_cycle = 1.6
                    n_cycles = 1
                elif name == "COSMOS":
                    n_per_cycle = 2
                    n_cycles = 3
                elif name == "M49":
                    n_per_cycle = 3
                    n_cycles = 3

            for cycle in range(n_cycles):
                for band in band_order:
                    if band in alert_field_bands[name]:
                        n_vis_band = int(nvis[band] * n_per_cycle)
                        obs = ScheduledObservationArray(n=n_vis_band)
                        obs["RA"] = np.radians(ra)
                        obs["dec"] = np.radians(dec)
                        # Set like this to ensure desired order of visits
                        # The mjd is not intended to be precise at all.
                        obs["mjd"] = c_mjd + np.arange(0, n_vis_band) * thirty_seconds
                        obs["flush_by_mjd"] = mjd + obs_params["flush_length"]
                        obs["exptime"] = expt[band]
                        obs["band"] = band
                        obs["nexp"] = 1
                        obs["scheduler_note"] = f"Alert:{name}"
                        obs["target_name"] = target_name
                        obs["science_program"] = science_program
                        obs["observation_reason"] = "alert_" + name.lower()
                        obs["mjd_tol"] = obs_params["mjd_tol"]
                        obs["dist_tol"] = obs_params["dist_tol"]
                        obs["HA_min"] = obs_params["HA_min"]
                        obs["HA_max"] = obs_params["HA_max"]
                        obs["alt_min"] = obs_params["alt_min"]
                        obs["alt_max"] = obs_params["alt_max"]
                        obs["sun_alt_max"] = obs_params["sun_alt_max"]
                        obs["moon_min_distance"] = obs_params["moon_min_distance"]
                        # I want the different fields to go sequentially
                        c_mjd += (n_vis_band + 2) * thirty_seconds
                        if name == "EDFS_a":
                            # Copy EDFS_a to EDFS_b, so they execute in tandem
                            obs_b = copy.deepcopy(obs)
                            name_b = "EDFS_b"
                            obs_b["RA"] = np.radians(alert_fields[name_b][0])
                            obs_b["dec"] = np.radians(alert_fields[name_b][1])
                            obs_b["scheduler_note"] = f"Alert:{name_b}"
                            obs_b["target_name"] = "ddf_" + name_b.lower()
                            obs_b["observation_reason"] = "alert_" + name_b.lower()
                            obs = np.concat([obs, obs_b])
                            obs = np.sort(obs, order="mjd")
                        sched_obs.append(obs)
    all_scheduled_obs = np.concatenate(sched_obs)
    # But now all_scheduled_obs is just an ndarray with the right fields
    # So .. convert it back.
    sched_obs_all = ScheduledObservationArray(n=len(all_scheduled_obs))
    for key in all_scheduled_obs.dtype.names:
        sched_obs_all[key] = all_scheduled_obs[key]
    return sched_obs_all


def get_scheduler() -> tuple[int, CoreScheduler]:
    """Construct the LSST survey scheduler.

    The parameters are not accessible when calling as 'config'.

    Returns
    -------
    nside : `int`
        Healpix map resolution.
    scheduler : `rubin_scheduler.scheduler.scheduler.CoreScheduler`
        Feature based scheduler.
    """
    nside = 32
    science_program = "BLOCK-421"
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
    # survey_start_mjd = 61100.5

    # Safety mask parameters - constraints on all survey pointings
    # Generally shadow_minutes value is set by the survey, but can
    # be set here as well (will be overwritten if too short for survey).
    safety_mask_params = {
        "nside": nside,
        "wind_speed_maximum": 40,
        "apply_time_limited_shadow": False,
        "time_to_sunrise": 3.0,
        "min_az_sunrise": 150,
        "max_az_sunrise": 250,
    }

    safety_mask_params_ddf = copy.deepcopy(safety_mask_params)
    safety_mask_params_ddf["shadow_minutes"] = 30

    # General parameters for standard pairs (-80/80 default)
    camera_rot_limits = (-60.0, 60.0)
    pair_time = 33
    # Adjust these as the expected timing updates.
    # -- sets the expected time and number of pointings in a 'blob'.
    blob_survey_params = {
        "slew_approx": 8,
        "band_change_approx": 320.0,
        "read_approx": 3.07,
        "flush_time": 30.0,
        "smoothing_kernel": None,
        "nside": nside,
        "seed": 42,
        "dither": "night",
        "twilight_scale": True,
    }
    # Seeing (FWHM in ") max for template
    fwhm_template_max = 1.4

    # Parameters for rolling cadence footprint definition
    nslice = 2  # N slices for rolling
    rolling_scale = 0.9  # Strength of rolling
    rolling_uniform = True  # Should we use the uniform rolling flag

    # Parameters for long-gaps survey
    nights_off = 3  # For long gaps

    # Parameters for near-sun twilight microsurvey
    ei_night_pattern = 4  # see pattern_dict below
    ei_bands = "riz"  # Bands to use for earth interior observations.
    ei_repeat = 4  # Number of times to repeat earth interior observations
    ei_am = 2.5  # Earth interior airmass limit
    ei_elong_req = 45.0  # Solar elongation required for inner solar system
    ei_area_req = 0.0  # Sky area required before attempting inner solar system

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

    # Define the long-gaps (triplets) survey.
    gaps_night_pattern = [True] + [False] * nights_off
    long_gaps = lsst_surveys.gen_long_gaps_survey(
        footprints=footprints,
        nside=nside,
        camera_rot_limits=camera_rot_limits,
        exptime=exptime,
        nexp=nexp,
        u_exptime=u_exptime,
        u_nexp=u_nexp,
        pair_time=pair_time,
        night_pattern=gaps_night_pattern,
        science_program=science_program,
        blob_survey_params=blob_survey_params,
        safety_mask_params=safety_mask_params,
    )

    # This hash is provided by the script that
    # generates the pre-computed data. Execute it and paste
    # the provided value here.
    expected_hex_digest = "70645ab"
    pre_comp_file = (
        pathlib.Path(get_data_dir())
        / "scheduler"
        / f"ts_ddf_array_{expected_hex_digest}.npz"
    )
    if os.path.exists(pre_comp_file):
        loaded = np.load(pre_comp_file, allow_pickle=True)
        hash_object = hashlib.sha256()
        hash_object.update(loaded["hash_digest"])
        hex_digest = hash_object.hexdigest()[:7]
        if hex_digest == expected_hex_digest:
            obs_array_loaded = loaded["obs_array"]
            obs_array = ScheduledObservationArray(obs_array_loaded.size)
            for key in obs_array_loaded.dtype.names:
                obs_array[key] = obs_array_loaded[key]
        else:
            raise RuntimeError(
                f"Provided hash {expected_hex_digest} does not match loaded file hash {hex_digest}. "
                "Reach out for support so they can help you generate the correct file."
            )
        loaded.close()
    else:
        raise RuntimeError(
            f"Pre-computed DDF files {pre_comp_file} not available. "
            "Reach out for support so they can help execute the script "
            "that will generate this file before proceeding."
        )

    # Parameters for  DDF dithers
    camera_ddf_rot_limit = 55  # Rotator limit for DDF (degrees) .. 75
    camera_ddf_rot_per_visit = 3.0  # small rotation per visit (degrees) .. 3
    max_dither = 0.2  # Max radial dither for DDF (degrees)
    per_night = False  # Dither DDF per night (True) or per visit (False)

    band_expt = {
        "u": u_exptime,
        "g": exptime,
        "r": exptime,
        "i": exptime,
        "z": exptime,
        "y": exptime,
    }
    band_nexp = {"u": u_nexp, "g": nexp, "r": nexp, "i": nexp, "z": nexp, "y": nexp}
    u_nexp = band_nexp["u"]
    u_exptime = band_expt["u"]

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
        detailers.BandNexp(bandname="u", nexp=u_nexp, exptime=u_exptime),
        detailers.BandSortDetailer(),
        detailers.LabelRegionsAndDDFs(),
        detailers.TruncatePreTwiDetailer(),
    ]

    ddfs = [
        ScriptedSurvey(
            safety_masks(**safety_mask_params_ddf),
            nside=nside,
            detailers=detailer_list,
            survey_name="deep drilling",
            before_twi_check=False,
        )
    ]
    ddfs[0].set_script(obs_array)

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
        detailers.BandNexp(bandname="u", nexp=u_nexp, exptime=u_exptime),
        # Purposefully excluded bandsort detailer
        detailers.LabelRegionsAndDDFs(),
        detailers.TruncatePreTwiDetailer(),
    ]
    first_alerts = [
        ScriptedSurvey(
            safety_masks(**safety_mask_params_ddf),
            nside=nside,
            detailers=detailer_list,
            survey_name="first alerts",
            before_twi_check=False,
            return_n_limit=15,
        )
    ]
    alert_obs_array = first_alert_obs(science_program)
    first_alerts[0].set_script(alert_obs_array)

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

    # Define the near-sun twilight microsurvey
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
        safety_mask_params=safety_mask_params,
    )

    # Define the alternate twilight (and other short time period)
    # short 15minute pairs
    short_blobs = lsst_surveys.generate_short_blobs(
        footprints=footprints,
        nside=nside,
        camera_rot_limits=camera_rot_limits,
        exptime=exptime,
        nexp=nexp,
        pair_time=15.0,
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

    # Define Roman scripted surveys
    roman_micro = [
        roman_surveys.gen_roman_on_season(
            nside=nside,
            max_dither=max_dither,
            per_night=per_night,
            camera_ddf_rot_limit=camera_ddf_rot_limit,
            camera_ddf_rot_per_visit=camera_ddf_rot_per_visit,
            exptimes=exptime,
            nexps=nexp,
            science_program=science_program,
            safety_mask_params=safety_mask_params,
        ),
        roman_surveys.gen_roman_off_season(
            nside=nside,
            max_dither=max_dither,
            per_night=per_night,
            camera_ddf_rot_limit=camera_ddf_rot_limit,
            camera_ddf_rot_per_visit=camera_ddf_rot_per_visit,
            exptimes=exptime,
            nexps=nexp,
            science_program=science_program,
            safety_mask_params=safety_mask_params,
        ),
    ]

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

    # Define ToO surveys
    too_detailers = []
    too_detailers.append(
        detailers.CameraRotDetailer(
            min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
        )
    )
    too_detailers.append(detailers.LabelRegionsAndDDFs())
    # Let's make a footprint to follow up ToO events
    too_footprint = np.where(footprints_hp["r"] > 0, 1.0, np.nan)

    toos = too_surveys.gen_too_surveys(
        nside=nside,
        detailer_list=too_detailers,
        too_footprint=too_footprint,
        # Re-evaluate this when on-sky
        split_long=False,
        n_snaps=nexp,
        science_program=science_program,
        safety_mask_params=safety_mask_params,
    )

    # Arrange the surveys in tiers.
    surveys = [
        toos,
        first_alerts,
        roman_micro,
        ddfs,
        template_surveys,
        long_gaps,
        blobs,
        short_blobs,
        neo_micro,
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
