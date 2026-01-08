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


from pathlib import Path

import numpy as np
import rubin_scheduler.scheduler.detailers as detailers
from astropy.time import Time
from lsst.ts.fbs.utils.maintel.lsst_surveys import safety_masks
from rubin_scheduler.scheduler.schedulers import CoreScheduler
from rubin_scheduler.scheduler.surveys import ScriptedSurvey
from rubin_scheduler.scheduler.utils import ScheduledObservationArray


def get_scheduler() -> tuple[int, CoreScheduler]:
    """Construct a scheduler for straylight templates creation.

    The parameters are not accessible when calling as 'config'.

    Returns
    -------
    nside : `int`
        Healpix map resolution.
    scheduler : `rubin_scheduler.scheduler.scheduler.CoreScheduler`
        Feature based scheduler.
    """
    nside = 32

    band_to_filter = {
        "u": "u_24",
        "g": "g_6",
        "r": "r_57",
        "i": "i_39",
        "z": "z_20",
        "y": "y_10",
    }

    band = "i"

    safety_mask_params = {
        "nside": nside,
        "wind_speed_maximum": 40,
        "apply_time_limited_shadow": False,
        "min_alt": 20,
        "shadow_minutes": 0,
    }

    # Get path for template pointing information
    config_dir = Path(__file__).parent
    template_pointing_file = Path.joinpath(
        config_dir, "small_data", "straylight_template_pointings.dat"
    )
    if not Path.exists(template_pointing_file):
        raise ValueError(
            f"Expected pointing file does not exist at {template_pointing_file}"
        )

    template_pointings = np.loadtxt(template_pointing_file, delimiter=",")
    tt = template_pointings.swapaxes(0, 1)
    ra_template = tt[1]
    dec_template = tt[2]
    ra_cen = np.mean(ra_template)
    dec_cen = np.mean(dec_template)

    template_sched_obs = ScheduledObservationArray(n=len(template_pointings))
    template_sched_obs["RA"] = np.radians(ra_template)
    template_sched_obs["dec"] = np.radians(dec_template)
    # To get the scheduled observation in the same order as the template_list
    fudge = np.arange(0, len(template_pointings), 1) * (10 / 60 / 60 / 24)
    template_sched_obs["mjd"] = Time("2025-12-06T12:00:00").mjd + fudge
    template_sched_obs["flush_by_mjd"] = Time("2025-12-06T12:00:00").mjd + 60
    template_sched_obs["mjd_tol"] = 1
    template_sched_obs["dist_tol"] = np.radians(2)
    template_sched_obs["HA_min"] = 18
    template_sched_obs["HA_max"] = 6
    template_sched_obs["alt_min"] = np.radians(20)
    template_sched_obs["alt_max"] = np.radians(85)
    template_sched_obs["sun_alt_max"] = np.radians(-12)
    template_sched_obs["moon_min_distance"] = np.radians(20)
    template_sched_obs["exptime"] = 30
    template_sched_obs["band"] = band
    template_sched_obs["nexp"] = 1
    template_sched_obs["scheduler_note"] = "straylight_templates_1"
    template_sched_obs["target_name"] = (
        f"straylight_templates_{ra_cen:.1f}_{dec_cen:.1f}"
    )
    template_sched_obs["observation_reason"] = "straylight_templates"
    template_sched_obs["science_program"] = "BLOCK-T506"

    detailer_list = [detailers.CameraRotDetailer(max_rot=10, min_rot=-10, dither="all")]

    basis_functions = safety_masks(**safety_mask_params)

    template_survey = ScriptedSurvey(
        basis_functions=basis_functions,
        nside=nside,
        detailers=detailer_list,
        return_n_limit=50,
        survey_name=f"StrayLightTemplate_{ra_cen:.1f}_{dec_cen:.1f}",
    )
    template_survey.set_script(template_sched_obs)

    # Add followup test visits

    test_pointing_file = Path.joinpath(
        config_dir, "small_data", "straylight_test_pointings.dat"
    )
    if not Path.exists(test_pointing_file):
        raise ValueError(
            f"Expected pointing file does not exist at {test_pointing_file}"
        )

    test_pointings = np.loadtxt(test_pointing_file, delimiter=",")
    tt = test_pointings.swapaxes(0, 1)
    ra_test = tt[1]
    dec_test = tt[2]
    ra_cen = np.mean(ra_test)
    dec_cen = np.mean(dec_test)

    test_sched_obs = ScheduledObservationArray(n=len(test_pointings))
    test_sched_obs["RA"] = np.radians(ra_test)
    test_sched_obs["dec"] = np.radians(dec_test)
    # To get the scheduled observation in the same order as the template_list
    fudge = np.arange(0, len(test_pointings), 1) * (10 / 60 / 60 / 24)
    test_sched_obs["mjd"] = Time("2025-12-06T12:00:00").mjd + fudge
    test_sched_obs["flush_by_mjd"] = Time("2025-12-06T12:00:00").mjd + 60
    test_sched_obs["mjd_tol"] = 1
    test_sched_obs["dist_tol"] = np.radians(2)
    test_sched_obs["HA_min"] = 18
    test_sched_obs["HA_max"] = 6
    test_sched_obs["alt_min"] = np.radians(20)
    test_sched_obs["alt_max"] = np.radians(85)
    test_sched_obs["sun_alt_max"] = np.radians(-12)
    test_sched_obs["moon_min_distance"] = np.radians(20)
    test_sched_obs["exptime"] = 30
    test_sched_obs["band"] = band
    test_sched_obs["nexp"] = 1
    test_sched_obs["scheduler_note"] = "straylight_test_1"
    test_sched_obs["target_name"] = f"straylight_test_{ra_cen:.1f}_{dec_cen:.1f}"
    test_sched_obs["observation_reason"] = "straylight_test"
    test_sched_obs["science_program"] = "BLOCK-T506"

    detailer_list = [detailers.CameraRotDetailer(max_rot=10, min_rot=-10, dither="all")]

    basis_functions = safety_masks(**safety_mask_params)

    test_survey = ScriptedSurvey(
        basis_functions=basis_functions,
        nside=nside,
        detailers=detailer_list,
        return_n_limit=50,
        survey_name=f"StrayLightTest_{ra_cen:.1f}_{dec_cen:.1f}",
    )
    test_survey.set_script(test_sched_obs)

    # Combine into CoreScheduler
    scheduler = CoreScheduler(
        [[template_survey], [test_survey]],
        nside=nside,
        band_to_filter=band_to_filter,
    )

    return nside, scheduler


if __name__ == "config":
    nside, scheduler = get_scheduler()
