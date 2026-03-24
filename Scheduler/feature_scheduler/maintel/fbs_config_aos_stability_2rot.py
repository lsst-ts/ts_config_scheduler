# This file is part of ts_config_scheduler.
#
# Developed for the Vera Rubin Observatory Telescope and Site System.
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

from lsst.ts.fbs.utils.maintel.stability_surveys import gen_az_el_rot_stability_survey
from rubin_scheduler.scheduler.schedulers import CoreScheduler


def get_scheduler() -> tuple[int, CoreScheduler]:
    """Construct feature based scheduler.

    Returns
    -------
    nside : int
        Healpix map resolution.
    scheduler : Core_scheduler
        Feature based scheduler.
    """

    nside = 32

    # Mapping from band to filter from
    # obs_lsst/python/lsst/obs/lsst/filters.py
    band_to_filter = {
        "u": "u_24",
        "g": "g_6",
        "r": "r_57",
        "i": "i_39",
        "z": "z_20",
        "y": "y_10",
    }

    safety_mask_params = {
        "nside": nside,
        "wind_speed_maximum": None,
        "shadow_minutes": 0,
        "apply_time_limited_shadow": False,
        "time_to_sunrise": 3.0,
        "min_az_sunrise": 144,
        "max_az_sunrise": 255,
    }

    az_values = [0.0, 0.0]
    el_values = [70.0, 70.0]
    rot_values = [0.0, 60.0]

    block_name = "BLOCK-T706"

    sequence = ["i"]
    nvis_per_cycle = 100

    survey_lists = gen_az_el_rot_stability_survey(
        az_values=az_values,
        el_values=el_values,
        rot_values=rot_values,
        science_program=block_name,
        observation_reason="fbs driven aos stability test",
        sequence=sequence,
        nvis_per_cycle=nvis_per_cycle,
        safety_mask_params=safety_mask_params,
    )

    return nside, CoreScheduler(
        survey_lists,
        nside=nside,
        band_to_filter=band_to_filter,
    )


if __name__ == "config":
    nside, scheduler = get_scheduler()
