# This file is part of ts_config_ocs.
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

import numpy as np
from astropy.time import Time
from lsst.ts.fbs.utils.maintel.lsst_surveys import safety_masks
from rubin_scheduler.scheduler.basis_functions import (
    BalanceVisits,
    RewardNObsSequence,
)
from rubin_scheduler.scheduler.detailers import (
    AltAz2RaDecDetailer,
    CopyValueDetailer,
    Rottep2RotspDesiredDetailer,
)
from rubin_scheduler.scheduler.schedulers import CoreScheduler
from rubin_scheduler.scheduler.surveys import FieldAltAzSurvey


def get_scheduler():
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
        "wind_speed_maximum": 20,
        "shadow_minutes": 0,
        "apply_time_limited_shadow": False,
        "time_to_sunrise": 3.0,
        "min_az_sunrise": 144,
        "max_az_sunrise": 255,
    }

    block_name = "BLOCK-T706"

    # For these alt-az-rotTel tests, should specify only a single visit
    # at a time, so sequence is 1 visit long. If multiple bands are
    # needed, the sequence will have to be longer but should be as short
    # as possible to avoid rotator and alt/az <-> RA/Dec drift.
    sequence = ["i"]
    nvisits = {"u": 1, "g": 1, "r": 1, "i": 1, "z": 1, "y": 1}
    exptimes = {"u": 38, "g": 30, "r": 30, "i": 30, "z": 30, "y": 30}
    # For each survey (alt/az/rotTelPos combo) how many visits each time
    # before going on to the next target?
    nvis_per_cycle = 100

    # Using the current time in the note
    # and requiring this to reset daily means that
    # this does require the FBS to reconfigured, from cold-start daily.
    day_obs = int(
        Time(int(Time.now().mjd - 0.5), format="mjd", scale="utc")
        .iso[0:10]
        .replace("-", "")
    )
    scheduler_root = f"{block_name} {day_obs}"

    # Set up target information
    target_dict = {}

    # Setup alt, az, rotTelPos values
    _az_values = [0.0]
    _alt_values = [70.0]
    _rot_values = [0.0, 60.0]

    # List of alt, az, rotTelPos to use -- in DEGREES.
    # Note order is important and will set the order of observation.
    for alt in _alt_values:
        for az in _az_values:
            for rotTelPos in _rot_values:
                name = f"alt:{alt:.1f} az:{az:.1f} rotTel:{rotTelPos:.0f}"
                target_dict[name] = {}
                target_dict[name]["note"] = f"{scheduler_root} {name}"
                target_dict[name]["alt"] = alt
                target_dict[name]["az"] = az
                target_dict[name]["rotTelPos"] = rotTelPos

    n_pointings = len(target_dict)

    detailers = [
        AltAz2RaDecDetailer(),
        Rottep2RotspDesiredDetailer(),
        CopyValueDetailer(source="rotSkyPos_desired", destination="rotSkyPos"),
    ]

    safety_masks_basis_functions = safety_masks(**safety_mask_params)
    # Remove avoid wind basis function
    safety_masks_basis_functions.pop(2)

    survey_lists = []
    for target in target_dict:
        tt = target_dict[target]
        survey = FieldAltAzSurvey(
            basis_functions=[
                RewardNObsSequence(
                    n_obs_survey=nvis_per_cycle,
                    note_survey=tt["note"],
                ),
                BalanceVisits(
                    nobs_reference=nvis_per_cycle * n_pointings,
                    note_survey=tt["note"],
                    note_interest=f"AOS {scheduler_root}",
                    nside=nside,
                ),
            ]
            + safety_masks_basis_functions,
            alt=tt["alt"],
            az=tt["az"],
            sequence=sequence,
            nvisits=nvisits,
            exptimes=exptimes,
            ignore_obs=None,
            survey_name=tt["note"],
            target_name=tt["note"],
            science_program=block_name,
            observation_reason="fbs driven aos stability test",
            scheduler_note=tt["note"],
            nside=nside,
            flush_pad=30.0,
            detailers=detailers,
        )
        survey.observations["rotTelPos"] = np.radians(tt["rotTelPos"])
        survey_lists.append(survey)

    return nside, CoreScheduler(
        [survey_lists],
        nside=nside,
        band_to_filter=band_to_filter,
    )


if __name__ == "config":
    nside, scheduler = get_scheduler()
