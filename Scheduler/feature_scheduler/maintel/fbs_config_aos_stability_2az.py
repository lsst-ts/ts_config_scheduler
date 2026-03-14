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

    # Change the alt/az here if necessary.
    alt = 70.0
    az_1 = 0.0
    az_2 = -90.0

    target_name_1 = f"AOS alt:{alt:.1f} az:{az_1:.1f}"
    target_name_2 = f"AOS alt:{alt:.1f} az:{az_2:.1f}"
    block_name = "BLOCK-T703"

    sequence = ["i"]
    nvisits = {"u": 1, "g": 1, "r": 1, "i": 1, "z": 1, "y": 1}
    exptimes = {"u": 38, "g": 30, "r": 30, "i": 30, "z": 30, "y": 30}

    detailers = [
        AltAz2RaDecDetailer(),
        Rottep2RotspDesiredDetailer(),
        CopyValueDetailer(source="rotSkyPos_desired", destination="rotSkyPos"),
    ]

    safety_masks_basis_functions = safety_masks(**safety_mask_params)
    # Remove avoid wind basis function
    safety_masks_basis_functions.pop(2)

    aos_stability_survey_1 = FieldAltAzSurvey(
        basis_functions=[
            RewardNObsSequence(n_obs_survey=100, note_survey=target_name_1),
            BalanceVisits(
                nobs_reference=240,
                note_survey=target_name_1,
                note_interest="AOS alt",
                nside=nside,
            ),
        ]
        + safety_masks_basis_functions,
        alt=alt,
        az=az_1,
        sequence=sequence,
        nvisits=nvisits,
        exptimes=exptimes,
        ignore_obs=None,
        survey_name=target_name_1,
        target_name=target_name_1,
        science_program=block_name,
        observation_reason="fbs driven aos stability test",
        scheduler_note=target_name_1,
        nside=nside,
        flush_pad=30.0,
        detailers=detailers,
    )

    aos_stability_survey_2 = FieldAltAzSurvey(
        basis_functions=[
            RewardNObsSequence(n_obs_survey=100, note_survey=target_name_2),
            BalanceVisits(
                nobs_reference=240,
                note_survey=target_name_2,
                note_interest="AOS alt",
                nside=nside,
            ),
        ]
        + safety_masks_basis_functions,
        alt=alt,
        az=az_2,
        sequence=sequence,
        nvisits=nvisits,
        exptimes=exptimes,
        ignore_obs=None,
        survey_name=target_name_2,
        target_name=target_name_2,
        science_program=block_name,
        observation_reason="fbs driven aos stability test",
        scheduler_note=target_name_2,
        nside=nside,
        flush_pad=30.0,
        detailers=detailers,
    )

    survey_lists = [
        [
            aos_stability_survey_1,
            aos_stability_survey_2,
        ],
    ]

    return nside, CoreScheduler(
        survey_lists,
        nside=nside,
        band_to_filter=band_to_filter,
    )


if __name__ == "config":
    nside, scheduler = get_scheduler()
