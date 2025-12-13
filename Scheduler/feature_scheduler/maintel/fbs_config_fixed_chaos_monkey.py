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

import lsst.ts.fbs.utils.maintel.lsst_surveys as lsst_surveys
from rubin_scheduler.scheduler.basis_functions import (
    VisitGap,
)
from rubin_scheduler.scheduler.detailers import AltAz2RaDecDetailer, ZeroRotDetailer
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

    science_program = "BLOCK-T644"

    # FYI: There is a bug in the calculation of
    # the sky angle for the FieldAltAzSurvey which
    # which forced me to hard code the sky angle in
    # the observation block. If you change the az
    # you will have to update the sky angle in the
    # BLOCK-T644 block.
    regular_images_survey_alt = 60
    regular_images_survey_az = 270

    regular_images_survey_basis_functions = lsst_surveys.safety_masks(
        apply_time_limited_shadow=False, shadow_minutes=10
    )
    regular_images_survey_sequence = ["r"]
    regular_images_survey_nvisits = dict(r=1)
    regular_images_survey_exptimes = dict(r=30.0)
    regular_images_survey_nexps = dict(r=1)
    regular_images_survey_target_name = (
        f"Field{regular_images_survey_az}_{regular_images_survey_alt}"
    )
    regular_images_survey_science_program = science_program
    regular_images_survey_detailers = [
        AltAz2RaDecDetailer(),
        ZeroRotDetailer(),
    ]

    regular_images_survey = FieldAltAzSurvey(
        basis_functions=regular_images_survey_basis_functions,
        alt=regular_images_survey_alt,
        az=regular_images_survey_az,
        sequence=regular_images_survey_sequence,
        nvisits=regular_images_survey_nvisits,
        exptimes=regular_images_survey_exptimes,
        nexps=regular_images_survey_nexps,
        ignore_obs=None,
        survey_name=regular_images_survey_target_name,
        target_name=regular_images_survey_target_name,
        science_program=regular_images_survey_science_program,
        observation_reason="FixedChaosMonkeyTest",
        scheduler_note=regular_images_survey_target_name,
        nside=nside,
        flush_pad=30.0,
        detailers=regular_images_survey_detailers,
    )

    perturbation_surveys = []

    for i in range(6):
        perturbation_survey_name = (
            f"AdditionalField{regular_images_survey_az}_"
            f"{regular_images_survey_alt}_{i+1}"
        )
        perturbation_basis_functions = lsst_surveys.safety_masks(
            apply_time_limited_shadow=False, shadow_minutes=10
        )
        perturbation_basis_functions.append(
            VisitGap(
                note=perturbation_survey_name,
                gap_min=(i + 1) * 30.0,
                band_names=["r"],
            )
        )
        perturbation_block = f"BLOCK_T644_{i+1}"
        perturbation_survey_detailers = [
            AltAz2RaDecDetailer(),
            ZeroRotDetailer(),
        ]
        survey = FieldAltAzSurvey(
            basis_functions=perturbation_basis_functions,
            alt=regular_images_survey_alt,
            az=regular_images_survey_az,
            sequence=regular_images_survey_sequence,
            nvisits=regular_images_survey_nvisits,
            exptimes=regular_images_survey_exptimes,
            nexps=regular_images_survey_nexps,
            ignore_obs=None,
            survey_name=perturbation_survey_name,
            target_name=regular_images_survey_target_name,
            science_program=perturbation_block,
            observation_reason=f"FixedChaosMonkeyTest_perturbation_{i+1}",
            scheduler_note=perturbation_survey_name,
            nside=nside,
            flush_pad=30.0,
            detailers=perturbation_survey_detailers,
        )
        perturbation_surveys.append(survey)

    survey_lists = [
        perturbation_surveys,
        [regular_images_survey],
    ]

    return nside, CoreScheduler(
        survey_lists,
        nside=nside,
        band_to_filter=band_to_filter,
    )


if __name__ == "config":
    nside, scheduler = get_scheduler()
