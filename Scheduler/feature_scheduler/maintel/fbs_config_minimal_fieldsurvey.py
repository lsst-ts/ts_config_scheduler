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

from pathlib import Path

from lsst.ts.fbs.utils.maintel.lsst_surveys import safety_masks
from lsst.ts.fbs.utils.maintel.make_fieldsurvey_scheduler import (
    MakeFieldSurveyScheduler,
    get_sv_targets,
)
from rubin_scheduler.scheduler import basis_functions, detailers


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

    # Get target coordinates - edit YAML file for updates
    # YAML file should be in the same directory as this .py config
    # ts_config_ocs/Scheduler/feature_scheduler/maintel/
    # fieldsurvey_centers.yaml
    target_dir = Path(__file__).parent
    target_file = Path.joinpath(target_dir, "fieldsurvey_centers.yaml")
    if not Path.exists(target_file):
        raise ValueError(f"Expected target yaml file does not exist at {target_file}")
    targets = get_sv_targets(target_file)

    make_scheduler = MakeFieldSurveyScheduler(
        targets=targets, nside=nside, ntiers=1, band_to_filter=band_to_filter
    )

    # Defaults

    config_basis_functions = safety_masks(
        nside=nside,
        min_alt=30.0,
        max_alt=83.0,
        shadow_minutes=10.0,
        wind_speed_maximum=40,
        apply_time_limited_shadow=False,
        time_to_sunrise=3.0,
        min_az_sunrise=144.0,
        max_az_sunrise=240.0,
    )
    config_basis_functions.append(
        basis_functions.SlewtimeBasisFunction(bandname=None, nside=nside)
    )

    observation_reason = "field_survey_science"
    science_program = "BLOCK-365"  # json BLOCK to be used

    nvisits = {"u": 30, "g": 30, "r": 30, "i": 30, "z": 30, "y": 30}
    sequence = ["y", "z", "i", "r"]
    # exposure time in seconds
    exptimes = {"u": 38.0, "g": 30.0, "r": 30.0, "i": 30.0, "z": 30.0, "y": 30.0}
    # 1 --> single 30 second exposure
    nexps = {"u": 1, "g": 1, "r": 1, "i": 1, "z": 1, "y": 1}
    field_survey_kwargs = {
        "nvisits": nvisits,
        "sequence": sequence,
        "exptimes": exptimes,
        "nexps": nexps,
    }

    # DENSELY DITHERED STAR FIELDS

    field_survey_kwargs["sequence"] = sequence

    # Densely dithered star fields
    # Spectrophotometric calibrators
    # Default for other fields
    # 2x raft scale dithers
    config_detailers = [
        detailers.DitherDetailer(max_dither=1.4, per_night=False),
        # Note: 1 center * 3 deg per visit * 30 visits = 90 deg
        detailers.CameraSmallRotPerObservationListDetailer(
            max_rot=45.0,
            min_rot=-45.0,
            per_visit_rot=3.0,
            # max_rot=15.0,
            # min_rot=-15.0,
            # per_visit_rot=1.0,
        ),
    ]
    tier = 0

    target_names = ["Rubin_SV_225_-40"]
    make_scheduler.add_field_surveys(
        tier,
        observation_reason,
        science_program,
        target_names,
        basis_functions=config_basis_functions,
        detailers=config_detailers,
        **field_survey_kwargs,
    )

    return make_scheduler.get_scheduler()


if __name__ == "config":
    nside, scheduler = get_scheduler()
