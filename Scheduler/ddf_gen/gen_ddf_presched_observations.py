#!/usr/bin/env python3
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

__all__ = ("define_ddf_seq", "gen_ddf_presched_observations")

import copy
import hashlib
from pathlib import Path

import lsst.ts.fbs.utils.maintel.lsst_ddf_presched as ddf_presched
import numpy as np
import pandas as pd
from lsst.ts.fbs.utils.maintel.lsst_surveys import (
    EXPTIME,
    NEXP,
    SCIENCE_PROGRAM,
    U_EXPTIME,
    U_NEXP,
)
from rubin_scheduler.data import get_data_dir
from rubin_scheduler.utils import SURVEY_START_MJD


def define_ddf_seq() -> pd.DataFrame:
    """Define the sequences for each field"""

    short_squences = [
        {
            "u": 3,
            "season_length": 225,
            "flush_length": 2.0,
            "g_depth_limit": 23.5,
            "n_sequences": 33,
        },
        {
            "y": 2,
            "season_length": 225,
            "flush_length": 2.0,
            "g_depth_limit": 23.5,
            "n_sequences": 33,
        },
        {
            "g": 2,
            "i": 2,
            "season_length": 225,
            "flush_length": 0.5,
            "g_depth_limit": 22.8,
            "n_sequences": 56,
            "even_odd": "even",
        },
        {
            "r": 2,
            "z": 2,
            "season_length": 225,
            "flush_length": 0.5,
            "g_depth_limit": 22.8,
            "n_sequences": 56,
            "even_odd": "odd",
        },
    ]

    shallow_squences = [
        {
            "u": 3,
            "season_length": 225,
            "flush_length": 2.0,
            "g_depth_limit": 23.5,
            "n_sequences": 75,
        },
        {
            "y": 2,
            "season_length": 225,
            "flush_length": 2.0,
            "g_depth_limit": 23.5,
            "n_sequences": 75,
        },
        {
            "g": 2,
            "i": 2,
            "season_length": 225,
            "flush_length": 0.5,
            "g_depth_limit": 22.8,
            "n_sequences": 100,
            "even_odd": "even",
        },
        {
            "r": 2,
            "z": 2,
            "season_length": 225,
            "flush_length": 0.5,
            "g_depth_limit": 22.8,
            "n_sequences": 100,
            "even_odd": "odd",
        },
    ]

    deep_sequences = [
        {
            "u": 8,
            "season_length": 225,
            "flush_length": 2.0,
            "g_depth_limit": 23.5,
            "n_sequences": 75,
        },
        {
            "y": 20,
            "season_length": 225,
            "flush_length": 2.0,
            "g_depth_limit": 23.5,
            "n_sequences": 75,
        },
        {
            "g": 2,
            "i": 2,
            "season_length": 225,
            "flush_length": 0.5,
            "g_depth_limit": 22.8,
            "n_sequences": 200,
            "even_odd": "even",
        },
        {
            "r": 2,
            "z": 2,
            "season_length": 225,
            "flush_length": 0.5,
            "g_depth_limit": 22.8,
            "n_sequences": 200,
            "even_odd": "odd",
        },
        {
            "g": 4,
            "r": 18,
            "i": 55,
            "z": 52,
            "season_length": 180,
            "flush_length": 2.0,
            "g_depth_limit": 23.5,
            "n_sequences": 110,
        },
    ]

    euclid_deep_seq = [
        {
            "u": 30,
            "season_length": 225,
            "flush_length": 2.0,
            "g_depth_limit": 23.5,
            "n_sequences": 75,
        },
        {
            "y": 40,
            "season_length": 225,
            "flush_length": 2.0,
            "g_depth_limit": 23.5,
            "n_sequences": 75,
        },
        {
            "g": 4,
            "i": 4,
            "season_length": 225,
            "flush_length": 0.5,
            "g_depth_limit": 22.8,
            "n_sequences": 200,
            "even_odd": "even",
        },
        {
            "r": 4,
            "z": 4,
            "season_length": 225,
            "flush_length": 0.5,
            "g_depth_limit": 22.8,
            "n_sequences": 200,
            "even_odd": "odd",
        },
        {
            "g": 8,
            "r": 36,
            "i": 110,
            "z": 104,
            "season_length": 125,
            "flush_length": 2.0,
            "g_depth_limit": 23.5,
            "n_sequences": 75,
        },
    ]

    short_seasons = {
        "XMM_LSS": [0, 10],
        "ELAISS1": [0, 10],
        "ECDFS": [0, 10],
        "EDFS_a": [0, 10],
    }

    shallow_seasons = {
        "COSMOS": [0, 4, 5, 6, 7, 8, 9, 10],
        "XMM_LSS": [1, 2, 3, 5, 6, 7, 8, 9],
        "ELAISS1": [1, 2, 3, 4, 6, 7, 8, 9],
        "ECDFS": [1, 2, 3, 4, 5, 7, 8, 9],
        "EDFS_a": [2, 3, 4, 5, 6, 7, 8, 9],
    }

    deep_seasons = {
        "COSMOS": [1, 2, 3],
        "XMM_LSS": [4],
        "ELAISS1": [5],
        "ECDFS": [6],
        "EDFS_a": [1],
    }

    dataframes = []

    for ddf_name in short_seasons:
        for season in short_seasons[ddf_name]:
            dict_for_df = {
                "ddf_name": ddf_name,
                "season": season,
                "even_odd": "None",
            }
            for key in "ugrizy":
                dict_for_df[key] = 0

            for seq in short_squences:
                row = copy.copy(dict_for_df)
                for key in seq:
                    row[key] = seq[key]
                dataframes.append(pd.DataFrame.from_dict(row, orient="index").T)

    for ddf_name in shallow_seasons:
        for season in shallow_seasons[ddf_name]:
            dict_for_df = {
                "ddf_name": ddf_name,
                "season": season,
                "even_odd": "None",
            }
            for key in "ugrizy":
                dict_for_df[key] = 0

            for seq in shallow_squences:
                row = copy.copy(dict_for_df)
                for key in seq:
                    row[key] = seq[key]
                dataframes.append(pd.DataFrame.from_dict(row, orient="index").T)

    for ddf_name in deep_seasons:
        for season in deep_seasons[ddf_name]:
            dict_for_df = {
                "ddf_name": ddf_name,
                "season": season,
                "even_odd": "None",
            }
            for key in "ugrizy":
                dict_for_df[key] = 0
            if ddf_name == "EDFS_a":
                for seq in euclid_deep_seq:
                    row = copy.copy(dict_for_df)
                    for key in seq:
                        row[key] = seq[key]
                    dataframes.append(pd.DataFrame.from_dict(row, orient="index").T)
            else:
                for seq in deep_sequences:
                    row = copy.copy(dict_for_df)
                    for key in seq:
                        row[key] = seq[key]

                    dataframes.append(pd.DataFrame.from_dict(row, orient="index").T)

    result = pd.concat(dataframes)

    return result


def gen_ddf_presched_observations(
    expt: dict | None = None,
    nexp: dict | None = None,
    survey_start: float = SURVEY_START_MJD,
    survey_length: int = 10,
    survey_name: str = "deep drilling",
    science_program: str = SCIENCE_PROGRAM,
    save_filename: str = "ts_ddf_array.npz",
    save_path: str = None,
    additional_hash_files: list[str] = [],
) -> None:
    """Generate surveys for DDF observations.

    Parameters
    ----------
    expt : `dict`  { `str` : `float` } or None
        Exposure time for DDF visits.
        Default of None uses defaults of EXPTIME/U_EXPTIME.
    nexp : `dict` { `str` : `int` } or None
        Number of exposures per visit.
        Default of None uses defaults of NEXP/U_NEXP.
    survey_start : `float`
        Start MJD of the survey. Used for prescheduling DDF visits.
    survey_length : `float`
        Length of the survey. Used for prescheduling DDF visits.
        In years.
    science_program : `str`
        Name of the science program for the Survey.
    save_filename : `str`
        Filename of the saved ddf array.
    save_path : `str`
        Path to saved DDF file. If none, uses get_data_dir to look for it.
    """
    if expt is None:
        expt = {
            "u": U_EXPTIME,
            "g": EXPTIME,
            "r": EXPTIME,
            "i": EXPTIME,
            "z": EXPTIME,
            "y": EXPTIME,
        }
    if nexp is None:
        nexp = {"u": U_NEXP, "g": NEXP, "r": NEXP, "i": NEXP, "z": NEXP, "y": NEXP}

    if save_path is None:
        save_path = Path(get_data_dir(), "scheduler")
    # Potetial pre-computed obs_array:
    root, ext = save_filename.rsplit(".", maxsplit=1)

    # Hash of the files that define the DDF sequences, to identify
    # if the saved file comes from the same sources.
    hash_files = additional_hash_files + [__file__, ddf_presched.__file__]
    print(f'hash files: {", ".join(hash_files)}')
    hash_digest = ddf_presched.calculate_checksum(hash_files)
    hash_object = hashlib.sha256()
    hash_object.update(hash_digest)
    hex_digest = hash_object.hexdigest()[:7]
    pre_comp_file = Path(save_path, f"{root}_{hex_digest}.{ext}")
    if pre_comp_file.exists():
        print(f"Pre computed DDF file {pre_comp_file} already exists. Skipping...")
        return
    print(f"Generating DDF array: {hex_digest}")
    ddf_dataframe = define_ddf_seq()

    obs_array = ddf_presched.generate_ddf_scheduled_obs(
        ddf_dataframe,
        expt=expt,
        nsnaps=nexp,
        survey_start_mjd=survey_start,
        survey_length=survey_length,
        science_program=science_program,
    )

    print(f"Saving DDF array to {pre_comp_file}.")
    np.savez(
        pre_comp_file,
        obs_array=obs_array.view(np.ndarray),
        hash_digest=hash_digest,
        expt=expt,
        nexp=nexp,
        survey_start=survey_start,
        survey_length=survey_length,
        science_program=science_program,
    )
