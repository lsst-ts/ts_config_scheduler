This directory contains an auxiliary file to define pre-scheduled sequences for FBS configurations.

The lsst_ddf_gen.py file contains the sequences for the fbs_config_lsst_survey.py configuration.
Changes to this file result in changes to the resulting array of pre-optimized, pre-scheduled DDF observations,
which are then used by the DDF ScriptedSurvey in the main survey configuration. These changes are tracked by
the hash of this file -- if the hash changes, recomputation of the prescheduled array is triggered.
