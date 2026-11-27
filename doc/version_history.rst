.. _Version_History:

===============
Version History
===============

.. towncrier release notes start

v0.5.4 (2026-11-27)
===================

Bug Fixes
---------

- In ``Scheduler/ddf_gen/gen_ddf_presched_observations.py``, explicitly passed exposure times and number of snaps when generating ddf pre-scheduled observations. This fixed an issue where 15s exposures with 2 snaps were generated for the DDF fields. (`SSW-3267 <https://rubinobs.atlassian.net//browse/SSW-3267>`_)


v0.5.3 (2026-11-19)
===================

Other Changes and Additions
---------------------------

- Updated maintel observatory model rotator configuration to use follow_sky=false. (`SSW-3218 <https://rubinobs.atlassian.net//browse/SSW-3218>`_)
- Updated ocs observatory model rotator configuration to use follow_sky=false. (`SSW-3218 <https://rubinobs.atlassian.net//browse/SSW-3218>`_)
- Updated auxtel observatory model rotator configuration to have the right limits and to use follow_sky=false. (`SSW-3218 <https://rubinobs.atlassian.net//browse/SSW-3218>`_)


v0.5.2 (2026-11-17)
===================

Bug Fixes
---------

- Patch detailer order for pairs_u<g,r> surveys. (`LSST-84 <https://rubinobs.atlassian.net//browse/LSST-84>`_)


Other Changes and Additions
---------------------------

- Added new BLOCK-T770 observing block for MainTel horseshoe stray light tests, using LSSTCam in the r_57 filter with a sequence of az/el offsets around the bright star Alnair. (`BLOCK-T770 <https://rubinobs.atlassian.net//browse/BLOCK-T770>`_)
- Reduced the MTScheduler ToO client ``delta_time`` from 691200 s (8 days) to 600 s (10 minutes) in the v10 base configuration. (`BLOCK-T770.1 <https://rubinobs.atlassian.net//browse/BLOCK-T770.1>`_)
- Move EDFS ultradeep ddf season to year 2, instead of year 1 of the survey.
  Update ddf hash accordingly in fbs_config_lsst_survey configurations. (`LSST-76 <https://rubinobs.atlassian.net//browse/LSST-76>`_)
- Removing Electrometer:103 from CalSys scripts while it is not in service. Also added FiberSpectrographs back into BLOCK-T399. (`RSO-958 <https://rubinobs.atlassian.net//browse/RSO-958>`_)
- In Block-T601, default of use_electrometer = False. (`RSO-958.1 <https://rubinobs.atlassian.net//browse/RSO-958.1>`_)
- Reduced predicted schedule window to 0.5 hour. (`SSW-3218 <https://rubinobs.atlassian.net//browse/SSW-3218>`_)


v0.5.1 (2026-09-23)
===================

Other Changes and Additions
---------------------------

- Added github workflow to enforce rebasing, linting and news fragment creation. (`OSW-2607 <https://rubinobs.atlassian.net//browse/OSW-2607>`_)
- Added fiber spectrographs and some extra steps to CBP mask changes to weekly test case. (`RSO-855 <https://rubinobs.atlassian.net//browse/RSO-855>`_)
- Update FBS configurations to add queue_manager where appropriate. (`SP-3402 <https://rubinobs.atlassian.net//browse/SP-3402>`_)


v0.5.0 (2026-08-31)
===================

Documentation
-------------

- Add support for towncrier. (`OSW-2607 <https://rubinobs.atlassian.net//browse/OSW-2607>`_)


Other Changes and Additions
---------------------------

- Removed Scheduler v7 configuration files. (`OSW-2607 <https://rubinobs.atlassian.net//browse/OSW-2607>`_)
- Removed Scheduler v8 configuration files. (`OSW-2607 <https://rubinobs.atlassian.net//browse/OSW-2607>`_)
- Created Scheduler v10 configuration files. (`OSW-2607 <https://rubinobs.atlassian.net//browse/OSW-2607>`_)
- Clean up v10 configuration files, lifting several parameters that were previously unnecessarily defined in the site and specific configurations to default status. (`OSW-2607 <https://rubinobs.atlassian.net//browse/OSW-2607>`_)
- Updated expected hash for scheduled observation arrays. (`OSW-2607 <https://rubinobs.atlassian.net//browse/OSW-2607>`_)
- Update license headers. (`OSW-2607 <https://rubinobs.atlassian.net//browse/OSW-2607>`_)
- Add license check. (`OSW-2607 <https://rubinobs.atlassian.net//browse/OSW-2607>`_)
