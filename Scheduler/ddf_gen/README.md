Deep Drilling Field Generators
==============================

This directory contains an auxiliary files to generate pre-scheduled Deep Drilling Fields sequences for the Feature Based Scheduler configurations.

The ``.py`` files here can be executed in an rubin-scheduler environment to produce the required \`\`.npz\`\` files.
Each file will have a hash in their name that is created from the source files used to calculate them.
When creating a new scheduler configuration (or updating an existing one), you can execute the script to retrieve the hash.
This hash should then be added to the scheduler configuration, which should ensure the configuration will fail if the file does not exists or if it doesn't match the hash.
This hash is also written in the file itself so just changing the file name will not work.

Changes to this file result in changes to the resulting array of pre-optimized, pre-scheduled DDF observations, which are then used by the DDF ScriptedSurvey in the main survey configuration.
These changes are tracked by the hash of this file -- if the hash changes, recomputation of the prescheduled array is triggered.

