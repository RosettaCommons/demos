#!/bin/bash

#minimal demonstration of RosettaDNA RosettaScripts functionalities. As of 2011-03, these are identicial to the integration tests provided with the Rosetta source package.

# set your paths to the Rosetta executable and rosetta database
rosettaDNA=./rosettaDNA.linuxiccrelease
database_path=$HOME/minirosetta_database

$rosettaDNA @design.flags -database $database_path 2>&1 > design.log
$rosettaDNA @resfile.flags -database $database_path 2>&1 > resfile.log
$rosettaDNA @predspec.flags -database $database_path 2>&1 > predspec.log
