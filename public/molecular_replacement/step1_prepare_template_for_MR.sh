#!/bin/sh

mkdir templates
cd templates
$ROSETTA3_SRC/src/apps/public/electron_density/prepare_template_for_MR.pl ../inputs/1crb.hhr
