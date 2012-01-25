#!/bin/bash

nice ~/rosetta/rosetta_source/bin/match.static.linuxiccrelease @general_matching.flags @scaf.flags @subs.flags -linmem_ig 10 -in:file::s $1
