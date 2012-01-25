#!/bin/bash

nice ~/rosetta/rosetta_source/bin/enzyme_design.static.linuxiccrelease @enzdes.flags -correct -linmem_ig 10 -in:file::s $1
