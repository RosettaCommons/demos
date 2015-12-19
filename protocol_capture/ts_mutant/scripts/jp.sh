#!/bin/bash
# wrapper to run java file converter
dn=$(dirname $0)
java -Dmodeller.root=$MODELLER_ROOT -cp ${dn}/javaparse app.TsParseCSV $*
