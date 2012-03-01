cat ../../PLATO/run_script |sed "s/???/$1/" >run_script
cat ../../PLATO/run_talos |sed "s/???/$1/" >run_talos
cp ../../PLATO/runCSRjob.com .
cp ../$1.talos .
