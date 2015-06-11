#!/usr/local/bin/bash
# Script file to execute all post-processing scripts at once
# Emily Koo

echo "Post-processing started."
echo "Refresh and wait for .png files to be generated."

# Surface contact map frequency plot (ads only)
nohup PlotSurfaceContactMap.py &

# Secondary structure histogram plot
nohup PlotSecStruct.sh &

# Contact map frequency plot
nohup PlotContactMap.sh &

