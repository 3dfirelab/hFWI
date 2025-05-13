#!/bin/bash
#
export srcDir=/home/paugam/Src/hFWI/src/
export dataDir=/mnt/data3/SILEX/AROME/
export logDir=/mnt/data3/SILEX/AROME/FWI/log
mkdir -p $logDir

source ~/.myKey.sh
source ~/Venv/arome/bin/activate
python $srcDir/fwi.py $dataDir >& $logDir/fwi.log
rm $dataDir/FWI/timeToCompute.txt

touch $dataDir/FWI/PLOT_latest/run_ffmpeg_on_moor
#python /home/paugam/bin/render_anim/render_anim.py -i /mnt/data3/SILEX/AROME/FWI/PLOT_latest/png/ -ext png -o /mnt/data3/SILEX/AROME/FWI/PLOT_latest/fwi_latest.avi -fr 1
#ffmpeg -i /mnt/data3/SILEX/AROME/FWI/PLOT_latest/fwi_latest.avi /mnt/data3/SILEX/AROME/FWI/PLOT_latest/fwi_latest.mp4

