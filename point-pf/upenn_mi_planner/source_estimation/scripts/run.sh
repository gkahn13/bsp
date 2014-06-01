#!/bin/bash
set -o nounset
set -o errexit
TODAY=$(date +"%Y-%m-%d--")
for fldr in expt1 expt2 expt3 expt4 ; do
    cd $fldr
    for bag in $(ls *.bag); do 
        echo $bag
        TARGET=$(echo $bag | awk -F'_' '{print $1}')
        echo $(pwd)/$bag
        roslaunch source_estimation bag.launch \
            bagfile:=$(pwd)/$bag target:=${TARGET}
    done
    cd ..

    dest_fldr="/tmp/new/$fldr"
    mkdir -p $dest_fldr
    mv $fldr/*${TODAY}* $dest_fldr

    cd $dest_fldr
    for bag in $(ls *.bag); do
        TARGET=$(echo $bag | awk -F'_' '{print $1}')
        rosrun source_estimation bagread $bag ${TARGET}
    done
    cd /tmp/min
done
