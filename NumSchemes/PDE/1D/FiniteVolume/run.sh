#!/bin/bash
schemes=("LW" "WB" "SOU" "Third")

for scheme in "${schemes[@]}"
do
    ./main.x -o data/"$scheme"s.h5 -n 3 "$scheme" "$scheme" "$scheme" no vanleer superbee
done

cd pproc || return

for scheme in "${schemes[@]}"
do
    echo "Postprocessing: " data/"$scheme"s.h5
    python h5pproc.py -d ../data/"$scheme"s.h5
done

cd ../ || return