#!/bin/bash

for line in $(seq 64 84)
do
    eval 'nohup /home/tsai/anaconda2/envs/p2/bin/python2.7 recompute.py $line >> log/log_line_s$1_l$line.out 2>&1 &'
    echo "Launched job $line"
    sleep 3s
done
