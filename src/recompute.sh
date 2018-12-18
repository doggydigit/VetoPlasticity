#!/bin/bash

for line in $(seq 1 24)
do
    eval 'nohup /home/tsai/anaconda2/envs/p2/bin/python2.7 recompute.py $line >> log/log_line_s$1_l$line.out 2>&1 &'
    echo "Launched job $job"
    sleep 2s
done
