#!/bin/bash

k="_r"

for job in $(seq 1 26)
do
    eval 'nohup /home/tsai/anaconda2/envs/p2/bin/python2.7 montesearch.py 0 >> log/log_s$1_j$job$k$2.out 2>&1 &'
    echo "Launched job $job"
    sleep 3s
done
