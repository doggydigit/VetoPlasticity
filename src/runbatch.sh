#!/bin/bash

k="_r"

for job in $(seq 1 16)
do
    eval 'nohup /home/tsai/anaconda2/envs/p2/bin/python2.7 montesearch.py >> log_s$1_j$job$k$2.out 2>&1 &'
    echo "Launched job $job"
    sleep 2s
done
