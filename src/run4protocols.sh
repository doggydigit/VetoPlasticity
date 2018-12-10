#!/bin/sh

py='/home/tsai/anaconda2/envs/p2/bin/python2.7 ziegler.py '

eval $py 'wTET >> log_wt.txt'
eval $py 'sTET >> log_st.txt'
eval $py 'wLFS >> log_wl.txt'
eval $py 'sLFS >> log_sl.txt'
