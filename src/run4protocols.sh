#!/bin/sh

python ziegler.py wTET >> log_wt.txt
python ziegler.py sTET >> log_st.txt
python ziegler.py wLFS >> log_wl.txt
python ziegler.py sLFS >> log_sl.txt
