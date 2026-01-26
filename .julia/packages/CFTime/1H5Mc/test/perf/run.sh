#!/bin/bash

# the correct answer should be
# (100*365+24)*24*60*60 = 3155673600

date
uname -a
cat /etc/issue
cat /proc/cpuinfo | grep 'model name' | head -n 1
echo "---------------------------"

julia benchmark-julia-CFTime.jl
echo "---------------------------"

julia benchmark-julia-Dates.jl
echo "---------------------------"

python3 benchmark-python-cftime.py
echo "---------------------------"

python3 benchmark-python-numpy.py
echo "---------------------------"

julia summary.jl
