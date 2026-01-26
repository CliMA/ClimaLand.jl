# Install dependencies via the shell commands:
#
# pip install cftime numpy

import numpy as np
import cftime
import timeit
import sys
import cftime
import datetime

def compute(n):
    t0 = cftime.num2date(np.arange(n),"milliseconds since 1900-01-1",calendar="proleptic_gregorian")
    t1 = cftime.num2date(np.arange(n),"milliseconds since 2000-01-1",calendar="proleptic_gregorian")

    diff = t1 - np.flip(t0)

    mean_total_seconds = np.vectorize(lambda d: d.total_seconds(),otypes = [int])(diff).mean()
    return mean_total_seconds


if __name__ == "__main__":
    print("python: ",sys.version)
    print("cftime: ",cftime.__version__)
    n = 1_000_000
#    n = 100_000
    mean_total_seconds = compute(n)
    print("mean_total_seconds: ", mean_total_seconds)

    setup = "from __main__ import compute"
    benchtime = timeit.repeat(lambda: compute(n), setup=setup,number = 1, repeat = 100)

    print("min time: ",min(benchtime))

    with open("python-cftime.txt","w") as f:
        for bt in benchtime:
            print(bt,file=f)
