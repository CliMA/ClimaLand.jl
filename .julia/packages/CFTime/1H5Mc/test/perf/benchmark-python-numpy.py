# Install dependencies via the shell command:
#
# pip install numpy

import numpy as np
import timeit
import sys


def compute(n):
    t0 = np.datetime64("1900-01-01") + np.arange(0,n).astype("timedelta64[s]")
    t1 = np.datetime64("2000-01-01") + np.arange(0,n).astype("timedelta64[s]")
    diff = t1 - np.flip(t0)

    return diff.astype("int64").mean()

if __name__ == "__main__":
    print("python: ",sys.version)
    print("numpy: ",np.__version__)

    n = 1_000_000
#    n = 100_000
    mean_total_seconds = compute(n)
    print("mean_total_seconds: ", mean_total_seconds)

    setup = "from __main__ import compute"
    benchtime = timeit.repeat(lambda: compute(n), setup=setup,number = 1, repeat = 100)

    print("min time: ",min(benchtime))

    with open("python-numpy.txt","w") as f:
        for bt in benchtime:
            print(bt,file=f)
