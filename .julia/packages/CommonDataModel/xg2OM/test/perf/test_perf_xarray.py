# Benchmark to be run on Linux as root
# dropping file cache is OS specific and requires root priviledges

import timeit
import xarray as xr
import numpy
import os
import sys


def mean_no_cache(ds):
    with open("/proc/sys/vm/drop_caches","w") as f:
        f.write("3")
    vm = ds["data"].groupby("time.month").mean().to_numpy();

def std_no_cache(ds):
    with open("/proc/sys/vm/drop_caches","w") as f:
        f.write("3")

    vm = ds["data"].groupby("time.month").std(ddof=1).to_numpy();

fname = os.path.expanduser("~/sample_perf2.nc")
ds = xr.open_dataset(fname)

print("python: ",sys.version)
print("xarray: ",xr.__version__)
print("numpy: ",numpy.__version__)


if __name__ == "__main__":
    print("runtime")

    for test_fun in [mean_no_cache, std_no_cache]:
        t = timeit.repeat(lambda: test_fun(ds),
                          setup="""from __main__ import ds""",
                          number=1,
                          repeat=30,
                          )

        print("  minimum time of ",test_fun,": ",min(t))


    month = ds["time.month"].to_numpy()

    print("accuracy")

    data_f64 = ds["data"].data.astype(dtype="f8")

    mean_ref = numpy.stack(
        [data_f64[(month == mm).nonzero()[0],:,:].mean(axis=0) for mm in range(1,13)],axis=0)

    std_ref = numpy.stack(
        [data_f64[(month == mm).nonzero()[0],:,:].std(axis=0,ddof=1) for mm in range(1,13)],axis=0)


    vm = ds["data"].groupby("time.month").mean().to_numpy();

    print("  accuracy of mean",
          numpy.sqrt(numpy.mean((mean_ref -  vm)**2)))

    vs = ds["data"].groupby("time.month").std(ddof=1).to_numpy()


    print("  accuracy of std",
          numpy.sqrt(numpy.mean((std_ref -  vs)**2)))
