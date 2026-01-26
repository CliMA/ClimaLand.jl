# Benchmark between xarray and groupby of CommonDataModel


Test case is 3D array of size 360x180x10959 (one degree resolution global dataset representing 30 years of daily data).
Normally distributed random data  (mean = 100, variance = 1) in single precision floats (`Float32`).
We compute the mean and standard deviation (std) of data grouped by month.

Accuracy is assessed by comparison with built-in functions (Statistics.jl or numpy) in double precision (using `Float64`).
Note that julia’s `mean`/`std` give exactly the same results as numpy’s equivalent.

Using 1 CPU core, xarray’s default implementations (i. e. no dask…)
30 trials, minimum time is reported here
Ubuntu 22.04, Julia 1.11, python 3.10.12, xarray 2024.12


Creation of the data file:

```bash
julia test_perf_init.jl
```

Get root priviledges (to drop file cache)

```
sudo -s
export HOME=/home/abarth
cd ~/.julia/dev/CommonDataModel/test/perf
```

## Laptop with a i5-1135G7 CPU and NVMe SSD WDC WDS100T2B0C


### CommonDataModel

```bash
~/.juliaup/bin/julia test_perf_cdm.jl
```

Output:

```
runtime of mean
  2.133 s (1686528 allocations: 2.71 GiB)
runtime of std
  2.574 s (1686525 allocations: 2.72 GiB)
accuracy
sqrt(mean((gm - mean_ref) .^ 2)) = 4.643795867042341e-5
sqrt(mean((gs - std_ref) .^ 2)) = 9.251281717683748e-7
```


### xarray

```bash
python3 test_perf_xarray.py
```

Output:

```
python:  3.10.12 (main, Sep 11 2024, 15:47:36) [GCC 11.4.0]
xarray:  2024.10.0
numpy:  1.26.1
runtime
  minimum time of  <function mean_no_cache at 0x741456563d90> :  4.260775363999983
  minimum time of  <function std_no_cache at 0x74143df69240> :  5.453749345000006
accuracy
  accuracy of mean 4.64379586704234e-05
  accuracy of std 2.1211715725139616e-07
```


# Workstation with i7-7700 CPU and SATA SSD (WD Green 120G)

```
~/.juliaup/bin/julia test_perf_init.jl
~/.juliaup/bin/julia test_perf_cdm.jl
python3 test_perf_xarray.py
```

Output:

```
runtime
  7.177 s (1686528 allocations: 2.71 GiB)
  8.090 s (1686525 allocations: 2.72 GiB)
accuracy
sqrt(mean((gm - mean_ref) .^ 2)) = 4.6300139982730906e-5
sqrt(mean((gs - std_ref) .^ 2)) = 9.268973317814482e-7
python:  3.10.12 (main, Jul 29 2024, 16:56:48) [GCC 11.4.0]
xarray:  2024.10.0
numpy:  1.26.2
runtime
  minimum time of  <function mean_no_cache at 0x7f54cb64bd90> :  8.740452307043597
  minimum time of  <function std_no_cache at 0x7f54b31a0a60> :  10.462690721964464
accuracy
  accuracy of mean 4.6300139982730906e-05
  accuracy of std 2.12226305970758e-07
```
