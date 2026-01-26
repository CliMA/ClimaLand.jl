# Benchmark Report for */home/vchuravy/src/HashArrayMappedTries*

## Job Properties
* Time of benchmark: 26 Aug 2023 - 16:12
* Package commit: dirty
* Julia commit: 661654
* Julia command flags: None
* Environment variables: None

## Results
Below is a table of this job's results, obtained by running the benchmarks.
The values listed in the `ID` column have the structure `[parent_group, child_group, ..., key]`, and can be used to
index into the BaseBenchmarks suite to retrieve the corresponding benchmarks.
The percentages accompanying time and memory values in the below table are noise tolerances. The "true"
time/memory value for a given benchmark is expected to fall within this percentage of the reported value.
An empty cell means that the value was zero.

| ID                                                   | time            | GC time   | memory          | allocations |
|------------------------------------------------------|----------------:|----------:|----------------:|------------:|
| `["Base.Dict", "creation (Persistent), size=0"]`     | 255.549 ns (5%) |           |  544 bytes (1%) |           4 |
| `["Base.Dict", "creation (Persistent), size=1"]`     | 375.539 ns (5%) |           |   1.06 KiB (1%) |           8 |
| `["Base.Dict", "creation (Persistent), size=1024"]`  |   3.237 ms (5%) |           |  37.36 MiB (1%) |        5005 |
| `["Base.Dict", "creation (Persistent), size=128"]`   |  63.660 μs (5%) |           | 450.48 KiB (1%) |         522 |
| `["Base.Dict", "creation (Persistent), size=16"]`    |   2.322 μs (5%) |           |  14.27 KiB (1%) |          71 |
| `["Base.Dict", "creation (Persistent), size=16384"]` | 973.768 ms (5%) | 88.025 ms |   7.94 GiB (1%) |      110830 |
| `["Base.Dict", "creation (Persistent), size=2"]`     | 480.256 ns (5%) |           |   1.59 KiB (1%) |          12 |
| `["Base.Dict", "creation (Persistent), size=2048"]`  |   8.968 ms (5%) |           | 105.71 MiB (1%) |       11149 |
| `["Base.Dict", "creation (Persistent), size=256"]`   | 246.730 μs (5%) |           |   2.10 MiB (1%) |        1037 |
| `["Base.Dict", "creation (Persistent), size=32"]`    |   4.483 μs (5%) |           |  35.52 KiB (1%) |         135 |
| `["Base.Dict", "creation (Persistent), size=4"]`     | 704.267 ns (5%) |           |   2.66 KiB (1%) |          20 |
| `["Base.Dict", "creation (Persistent), size=4096"]`  |  43.067 ms (5%) |  1.887 ms | 514.53 MiB (1%) |       24808 |
| `["Base.Dict", "creation (Persistent), size=512"]`   | 664.730 μs (5%) |           |   6.44 MiB (1%) |        2062 |
| `["Base.Dict", "creation (Persistent), size=64"]`    |  20.129 μs (5%) |           | 152.48 KiB (1%) |         266 |
| `["Base.Dict", "creation (Persistent), size=8"]`     |   1.106 μs (5%) |           |   4.78 KiB (1%) |          36 |
| `["Base.Dict", "creation (Persistent), size=8192"]`  | 132.564 ms (5%) |  6.889 ms |   1.57 GiB (1%) |       53480 |
| `["Base.Dict", "creation, size=0"]`                  | 249.415 ns (5%) |           |  544 bytes (1%) |           4 |
| `["Base.Dict", "creation, size=1"]`                  | 284.456 ns (5%) |           |  544 bytes (1%) |           4 |
| `["Base.Dict", "creation, size=1024"]`               |  17.830 μs (5%) |           |  91.97 KiB (1%) |          19 |
| `["Base.Dict", "creation, size=128"]`                |   2.428 μs (5%) |           |   6.36 KiB (1%) |          10 |
| `["Base.Dict", "creation, size=16"]`                 | 611.899 ns (5%) |           |   1.78 KiB (1%) |           7 |
| `["Base.Dict", "creation, size=16384"]`              | 364.750 μs (5%) |           |   1.42 MiB (1%) |          31 |
| `["Base.Dict", "creation, size=2"]`                  | 293.381 ns (5%) |           |  544 bytes (1%) |           4 |
| `["Base.Dict", "creation, size=2048"]`               |  31.470 μs (5%) |           |  91.97 KiB (1%) |          19 |
| `["Base.Dict", "creation, size=256"]`                |   5.704 μs (5%) |           |  23.67 KiB (1%) |          13 |
| `["Base.Dict", "creation, size=32"]`                 | 726.554 ns (5%) |           |   1.78 KiB (1%) |           7 |
| `["Base.Dict", "creation, size=4"]`                  | 302.078 ns (5%) |           |  544 bytes (1%) |           4 |
| `["Base.Dict", "creation, size=4096"]`               |  75.530 μs (5%) |           | 364.17 KiB (1%) |          25 |
| `["Base.Dict", "creation, size=512"]`                |   8.437 μs (5%) |           |  23.69 KiB (1%) |          14 |
| `["Base.Dict", "creation, size=64"]`                 |   1.835 μs (5%) |           |   6.36 KiB (1%) |          10 |
| `["Base.Dict", "creation, size=8"]`                  | 331.920 ns (5%) |           |  544 bytes (1%) |           4 |
| `["Base.Dict", "creation, size=8192"]`               | 139.200 μs (5%) |           | 364.17 KiB (1%) |          25 |
| `["Base.Dict", "delete, size=1"]`                    |  94.168 ns (5%) |           |  544 bytes (1%) |           4 |
| `["Base.Dict", "delete, size=1024"]`                 |   4.836 μs (5%) |           |  68.36 KiB (1%) |           6 |
| `["Base.Dict", "delete, size=128"]`                  | 651.962 ns (5%) |           |   4.66 KiB (1%) |           4 |
| `["Base.Dict", "delete, size=16"]`                   | 113.181 ns (5%) |           |   1.33 KiB (1%) |           4 |
| `["Base.Dict", "delete, size=16384"]`                |  74.580 μs (5%) |           |   1.06 MiB (1%) |           7 |
| `["Base.Dict", "delete, size=2"]`                    | 101.113 ns (5%) |           |  544 bytes (1%) |           4 |
| `["Base.Dict", "delete, size=2048"]`                 |   4.961 μs (5%) |           |  68.36 KiB (1%) |           6 |
| `["Base.Dict", "delete, size=256"]`                  |   1.446 μs (5%) |           |  17.39 KiB (1%) |           4 |
| `["Base.Dict", "delete, size=32"]`                   | 132.825 ns (5%) |           |   1.33 KiB (1%) |           4 |
| `["Base.Dict", "delete, size=4"]`                    |  99.864 ns (5%) |           |  544 bytes (1%) |           4 |
| `["Base.Dict", "delete, size=4096"]`                 |  18.209 μs (5%) |           | 272.28 KiB (1%) |           7 |
| `["Base.Dict", "delete, size=512"]`                  |   1.499 μs (5%) |           |  17.39 KiB (1%) |           4 |
| `["Base.Dict", "delete, size=64"]`                   | 644.902 ns (5%) |           |   4.66 KiB (1%) |           4 |
| `["Base.Dict", "delete, size=8"]`                    |  99.380 ns (5%) |           |  544 bytes (1%) |           4 |
| `["Base.Dict", "delete, size=8192"]`                 |  18.710 μs (5%) |           | 272.28 KiB (1%) |           7 |
| `["Base.Dict", "getindex, size=1"]`                  |   6.450 ns (5%) |           |                 |             |
| `["Base.Dict", "getindex, size=1024"]`               |   6.450 ns (5%) |           |                 |             |
| `["Base.Dict", "getindex, size=128"]`                |   7.317 ns (5%) |           |                 |             |
| `["Base.Dict", "getindex, size=16"]`                 |   6.460 ns (5%) |           |                 |             |
| `["Base.Dict", "getindex, size=16384"]`              |   6.910 ns (5%) |           |                 |             |
| `["Base.Dict", "getindex, size=2"]`                  |   6.450 ns (5%) |           |                 |             |
| `["Base.Dict", "getindex, size=2048"]`               |   6.900 ns (5%) |           |                 |             |
| `["Base.Dict", "getindex, size=256"]`                |   6.460 ns (5%) |           |                 |             |
| `["Base.Dict", "getindex, size=32"]`                 |   7.688 ns (5%) |           |                 |             |
| `["Base.Dict", "getindex, size=4"]`                  |   6.449 ns (5%) |           |                 |             |
| `["Base.Dict", "getindex, size=4096"]`               |   6.460 ns (5%) |           |                 |             |
| `["Base.Dict", "getindex, size=512"]`                |   6.920 ns (5%) |           |                 |             |
| `["Base.Dict", "getindex, size=64"]`                 |   6.460 ns (5%) |           |                 |             |
| `["Base.Dict", "getindex, size=8"]`                  |   6.460 ns (5%) |           |                 |             |
| `["Base.Dict", "getindex, size=8192"]`               |   7.257 ns (5%) |           |                 |             |
| `["Base.Dict", "insert, size=0"]`                    |  96.143 ns (5%) |           |  544 bytes (1%) |           4 |
| `["Base.Dict", "insert, size=1"]`                    | 102.561 ns (5%) |           |  544 bytes (1%) |           4 |
| `["Base.Dict", "insert, size=1024"]`                 |   4.547 μs (5%) |           |  68.36 KiB (1%) |           6 |
| `["Base.Dict", "insert, size=128"]`                  | 655.882 ns (5%) |           |   4.66 KiB (1%) |           4 |
| `["Base.Dict", "insert, size=16"]`                   | 118.819 ns (5%) |           |   1.33 KiB (1%) |           4 |
| `["Base.Dict", "insert, size=16384"]`                |  73.600 μs (5%) |           |   1.06 MiB (1%) |           7 |
| `["Base.Dict", "insert, size=2"]`                    | 102.608 ns (5%) |           |  544 bytes (1%) |           4 |
| `["Base.Dict", "insert, size=2048"]`                 |   4.816 μs (5%) |           |  68.36 KiB (1%) |           6 |
| `["Base.Dict", "insert, size=256"]`                  |   1.475 μs (5%) |           |  17.39 KiB (1%) |           4 |
| `["Base.Dict", "insert, size=32"]`                   | 117.118 ns (5%) |           |   1.33 KiB (1%) |           4 |
| `["Base.Dict", "insert, size=4"]`                    | 103.010 ns (5%) |           |  544 bytes (1%) |           4 |
| `["Base.Dict", "insert, size=4096"]`                 |  16.320 μs (5%) |           | 272.28 KiB (1%) |           7 |
| `["Base.Dict", "insert, size=512"]`                  |   1.515 μs (5%) |           |  17.39 KiB (1%) |           4 |
| `["Base.Dict", "insert, size=64"]`                   | 653.087 ns (5%) |           |   4.66 KiB (1%) |           4 |
| `["Base.Dict", "insert, size=8"]`                    | 103.189 ns (5%) |           |  544 bytes (1%) |           4 |
| `["Base.Dict", "insert, size=8192"]`                 |  18.520 μs (5%) |           | 272.28 KiB (1%) |           7 |
| `["HAMT", "creation (Persistent), size=0"]`          | 195.119 ns (5%) |           |   80 bytes (1%) |           2 |
| `["HAMT", "creation (Persistent), size=1"]`          | 273.527 ns (5%) |           |  272 bytes (1%) |           6 |
| `["HAMT", "creation (Persistent), size=1024"]`       | 169.620 μs (5%) |           | 820.34 KiB (1%) |        6883 |
| `["HAMT", "creation (Persistent), size=128"]`        |  15.640 μs (5%) |           |  71.72 KiB (1%) |         730 |
| `["HAMT", "creation (Persistent), size=16"]`         |   1.576 μs (5%) |           |   5.77 KiB (1%) |          75 |
| `["HAMT", "creation (Persistent), size=16384"]`      |   3.714 ms (5%) |           |  16.18 MiB (1%) |      136157 |
| `["HAMT", "creation (Persistent), size=2"]`          | 340.601 ns (5%) |           |  480 bytes (1%) |          10 |
| `["HAMT", "creation (Persistent), size=2048"]`       | 380.921 μs (5%) |           |   1.75 MiB (1%) |       14694 |
| `["HAMT", "creation (Persistent), size=256"]`        |  32.240 μs (5%) |           | 150.27 KiB (1%) |        1545 |
| `["HAMT", "creation (Persistent), size=32"]`         |   3.956 μs (5%) |           |  15.91 KiB (1%) |         158 |
| `["HAMT", "creation (Persistent), size=4"]`          | 463.418 ns (5%) |           |  912 bytes (1%) |          18 |
| `["HAMT", "creation (Persistent), size=4096"]`       | 782.240 μs (5%) |           |   3.56 MiB (1%) |       31229 |
| `["HAMT", "creation (Persistent), size=512"]`        |  73.050 μs (5%) |           | 347.98 KiB (1%) |        3249 |
| `["HAMT", "creation (Persistent), size=64"]`         |   7.907 μs (5%) |           |  34.41 KiB (1%) |         340 |
| `["HAMT", "creation (Persistent), size=8"]`          | 743.448 ns (5%) |           |   1.83 KiB (1%) |          34 |
| `["HAMT", "creation (Persistent), size=8192"]`       |   1.627 ms (5%) |           |   7.33 MiB (1%) |       65424 |
| `["HAMT", "creation, size=0"]`                       | 186.677 ns (5%) |           |   80 bytes (1%) |           2 |
| `["HAMT", "creation, size=1"]`                       | 256.185 ns (5%) |           |  192 bytes (1%) |           4 |
| `["HAMT", "creation, size=1024"]`                    |  51.960 μs (5%) |           |  95.03 KiB (1%) |        2060 |
| `["HAMT", "creation, size=128"]`                     |   6.004 μs (5%) |           |  11.05 KiB (1%) |         258 |
| `["HAMT", "creation, size=16"]`                      | 844.918 ns (5%) |           |   1.61 KiB (1%) |          32 |
| `["HAMT", "creation, size=16384"]`                   |   1.040 ms (5%) |           |   1.45 MiB (1%) |       29653 |
| `["HAMT", "creation, size=2"]`                       | 269.744 ns (5%) |           |  224 bytes (1%) |           5 |
| `["HAMT", "creation, size=2048"]`                    | 113.780 μs (5%) |           | 185.78 KiB (1%) |        4212 |
| `["HAMT", "creation, size=256"]`                     |  10.850 μs (5%) |           |  21.92 KiB (1%) |         465 |
| `["HAMT", "creation, size=32"]`                      |   1.695 μs (5%) |           |   3.36 KiB (1%) |          72 |
| `["HAMT", "creation, size=4"]`                       | 318.025 ns (5%) |           |  288 bytes (1%) |           7 |
| `["HAMT", "creation, size=4096"]`                    | 221.360 μs (5%) |           | 338.72 KiB (1%) |        7807 |
| `["HAMT", "creation, size=512"]`                     |  22.600 μs (5%) |           |  46.30 KiB (1%) |         946 |
| `["HAMT", "creation, size=64"]`                      |   3.085 μs (5%) |           |   5.77 KiB (1%) |         131 |
| `["HAMT", "creation, size=8"]`                       | 395.990 ns (5%) |           |  416 bytes (1%) |          11 |
| `["HAMT", "creation, size=8192"]`                    | 455.040 μs (5%) |           | 688.94 KiB (1%) |       14410 |
| `["HAMT", "delete, size=1"]`                         |  35.479 ns (5%) |           |   96 bytes (1%) |           2 |
| `["HAMT", "delete, size=1024"]`                      | 108.576 ns (5%) |           |  672 bytes (1%) |           6 |
| `["HAMT", "delete, size=128"]`                       | 102.283 ns (5%) |           |  544 bytes (1%) |           6 |
| `["HAMT", "delete, size=16"]`                        |  47.065 ns (5%) |           |  192 bytes (1%) |           2 |
| `["HAMT", "delete, size=16384"]`                     | 150.095 ns (5%) |           |  944 bytes (1%) |           8 |
| `["HAMT", "delete, size=2"]`                         |  35.146 ns (5%) |           |   96 bytes (1%) |           2 |
| `["HAMT", "delete, size=2048"]`                      | 113.906 ns (5%) |           |  768 bytes (1%) |           6 |
| `["HAMT", "delete, size=256"]`                       |  76.726 ns (5%) |           |  480 bytes (1%) |           4 |
| `["HAMT", "delete, size=32"]`                        |  69.162 ns (5%) |           |  352 bytes (1%) |           4 |
| `["HAMT", "delete, size=4"]`                         |  41.848 ns (5%) |           |  112 bytes (1%) |           2 |
| `["HAMT", "delete, size=4096"]`                      | 119.226 ns (5%) |           |  816 bytes (1%) |           6 |
| `["HAMT", "delete, size=512"]`                       |  78.642 ns (5%) |           |  528 bytes (1%) |           4 |
| `["HAMT", "delete, size=64"]`                        |  76.012 ns (5%) |           |  416 bytes (1%) |           4 |
| `["HAMT", "delete, size=8"]`                         |  46.559 ns (5%) |           |  144 bytes (1%) |           2 |
| `["HAMT", "delete, size=8192"]`                      | 119.310 ns (5%) |           |  848 bytes (1%) |           6 |
| `["HAMT", "getindex, size=1"]`                       |   5.220 ns (5%) |           |                 |             |
| `["HAMT", "getindex, size=1024"]`                    |   8.669 ns (5%) |           |                 |             |
| `["HAMT", "getindex, size=128"]`                     |   8.669 ns (5%) |           |                 |             |
| `["HAMT", "getindex, size=16"]`                      |   5.209 ns (5%) |           |                 |             |
| `["HAMT", "getindex, size=16384"]`                   |  11.572 ns (5%) |           |                 |             |
| `["HAMT", "getindex, size=2"]`                       |   5.210 ns (5%) |           |                 |             |
| `["HAMT", "getindex, size=2048"]`                    |   8.669 ns (5%) |           |                 |             |
| `["HAMT", "getindex, size=256"]`                     |   6.710 ns (5%) |           |                 |             |
| `["HAMT", "getindex, size=32"]`                      |   6.760 ns (5%) |           |                 |             |
| `["HAMT", "getindex, size=4"]`                       |   5.210 ns (5%) |           |                 |             |
| `["HAMT", "getindex, size=4096"]`                    |   8.669 ns (5%) |           |                 |             |
| `["HAMT", "getindex, size=512"]`                     |   6.860 ns (5%) |           |                 |             |
| `["HAMT", "getindex, size=64"]`                      |   6.780 ns (5%) |           |                 |             |
| `["HAMT", "getindex, size=8"]`                       |   5.209 ns (5%) |           |                 |             |
| `["HAMT", "getindex, size=8192"]`                    |   8.678 ns (5%) |           |                 |             |
| `["HAMT", "insert, size=0"]`                         |  65.594 ns (5%) |           |  192 bytes (1%) |           4 |
| `["HAMT", "insert, size=1"]`                         |  65.015 ns (5%) |           |  208 bytes (1%) |           4 |
| `["HAMT", "insert, size=1024"]`                      | 108.324 ns (5%) |           |   1.31 KiB (1%) |           6 |
| `["HAMT", "insert, size=128"]`                       | 108.866 ns (5%) |           |  576 bytes (1%) |           6 |
| `["HAMT", "insert, size=16"]`                        |  69.619 ns (5%) |           |  624 bytes (1%) |           4 |
| `["HAMT", "insert, size=16384"]`                     | 158.961 ns (5%) |           |   1.22 KiB (1%) |           8 |
| `["HAMT", "insert, size=2"]`                         |  63.245 ns (5%) |           |  208 bytes (1%) |           4 |
| `["HAMT", "insert, size=2048"]`                      | 197.957 ns (5%) |           |   1.45 KiB (1%) |           6 |
| `["HAMT", "insert, size=256"]`                       | 105.794 ns (5%) |           |  592 bytes (1%) |           6 |
| `["HAMT", "insert, size=32"]`                        | 101.133 ns (5%) |           |  464 bytes (1%) |           6 |
| `["HAMT", "insert, size=4"]`                         |  63.908 ns (5%) |           |  224 bytes (1%) |           4 |
| `["HAMT", "insert, size=4096"]`                      | 145.598 ns (5%) |           |  912 bytes (1%) |           8 |
| `["HAMT", "insert, size=512"]`                       | 153.153 ns (5%) |           |  672 bytes (1%) |           8 |
| `["HAMT", "insert, size=64"]`                        | 102.006 ns (5%) |           |  528 bytes (1%) |           6 |
| `["HAMT", "insert, size=8"]`                         |  69.670 ns (5%) |           |  512 bytes (1%) |           4 |
| `["HAMT", "insert, size=8192"]`                      | 166.904 ns (5%) |           |   1.20 KiB (1%) |           8 |


## Benchmark Group List
Here's a list of all the benchmark groups executed by this job:

- `["Base.Dict"]`
- `["HAMT"]`

## Julia versioninfo
```
Julia Version 1.10.0-beta1
Commit 6616549950e (2023-07-25 17:43 UTC)
Platform Info:
  OS: Linux (x86_64-linux-gnu)
      "Arch Linux"
  uname: Linux 6.3.2-arch1-1 #1 SMP PREEMPT_DYNAMIC Thu, 11 May 2023 16:40:42 +0000 x86_64 unknown
  CPU: AMD Ryzen 7 3700X 8-Core Processor: 
                 speed         user         nice          sys         idle          irq
       #1-16  2200 MHz    1214558 s       2006 s      98149 s   12788419 s      20501 s
  Memory: 125.69889831542969 GB (83676.40234375 MB free)
  Uptime: 367493.26 sec
  Load Avg:  1.12  1.03  0.83
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-15.0.7 (ORCJIT, znver2)
  Threads: 1 on 16 virtual cores
```