# Benchmark Report for */home/vchuravy/src/ScopedVariables*

## Job Properties
* Time of benchmark: 26 Aug 2023 - 11:10
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

| ID                                          | time            | GC time | memory         | allocations |
|---------------------------------------------|----------------:|--------:|---------------:|------------:|
| `["BASIC", "scoped with assignment & ref"]` | 202.131 ns (5%) |         | 368 bytes (1%) |          10 |
| `["BASIC", "scoped with assignment"]`       | 187.893 ns (5%) |         | 368 bytes (1%) |          10 |
| `["BASIC", "scoped"]`                       |   1.249 ns (5%) |         |                |             |
| `["BASIC", "unscoped"]`                     |   4.469 ns (5%) |         |                |             |
| `["DEPTH", "access 2^i, depth=1"]`          |   9.660 ns (5%) |         |                |             |
| `["DEPTH", "access 2^i, depth=1024"]`       |  12.900 ns (5%) |         |                |             |
| `["DEPTH", "access 2^i, depth=128"]`        |   9.663 ns (5%) |         |                |             |
| `["DEPTH", "access 2^i, depth=16"]`         |   9.620 ns (5%) |         |                |             |
| `["DEPTH", "access 2^i, depth=2"]`          |   9.660 ns (5%) |         |                |             |
| `["DEPTH", "access 2^i, depth=2048"]`       |  13.333 ns (5%) |         |                |             |
| `["DEPTH", "access 2^i, depth=256"]`        |   9.721 ns (5%) |         |                |             |
| `["DEPTH", "access 2^i, depth=32"]`         |   9.659 ns (5%) |         |                |             |
| `["DEPTH", "access 2^i, depth=4"]`          |   9.660 ns (5%) |         |                |             |
| `["DEPTH", "access 2^i, depth=4096"]`       |  14.286 ns (5%) |         |                |             |
| `["DEPTH", "access 2^i, depth=512"]`        |   9.792 ns (5%) |         |                |             |
| `["DEPTH", "access 2^i, depth=64"]`         |   9.624 ns (5%) |         |                |             |
| `["DEPTH", "access 2^i, depth=8"]`          |   9.660 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=1"]`              |  10.640 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=10"]`             |   9.660 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=11"]`             |   9.659 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=12"]`             |   9.659 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=13"]`             |   9.590 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=14"]`             |   9.620 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=15"]`             |   9.630 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=16"]`             |   9.610 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=17"]`             |  11.872 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=18"]`             |   9.630 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=19"]`             |   9.659 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=2"]`              |   9.660 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=20"]`             |   9.659 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=21"]`             |   9.659 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=22"]`             |   9.609 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=23"]`             |   9.599 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=24"]`             |   9.658 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=25"]`             |   9.658 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=26"]`             |   9.609 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=27"]`             |   9.658 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=28"]`             |   9.659 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=29"]`             |   9.659 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=3"]`              |   9.660 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=30"]`             |   9.629 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=31"]`             |   9.659 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=32"]`             |   9.568 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=33"]`             |   9.657 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=34"]`             |   9.658 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=35"]`             |   9.657 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=36"]`             |   9.638 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=37"]`             |   9.658 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=38"]`             |   9.638 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=39"]`             |   9.638 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=4"]`              |   9.620 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=40"]`             |   9.636 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=41"]`             |   9.657 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=42"]`             |   9.556 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=43"]`             |   9.617 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=44"]`             |   9.617 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=45"]`             |   9.656 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=46"]`             |   9.656 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=47"]`             |   9.627 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=48"]`             |   9.637 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=49"]`             |   9.616 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=5"]`              |   9.660 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=50"]`             |   9.606 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=51"]`             |   9.625 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=52"]`             |   9.655 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=53"]`             |   9.655 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=54"]`             |   9.656 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=55"]`             |   9.605 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=56"]`             |   9.615 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=57"]`             |   9.655 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=58"]`             |   9.664 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=59"]`             |   9.655 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=6"]`              |   9.660 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=60"]`             |   9.634 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=61"]`             |   9.553 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=62"]`             |   9.664 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=63"]`             |   9.654 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=64"]`             |   9.664 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=7"]`              |   9.610 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=8"]`              |   9.660 ns (5%) |         |                |             |
| `["DEPTH", "access, depth=9"]`              |   9.660 ns (5%) |         |                |             |
| `["DEPTH", "emtpyf, depth=1"]`              | 137.242 ns (5%) |         | 240 bytes (1%) |           8 |
| `["DEPTH", "emtpyf, depth=10"]`             |   1.505 μs (5%) |         |  2.34 KiB (1%) |          80 |
| `["DEPTH", "emtpyf, depth=11"]`             |   1.643 μs (5%) |         |  2.58 KiB (1%) |          88 |
| `["DEPTH", "emtpyf, depth=12"]`             |   1.811 μs (5%) |         |  2.81 KiB (1%) |          96 |
| `["DEPTH", "emtpyf, depth=13"]`             |   1.978 μs (5%) |         |  3.05 KiB (1%) |         104 |
| `["DEPTH", "emtpyf, depth=14"]`             |   2.144 μs (5%) |         |  3.28 KiB (1%) |         112 |
| `["DEPTH", "emtpyf, depth=15"]`             |   2.260 μs (5%) |         |  3.52 KiB (1%) |         120 |
| `["DEPTH", "emtpyf, depth=16"]`             |   2.511 μs (5%) |         |  3.75 KiB (1%) |         128 |
| `["DEPTH", "emtpyf, depth=17"]`             |   2.624 μs (5%) |         |  3.98 KiB (1%) |         136 |
| `["DEPTH", "emtpyf, depth=18"]`             |   2.801 μs (5%) |         |  4.22 KiB (1%) |         144 |
| `["DEPTH", "emtpyf, depth=19"]`             |   2.956 μs (5%) |         |  4.45 KiB (1%) |         152 |
| `["DEPTH", "emtpyf, depth=2"]`              | 267.079 ns (5%) |         | 480 bytes (1%) |          16 |
| `["DEPTH", "emtpyf, depth=20"]`             |   3.159 μs (5%) |         |  4.69 KiB (1%) |         160 |
| `["DEPTH", "emtpyf, depth=21"]`             |   3.210 μs (5%) |         |  4.92 KiB (1%) |         168 |
| `["DEPTH", "emtpyf, depth=22"]`             |   3.468 μs (5%) |         |  5.16 KiB (1%) |         176 |
| `["DEPTH", "emtpyf, depth=23"]`             |   3.623 μs (5%) |         |  5.39 KiB (1%) |         184 |
| `["DEPTH", "emtpyf, depth=24"]`             |   3.629 μs (5%) |         |  5.62 KiB (1%) |         192 |
| `["DEPTH", "emtpyf, depth=25"]`             |   3.937 μs (5%) |         |  5.86 KiB (1%) |         200 |
| `["DEPTH", "emtpyf, depth=26"]`             |   4.100 μs (5%) |         |  6.09 KiB (1%) |         208 |
| `["DEPTH", "emtpyf, depth=27"]`             |   4.136 μs (5%) |         |  6.33 KiB (1%) |         216 |
| `["DEPTH", "emtpyf, depth=28"]`             |   4.359 μs (5%) |         |  6.56 KiB (1%) |         224 |
| `["DEPTH", "emtpyf, depth=29"]`             |   4.549 μs (5%) |         |  6.80 KiB (1%) |         232 |
| `["DEPTH", "emtpyf, depth=3"]`              | 394.304 ns (5%) |         | 720 bytes (1%) |          24 |
| `["DEPTH", "emtpyf, depth=30"]`             |   4.765 μs (5%) |         |  7.03 KiB (1%) |         240 |
| `["DEPTH", "emtpyf, depth=31"]`             |   4.964 μs (5%) |         |  7.27 KiB (1%) |         248 |
| `["DEPTH", "emtpyf, depth=32"]`             |   5.099 μs (5%) |         |  7.50 KiB (1%) |         256 |
| `["DEPTH", "emtpyf, depth=33"]`             |   5.324 μs (5%) |         |  7.73 KiB (1%) |         264 |
| `["DEPTH", "emtpyf, depth=34"]`             |   5.315 μs (5%) |         |  7.97 KiB (1%) |         272 |
| `["DEPTH", "emtpyf, depth=35"]`             |   5.604 μs (5%) |         |  8.20 KiB (1%) |         280 |
| `["DEPTH", "emtpyf, depth=36"]`             |   5.511 μs (5%) |         |  8.44 KiB (1%) |         288 |
| `["DEPTH", "emtpyf, depth=37"]`             |   6.017 μs (5%) |         |  8.67 KiB (1%) |         296 |
| `["DEPTH", "emtpyf, depth=38"]`             |   6.097 μs (5%) |         |  8.91 KiB (1%) |         304 |
| `["DEPTH", "emtpyf, depth=39"]`             |   6.271 μs (5%) |         |  9.14 KiB (1%) |         312 |
| `["DEPTH", "emtpyf, depth=4"]`              | 556.000 ns (5%) |         | 960 bytes (1%) |          32 |
| `["DEPTH", "emtpyf, depth=40"]`             |   6.439 μs (5%) |         |  9.38 KiB (1%) |         320 |
| `["DEPTH", "emtpyf, depth=41"]`             |   6.719 μs (5%) |         |  9.61 KiB (1%) |         328 |
| `["DEPTH", "emtpyf, depth=42"]`             |   6.699 μs (5%) |         |  9.84 KiB (1%) |         336 |
| `["DEPTH", "emtpyf, depth=43"]`             |   6.914 μs (5%) |         | 10.08 KiB (1%) |         344 |
| `["DEPTH", "emtpyf, depth=44"]`             |   6.984 μs (5%) |         | 10.31 KiB (1%) |         352 |
| `["DEPTH", "emtpyf, depth=45"]`             |   7.345 μs (5%) |         | 10.55 KiB (1%) |         360 |
| `["DEPTH", "emtpyf, depth=46"]`             |   7.252 μs (5%) |         | 10.78 KiB (1%) |         368 |
| `["DEPTH", "emtpyf, depth=47"]`             |   7.582 μs (5%) |         | 11.02 KiB (1%) |         376 |
| `["DEPTH", "emtpyf, depth=48"]`             |   7.747 μs (5%) |         | 11.25 KiB (1%) |         384 |
| `["DEPTH", "emtpyf, depth=49"]`             |   8.000 μs (5%) |         | 11.48 KiB (1%) |         392 |
| `["DEPTH", "emtpyf, depth=5"]`              | 709.495 ns (5%) |         |  1.17 KiB (1%) |          40 |
| `["DEPTH", "emtpyf, depth=50"]`             |   8.145 μs (5%) |         | 11.72 KiB (1%) |         400 |
| `["DEPTH", "emtpyf, depth=51"]`             |   8.347 μs (5%) |         | 11.95 KiB (1%) |         408 |
| `["DEPTH", "emtpyf, depth=52"]`             |   8.432 μs (5%) |         | 12.19 KiB (1%) |         416 |
| `["DEPTH", "emtpyf, depth=53"]`             |   8.614 μs (5%) |         | 12.42 KiB (1%) |         424 |
| `["DEPTH", "emtpyf, depth=54"]`             |   8.670 μs (5%) |         | 12.66 KiB (1%) |         432 |
| `["DEPTH", "emtpyf, depth=55"]`             |   8.670 μs (5%) |         | 12.89 KiB (1%) |         440 |
| `["DEPTH", "emtpyf, depth=56"]`             |   9.260 μs (5%) |         | 13.12 KiB (1%) |         448 |
| `["DEPTH", "emtpyf, depth=57"]`             |   9.212 μs (5%) |         | 13.36 KiB (1%) |         456 |
| `["DEPTH", "emtpyf, depth=58"]`             |   9.556 μs (5%) |         | 13.59 KiB (1%) |         464 |
| `["DEPTH", "emtpyf, depth=59"]`             |   9.564 μs (5%) |         | 13.83 KiB (1%) |         472 |
| `["DEPTH", "emtpyf, depth=6"]`              | 864.339 ns (5%) |         |  1.41 KiB (1%) |          48 |
| `["DEPTH", "emtpyf, depth=60"]`             |   9.216 μs (5%) |         | 14.06 KiB (1%) |         480 |
| `["DEPTH", "emtpyf, depth=61"]`             |   9.572 μs (5%) |         | 14.30 KiB (1%) |         488 |
| `["DEPTH", "emtpyf, depth=62"]`             |   9.718 μs (5%) |         | 14.53 KiB (1%) |         496 |
| `["DEPTH", "emtpyf, depth=63"]`             |  10.198 μs (5%) |         | 14.77 KiB (1%) |         504 |
| `["DEPTH", "emtpyf, depth=64"]`             |  10.328 μs (5%) |         | 15.00 KiB (1%) |         512 |
| `["DEPTH", "emtpyf, depth=7"]`              |   1.037 μs (5%) |         |  1.64 KiB (1%) |          56 |
| `["DEPTH", "emtpyf, depth=8"]`              |   1.194 μs (5%) |         |  1.88 KiB (1%) |          64 |
| `["DEPTH", "emtpyf, depth=9"]`              |   1.340 μs (5%) |         |  2.11 KiB (1%) |          72 |
| `["DEPTH", "scope + access, depth=1"]`      | 150.230 ns (5%) |         | 240 bytes (1%) |           8 |
| `["DEPTH", "scope + access, depth=10"]`     |   1.602 μs (5%) |         |  2.34 KiB (1%) |          80 |
| `["DEPTH", "scope + access, depth=11"]`     |   1.840 μs (5%) |         |  2.58 KiB (1%) |          88 |
| `["DEPTH", "scope + access, depth=12"]`     |   2.152 μs (5%) |         |  2.81 KiB (1%) |          96 |
| `["DEPTH", "scope + access, depth=13"]`     |   2.198 μs (5%) |         |  3.05 KiB (1%) |         104 |
| `["DEPTH", "scope + access, depth=14"]`     |   2.361 μs (5%) |         |  3.28 KiB (1%) |         112 |
| `["DEPTH", "scope + access, depth=15"]`     |   2.392 μs (5%) |         |  3.52 KiB (1%) |         120 |
| `["DEPTH", "scope + access, depth=16"]`     |   2.604 μs (5%) |         |  3.75 KiB (1%) |         128 |
| `["DEPTH", "scope + access, depth=17"]`     |   2.807 μs (5%) |         |  3.98 KiB (1%) |         136 |
| `["DEPTH", "scope + access, depth=18"]`     |   3.061 μs (5%) |         |  4.22 KiB (1%) |         144 |
| `["DEPTH", "scope + access, depth=19"]`     |   3.036 μs (5%) |         |  4.45 KiB (1%) |         152 |
| `["DEPTH", "scope + access, depth=2"]`      | 288.401 ns (5%) |         | 480 bytes (1%) |          16 |
| `["DEPTH", "scope + access, depth=20"]`     |   3.290 μs (5%) |         |  4.69 KiB (1%) |         160 |
| `["DEPTH", "scope + access, depth=21"]`     |   3.470 μs (5%) |         |  4.92 KiB (1%) |         168 |
| `["DEPTH", "scope + access, depth=22"]`     |   3.816 μs (5%) |         |  5.16 KiB (1%) |         176 |
| `["DEPTH", "scope + access, depth=23"]`     |   3.704 μs (5%) |         |  5.39 KiB (1%) |         184 |
| `["DEPTH", "scope + access, depth=24"]`     |   3.956 μs (5%) |         |  5.62 KiB (1%) |         192 |
| `["DEPTH", "scope + access, depth=25"]`     |   3.984 μs (5%) |         |  5.86 KiB (1%) |         200 |
| `["DEPTH", "scope + access, depth=26"]`     |   4.427 μs (5%) |         |  6.09 KiB (1%) |         208 |
| `["DEPTH", "scope + access, depth=27"]`     |   4.633 μs (5%) |         |  6.33 KiB (1%) |         216 |
| `["DEPTH", "scope + access, depth=28"]`     |   4.711 μs (5%) |         |  6.56 KiB (1%) |         224 |
| `["DEPTH", "scope + access, depth=29"]`     |   4.721 μs (5%) |         |  6.80 KiB (1%) |         232 |
| `["DEPTH", "scope + access, depth=3"]`      | 439.271 ns (5%) |         | 720 bytes (1%) |          24 |
| `["DEPTH", "scope + access, depth=30"]`     |   5.253 μs (5%) |         |  7.03 KiB (1%) |         240 |
| `["DEPTH", "scope + access, depth=31"]`     |   5.531 μs (5%) |         |  7.27 KiB (1%) |         248 |
| `["DEPTH", "scope + access, depth=32"]`     |   5.277 μs (5%) |         |  7.50 KiB (1%) |         256 |
| `["DEPTH", "scope + access, depth=33"]`     |   5.617 μs (5%) |         |  7.73 KiB (1%) |         264 |
| `["DEPTH", "scope + access, depth=34"]`     |   5.880 μs (5%) |         |  7.97 KiB (1%) |         272 |
| `["DEPTH", "scope + access, depth=35"]`     |   5.960 μs (5%) |         |  8.20 KiB (1%) |         280 |
| `["DEPTH", "scope + access, depth=36"]`     |   6.178 μs (5%) |         |  8.44 KiB (1%) |         288 |
| `["DEPTH", "scope + access, depth=37"]`     |   6.184 μs (5%) |         |  8.67 KiB (1%) |         296 |
| `["DEPTH", "scope + access, depth=38"]`     |   6.574 μs (5%) |         |  8.91 KiB (1%) |         304 |
| `["DEPTH", "scope + access, depth=39"]`     |   6.543 μs (5%) |         |  9.14 KiB (1%) |         312 |
| `["DEPTH", "scope + access, depth=4"]`      | 589.799 ns (5%) |         | 960 bytes (1%) |          32 |
| `["DEPTH", "scope + access, depth=40"]`     |   6.906 μs (5%) |         |  9.38 KiB (1%) |         320 |
| `["DEPTH", "scope + access, depth=41"]`     |   7.049 μs (5%) |         |  9.61 KiB (1%) |         328 |
| `["DEPTH", "scope + access, depth=42"]`     |   7.067 μs (5%) |         |  9.84 KiB (1%) |         336 |
| `["DEPTH", "scope + access, depth=43"]`     |   7.290 μs (5%) |         | 10.08 KiB (1%) |         344 |
| `["DEPTH", "scope + access, depth=44"]`     |   7.317 μs (5%) |         | 10.31 KiB (1%) |         352 |
| `["DEPTH", "scope + access, depth=45"]`     |   7.962 μs (5%) |         | 10.55 KiB (1%) |         360 |
| `["DEPTH", "scope + access, depth=46"]`     |   7.738 μs (5%) |         | 10.78 KiB (1%) |         368 |
| `["DEPTH", "scope + access, depth=47"]`     |   8.062 μs (5%) |         | 11.02 KiB (1%) |         376 |
| `["DEPTH", "scope + access, depth=48"]`     |   8.220 μs (5%) |         | 11.25 KiB (1%) |         384 |
| `["DEPTH", "scope + access, depth=49"]`     |   8.415 μs (5%) |         | 11.48 KiB (1%) |         392 |
| `["DEPTH", "scope + access, depth=5"]`      | 797.249 ns (5%) |         |  1.17 KiB (1%) |          40 |
| `["DEPTH", "scope + access, depth=50"]`     |   8.624 μs (5%) |         | 11.72 KiB (1%) |         400 |
| `["DEPTH", "scope + access, depth=51"]`     |   8.757 μs (5%) |         | 11.95 KiB (1%) |         408 |
| `["DEPTH", "scope + access, depth=52"]`     |   8.990 μs (5%) |         | 12.19 KiB (1%) |         416 |
| `["DEPTH", "scope + access, depth=53"]`     |   8.740 μs (5%) |         | 12.42 KiB (1%) |         424 |
| `["DEPTH", "scope + access, depth=54"]`     |   9.172 μs (5%) |         | 12.66 KiB (1%) |         432 |
| `["DEPTH", "scope + access, depth=55"]`     |   9.374 μs (5%) |         | 12.89 KiB (1%) |         440 |
| `["DEPTH", "scope + access, depth=56"]`     |   9.476 μs (5%) |         | 13.12 KiB (1%) |         448 |
| `["DEPTH", "scope + access, depth=57"]`     |   9.802 μs (5%) |         | 13.36 KiB (1%) |         456 |
| `["DEPTH", "scope + access, depth=58"]`     |   9.912 μs (5%) |         | 13.59 KiB (1%) |         464 |
| `["DEPTH", "scope + access, depth=59"]`     |   9.690 μs (5%) |         | 13.83 KiB (1%) |         472 |
| `["DEPTH", "scope + access, depth=6"]`      | 941.471 ns (5%) |         |  1.41 KiB (1%) |          48 |
| `["DEPTH", "scope + access, depth=60"]`     |  10.055 μs (5%) |         | 14.06 KiB (1%) |         480 |
| `["DEPTH", "scope + access, depth=61"]`     |  10.475 μs (5%) |         | 14.30 KiB (1%) |         488 |
| `["DEPTH", "scope + access, depth=62"]`     |  10.730 μs (5%) |         | 14.53 KiB (1%) |         496 |
| `["DEPTH", "scope + access, depth=63"]`     |  10.742 μs (5%) |         | 14.77 KiB (1%) |         504 |
| `["DEPTH", "scope + access, depth=64"]`     |  10.915 μs (5%) |         | 15.00 KiB (1%) |         512 |
| `["DEPTH", "scope + access, depth=7"]`      |   1.085 μs (5%) |         |  1.64 KiB (1%) |          56 |
| `["DEPTH", "scope + access, depth=8"]`      |   1.250 μs (5%) |         |  1.88 KiB (1%) |          64 |
| `["DEPTH", "scope + access, depth=9"]`      |   1.395 μs (5%) |         |  2.11 KiB (1%) |          72 |
| `["WIDTH", "access 2^i, width=1"]`          |   9.660 ns (5%) |         |                |             |
| `["WIDTH", "access 2^i, width=1024"]`       |  16.000 ns (5%) |         |                |             |
| `["WIDTH", "access 2^i, width=128"]`        |   9.663 ns (5%) |         |                |             |
| `["WIDTH", "access 2^i, width=16"]`         |   9.599 ns (5%) |         |                |             |
| `["WIDTH", "access 2^i, width=2"]`          |   9.660 ns (5%) |         |                |             |
| `["WIDTH", "access 2^i, width=2048"]`       |  17.778 ns (5%) |         |                |             |
| `["WIDTH", "access 2^i, width=256"]`        |   9.939 ns (5%) |         |                |             |
| `["WIDTH", "access 2^i, width=32"]`         |   9.648 ns (5%) |         |                |             |
| `["WIDTH", "access 2^i, width=4"]`          |   9.630 ns (5%) |         |                |             |
| `["WIDTH", "access 2^i, width=4096"]`       |  20.000 ns (5%) |         |                |             |
| `["WIDTH", "access 2^i, width=512"]`        |  11.771 ns (5%) |         |                |             |
| `["WIDTH", "access 2^i, width=64"]`         |   9.674 ns (5%) |         |                |             |
| `["WIDTH", "access 2^i, width=8"]`          |   9.660 ns (5%) |         |                |             |

## Benchmark Group List
Here's a list of all the benchmark groups executed by this job:

- `["BASIC"]`
- `["DEPTH"]`
- `["WIDTH"]`

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
       #1-16  4049 MHz    1182978 s       1870 s      89439 s    9926959 s      19717 s
  Memory: 125.69889831542969 GB (82936.19921875 MB free)
  Uptime: 349340.11 sec
  Load Avg:  1.06  0.69  0.3
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-15.0.7 (ORCJIT, znver2)
  Threads: 1 on 16 virtual cores
```