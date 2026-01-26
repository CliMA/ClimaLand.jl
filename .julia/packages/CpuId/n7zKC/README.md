# CPU, what can you do?

*CpuId* is a package for the Julia programming language that enables you to
query the availability of specific CPU features with low run-time cost
using the assembly instruction `cpuid`.

Test                        | Status
----------------------------|----------------------
Windows, Linux & Mac Build  | [![Build Status](https://travis-ci.org/m-j-w/CpuId.jl.svg?branch=master)](https://travis-ci.org/m-j-w/CpuId.jl)
Code Coverage               | [![codecov](https://codecov.io/gh/m-j-w/CpuId.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/m-j-w/CpuId.jl)

_Status: considered a beta version for the core functionality, ready for you to try out._

Works on Julia 1.0 and later, on Linux, Mac and Windows with Intel CPUs
and AMD CPUs.  Other processor types like ARM are _not_ supported.


## Motivation

Besides the obvious reason to gather information for diagnostics, the CPU
provides valuable information when aiming at increasing the efficiency of code.
Such use-cases could be to tailor the size of working sets of data according to
the available cache sizes, to detect when the code is executed in a virtual
machine (hypervisor), or to determine the size of the largest SIMD registers
available to choose the best algorithm for the current hardware.

This information is obtained by directly querying the CPU through the `cpuid`
assembly instruction which operates only using CPU registers, and provides
a portable way to adapt code to specific hardware.

Same information may of course be collected from various sources, from Julia
itself or from the operating system, e.g. on Linux from `/proc/cpuinfo`.  See
below for a few [alternatives](#alternatives).  However, the `cpuid` instruction
is portable in the sense that it doesn't rely on other external dependencies.

The full documentation is found in Intel's 4670 page combined [Architectures
Software Devleoper Manual](
http://www.intel.com/content/www/us/en/architecture-and-technology/64-ia-32-architectures-software-developer-manual-325462.html).
A more practical and concise overview of the `cpuid` instruction is found at
[sandpile.org](http://www.sandpile.org/x86/cpuid.htm).


## Installation and Usage

*CpuId* is a registered Julia package; use the package manager to install:

    Julia> Pkg.add("CpuId")

Or, if you're keen to get some intermediate updates, clone from GitHub
[master branch](https://github.com/m-j-w/CpuId.jl/tree/master):

    Julia> Pkg.clone("https://github.com/m-j-w/CpuId.jl")

See the diagnostic summary on your CPU by typing

```
julia> using CpuId
julia> cpuinfo()

   Cpu Property         Value
  ╾───────────────────╌───────────────────────────────────────────────────────────╼
   Brand                Intel(R) Xeon(R) CPU E3-1225 v5 @ 3.30GHz
   Vendor               :Intel
   Architecture         :Skylake
   Model                Family: 6, Model: 94, Stepping: 3, Type: 0
   Cores                4 physical cores, 4 logical cores (on executing CPU)
                        No Hyperthreading detected
   Clock Frequencies    3300 / 3700 MHz (base/max), 100 MHz bus
   Data Cache           Level 1:3 : (32, 256, 8192) kbytes
                        64 byte cache line size
   Address Size         48 bits virtual, 39 bits physical
   SIMD                 256 bit = 32 byte max. SIMD vector size
   Time Stamp Counter   TSC is accessible via `rdtsc`
                        TSC runs at constant rate (invariant from clock frequency)
   Perf. Monitoring     Performance Monitoring Counters (PMC) revision 4
                        Available hardware counters per logical core:
                        3 fixed-function counters of 48 bit width
                        8 general-purpose counters of 48 bit width
   Hypervisor           No
```

Or get a list of the feature flags of your CPU with
```
julia> cpufeaturetable()

   Cpu Feature   Description
  ╾────────────╌───────────────────────────────────────────────────────────────╼
   3DNowP        3D Now PREFETCH and PREFETCHW instructions
   ACPI          Thermal monitor and software controlled clock facilities (MSR)
   ADX           Intel ADX (Multi-Precision Add-Carry Instruction Extensions)
   AES           AES encryption instruction set
   AHF64         LAHF and SAHF in PM64
   APIC          APIC on-chip (Advanced Programmable Interrupt Controller)
   AVX           256bit Advanced Vector Extensions, AVX
   AVX2          SIMD 256bit Advanced Vector Extensions 2
   BMI1          Bit Manipulation Instruction Set 1
   BMI2          Bit Manipulation Instruction Set 2
   CLFLUSH       CLFLUSHOPT Instructions
   CLFSH         CLFLUSH instruction (SSE2)
   CMOV          Conditional move CMOV and FCMOV instructions
   CX16          CMPXCHG16B instruction
   CX8           CMPXCHG8 instruction (64bit compare and exchange)
   ...
```

## Features

This release covers a selection of fundamental and higher level functionality:

 - `cpuinfo()` generates the summary shown above (markdown string).
 - `cpubrand()`, `cpumodel()`, `cpuvendor()` allow the identification of the
     CPU.
 - `cpuarchitecture()` tries to infer the microarchitecture, currently only of
     Intel CPUs.
 - `cpucores()` and `cputhreads()` to determine the number of physical and
     logical cores on the currently executing CPU, which typically share L3
     caches and main memory bandwidth.  If the result of both functions is
     equal, then the CPU does *not* use of hyperthreading.
 - `address_size()` and `physical_address_size()` return the number of bits used
     in pointers.  Useful when stealing a few bits from a pointer.
 - `cachelinesize()` gives the size in bytes of one cache line, which is
     typically 64 bytes.
 - `cachesize()` returns a tuple with the sizes of the data caches in bytes.
 - `cacheinclusive()` returns a tuple indicating lower cache levels being
    included in the data cache sizes reported by `cachesize()`.
 - `cpu_base_frequency()`, `cpu_max_frequency()`, `cpu_bus_frequency()` give -
     if supported by the CPU, the base, maximum and bus clock frequencies.
     Use `has_cpu_frequencies()` to check whether this property is supported.
 - `hypervised()` returns true when the CPU indicates that a hypervisor is
     running the operating system, aka a virtual machine.  In that case,
     `hvvendor()` may be invoked to get the, well, hypervisor vendor, and
     `hvversion()` returns a dictionary of additional version tags.
     `hvinfo()` generates a markdown summary of same dictionary.
 - `simdbits()` and `simdbytes()` return the size of the largest SIMD register
     available on the executing CPU.
 - `perf_revision()` to query the revision number of hardware performance
     monitoring counters, along with `perf_fix_counters()`, `perf_gen_counters()`,
     `perf_fix_bits()`, `perf_gen_bits()` to determine the number and bit width
     of available fixed-function and general purpose counters per logical core.
 - `cpucycle()` and `cpucycle_id()` let you directly get the CPU's time stamp
     counter, which is increased for every CPU clock cycle. This is a method to
     perform low overhead micro-benchmarking; though, technically, this uses the
     `rdtsc` and `rdtscp` instructions rather than `cpuid`.
 - `cpufeature(::Symbol)` permits asking for the availability of a specific
     feature, and `cpufeaturetable()` gives a complete overview of all detected
     features with a brief explanation, as shown above.


## Background

The `cpuid` instruction is a generic way provided by the CPU vendor to obtain
basic hardware information.  It provides data in form of boolean bit fields,
integer fields and strings, all packed in the returned CPU registers EAX, EBX,
ECX and EDX. Which information is returned is determined by the so called leaf,
which is defined by setting the input register EAX to a specific 32 bit integer
value before executing the instruction.  The extent and kind of information
obtainable via this instruction has changed quite a lot over the past decade and
still evolves with every CPU generation.  Thus, not all information is available
on every CPU model, and certainly everything is vendor dependent.

This Julia package also provides the `cpucycle()` function which allows getting
the current time stamp counter (TSC) by emitting a `rdtsc` instruction.
Similarly to `cpuid`, it only requires CPU registers and is usable in user-land
code and facilitates an alternative approach to micro-benchmarking.


## Limitations

The behaviour on non-Intel CPUs is currently unknown; though technically a crash
of Julia could be expected, theoretically, a rather large list of CPUs support
the `cpuid` instruction. Tip: Just try it and report back.

There are plenty of different CPUs, and in particular the `cpuid` instruction
has numerous corner cases, which this package does not address, yet.

Moreover, the `cpuid` instruction can only provide information for the executing
physical CPU, called a package.  To obtain information on all packages, and all
physical and logical cores, the executing program must be pinned sequentially to
each and every core, and gather its properties. This is how `libuv`, `hwloc` or
the operating system obtain that kind information.  However, this would require
additional external or operating system dependent code which is not the scope of
this package.

#### Specific limitations and caveats

- The number of physical cores and logical cores reported by `CpuId` seems wrong!
    If you have multiple processors on your motherboard, then `CpuId` will
    always only give you information for the processor the current task is running
    on.  For example: You have 2 processors, each with 12 physical cores and 24
    logical cores (thus with hyperthreading).  While you have in total 48 logical
    cores on both processors, `CpuId` will only give you 24 logical and 12 physical
    cores from the one it is running on.
    Resolving this is outside the scope of this Julia module, since it requires
    additional other operating system dependent functions, pinning the current
    task to a specific CPU, or querying other BIOS related functions.

- Why aren't all infos available that are seen e.g. in `/proc/cpuinfo`?
    Many of those features, flags and properties reside in the so called machine
    specific registers (MSR), which are only accessible to privileged programs
    running in the so called *ring 0*, such as the Linux kernel itself. Thus,
    short answer: You're not allowed to.

- The results obtained by `CpuId` functions are inconsistent with my hardware!
    Try other programs whether they give the same information or differ. If they
    differ, then you found a bug.  See below for some
    [alternatives](#alternatives), in particular the Linux command line tool
    *cpuid*.

- When running a hypervisor (virtual machine) the presented information is wrong!
    Hypervisor vendors are free to provide the `cpuid` information
    by intercepting calls to that instruction.  Not all vendors comply, and some
    even permit the user to change what is reported.  Thus, expect some
    surprises when a hypervisor is detected.

- My hypervisor is not detected!
    Either you're not really running a hypervisor, e.g. *Bash on Windows* is
    _not_ a virtual machine, or there is a feature missing. Raise an issue on
    GitHub.

- `cpucycles()` invokes `rdtsc`; that is not `cpuid`!
    True. However, both are valuable when diagnosing performance issues
    and trying to perform micro benchmarks on specific hardware.


## Alternatives

On Linux, most of the information may be obtained by reading from the `/proc`
tree, in particular `/proc/cpuinfo`, which eventually also invokes the `cpuid`
instruction.  Type `man 4 cpuid` to get a brief description of this kernel
driver.

On many Linux distributions, there is also the command line tool [cpuid](
http://www.etallen.com/cpuid.html), which essentially does exactly the same.  On
Ubuntu, you would install it using `sudo apt install cpuid`, then use it to show
a summary by simply typing `cpuid`.

Then, of course, there are a few functions in Julia Base. These are
`Base.Sys.cpu_info()`, and `Base.Sys.cpu_summary()`, as well as the global
variable `Base.Sys.CPU_THREADS`.  These are mostly provided by wrapping *libuv*.
In particular `CPU_THREADS` is the reason for this module: This reports the
number of logical cores, but how many physical cores do you have that you would
want to run your code on?

The Julia package [Hwloc.jl](https://github.com/JuliaParallel/Hwloc.jl) provides
similar and more information primarily directed towards the topology of your
CPUs, viz.  number of CPU packages, physical & logical cores and associated
caches, along with a number of features to deal with thread affinity. However,
it also pulls in additional external binary dependencies in that it relies on
[hwloc](https://www.open-mpi.org/projects/hwloc/), which also implies quite some
computational overhead. Whether this is an issue in the first place depends much
on your use-case.


## Terms of usage

This Julia package *CpuId* is published as open source and licensed under the
[MIT "Expat" License](./LICENSE.md).


**Contributions welcome!**

You're welcome to report successful usage or any issues via GitHub, and to open
pull requests to extend the current functionality.

