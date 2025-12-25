# Parallel-Cahn-Hilliard-PDE-Solver
A parallel numerical solver for the Cahn–Hilliard partial differential equation. 

- Contributors: Anqi Chen, Zhang Wu, Raphael Rozenberg, Matias Molinolo De Ferrari

<i>Note: This project is migrated from code.harvard.edu - a Harvard-owned version of GitHub where all contributors originally worked on. </i>

## Directory Structure of Code and Application Files
```
Parallel-Cahn-Hilliard-PDE-Solver/
└── sequential_baseline/
    ── test_data   
        ├── frames
        ├── raw
        ├── animation.gif      
    ── Makefile
    ── gmon.out
    ── plot.py
    ── profiling.txt
    ── python_prototyping.ipynb
    ── sequential.cpp
|    
└── code/
    ├── measure_sequential_portion
       ├── CahnHilliard.cpp
       ├── CahnHilliard.h
       ├── Makefile
       ├── main.cpp
       ├── time_parallel.sh
    ├── performance_analysis
       ├── performance_analysis
          ├── N100
            ── Result_postprocess.ipynb
          ├── N350
            ── Result_postprocess.ipynb
       ├── CahnHilliard.cpp
       ├── CahnHilliard.h
       ├── Makefile
       ├── main.cpp
       ├── performance_multiple_both.sh
       ├── performance_multiple_ranks.sh
       ├── performance_multiple_threads.sh
    ├── CahnHilliard.cpp
    ├── CahnHilliard.h
    ├── Makefile
    ├── example_job.sh
    ├── main.cpp
    ├── plot.py
│ 
├── LICENSE
├── .gitignore
├── README.md  

```

## How to Compile and Run Our Code
- sp1: clone our git repo on cluster and navigate into the 'code' folder:
```bash
$ module load git/2.17.0-fasrc01
$ git clone https://github.com/zhangwu01/Parallel-Cahn-Hilliard-PDE-Solver.git
$ cd ./Parallel-Cahn-Hilliard-PDE-Solver/code
```

- sp2: compile the code with the command:
```bash
$ make
```
- sp3: run the parallized code with the command:
```bash
$ ./main
```
- sp4: reproduce the performance results with the command:
```bash
$ sbatch performance_multiple_both.sh
$ sbatch performance_multiple_ranks.sh
$ sbatch performance_multiple_threads.sh
```

## Information of the System and the Enviroment

### Cluster Architecture

- Architecture:          x86\_64
- CPU op-mode(s):        32-bit, 64-bit
- Byte Order:            Little Endian
- CPU(s):                32
- On-line CPU(s) list:   0-31
- Thread(s) per core:    1
- Core(s) per socket:    16
- Socket(s):             2
- NUMA node(s):          2
- Vendor ID:             GenuineIntel
- CPU family:            6
- Model:                 79
- Model name:            Intel(R) Xeon(R) CPU E5-2683 v4 @ 2.10GHz
- Stepping:              1
- CPU MHz:               1200.091
- CPU max MHz:           3000.0000
- CPU min MHz:           1200.0000
- BogoMIPS:              4190.17
- Virtualization:        VT-x
- L1d cache:             32K
- L1i cache:             32K
- L2 cache:              256K
- L3 cache:              40960K
- NUMA node0 CPU(s):     0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30
- NUMA node1 CPU(s):     1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31
- Flags: fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe syscall nx pdpe1gb rdtscp lm constant_tsc arch_perfmon pebs bts rep_good nopl xtopology nonstop_tsc aperfmperf eagerfpu pni pclmulqdq dtes64 monitor ds_cpl vmx smx est tm2 ssse3 sdbg fma cx16 xtpr pdcm pcid dca sse4_1 sse4_2 x2apic movbe popcnt tsc_deadline_timer aes xsave avx f16c rdrand lahf_lm abm 3dnowprefetch epb cat_l3 cdp_l3 invpcid_single intel_pt tpr_shadow vnmi flexpriority ept vpid fsgsbase tsc_adjust bmi1 hle avx2 smep bmi2 erms invpcid rtm cqm rdt_a rdseed adx smap xsaveopt cqm_llc cqm_occup_llc cqm_mbm_total cqm_mbm_local dtherm ida arat pln pts
