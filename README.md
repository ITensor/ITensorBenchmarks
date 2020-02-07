# ITensorBenchmarks
Benchmarks for the ITensor library.

Installation instructions

First, download and install the ITensor library.

Clone from Github:

 ~ cd ~/software
 ~ git clone git@github.com:ITensor/ITensor.git

Copy the options.mk from this repository into the ITensor repository:

 ~ git clone git@github.com:ITensor/ITensorBenchmarks.git
 ~ cp ~/software/ITensorBenchmarks/options.mk ~/software/ITensor

Build the ITensor library, first adding the modules:

 ~ module load gcc
 ~ module load llvm
 ~ module load intel/mkl
 ~ cd ~/software/ITensor
 ~ make

You can run the unittests to check that the installation is correct:

 ~ cd ~/software/ITensor/unittest
 ~ make

Now to run benchmark 1, enter the example1 directory:

 ~ cd ~/software/ITensorBenchmark/example1

First modify the Makefile to point to the location of the ITensor libary (in this example, ~/software/ITensor), and then build:

 ~ make

You can try setting the MKL threads, but I have not found that it helps so I usually just set it to 1 (export MKL_NUM_THREADS=1).

