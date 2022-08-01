#!/bin/sh

#only output cpu version, 
./lpp  simple.loci  -o simple_cpu.cc -p

#output cpu and cuda version in differdnt files
./lpp  simple.loci -d cpu,cuda -o simple.cc,simple_cuda.cc -p

#only output cuda version
./lpp  simple.loci  -d cuda -o simple_cuda1.cc -p

#output cpu and cuda version in one file
./lpp  simple.loci  -d cpu,cuda -o simple1.cc -p
