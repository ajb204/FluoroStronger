#!/bin/csh
nmrPipe -in test.fid \
| nmrPipe  -fn EM  -lb 15.000000 -c 0.5              \
| nmrPipe  -fn ZF -auto                        \
| nmrPipe  -fn FT -auto                        \
| nmrPipe  -fn PS -p0 173.798760 -p1 96.768563 -di -verb  \
| nmrPipe -fn EXT -x1 -111.000000ppm -xn -65.000000ppm -sw     \
| nmrPipe -fn BASE -nw 5 -nl -65.000000ppm -70.000000ppm -85.000000ppm -90.000000ppm -108.000000ppm -111.000000ppm   \
   -ov -out test.ft2
