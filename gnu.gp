set term post eps enh color solid
set output 'fit.CF2lysMe3Gd.eps'
set title 'CF2lysMe3Gd'
set xrange[-70.000000:-105.000000]
unset ytics
set xlabel'ppm'
set format y ""
set border 1
set label sprintf("ppm: -99.65,-98.01 R2: 313 J: 247 I: 1.00") at graph 0.02,0.350000 font "Arial,10"
set label sprintf("ppm: -98.62,-97.99 R2: 298 J: 239 I: 0.59") at graph 0.02,0.300000 font "Arial,10"
set label sprintf("ppm: -75.37 R2: 216 I: 0.09") at graph 0.02,0.250000 font "Arial,10"
set label sprintf("ppm: -101.27 R2: 505 I: 0.12") at graph 0.02,0.200000 font "Arial,10"
set label sprintf("ppm: -98.95 R2: 528 I: 0.11") at graph 0.02,0.150000 font "Arial,10"
plot 'outy.out' i 0 ti 'raw' w li,'' i 0 u 1:3 ti 'fit'w li,'' i 1 u 1:($3-50000000000.000000) noti w li lc 3,'' i 2 u 1:($3-100000000000.000000) noti w li lc 3,'' i 3 u 1:($3-150000000000.000000) noti w li lc 3,'' i 4 u 1:($3-200000000000.000000) noti w li lc 3,'' i 5 u 1:($3-250000000000.000000) noti w li lc 3