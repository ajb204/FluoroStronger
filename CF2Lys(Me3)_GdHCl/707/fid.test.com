#!/bin/csh
bruk2pipe -in ./fid  \
-bad 0.0 -ext -aswap -AMX -decim 152 -dspfvs 21 -grpdly 76 -ws 8 -noi2f \
 -xN    131072     \
 -xT    65536      \
 -xMODE DQD        \
 -xSW   131578.947368 \
 -xOBS  564.491615211 \
 -xCAR  -99.9999999968 \
 -xLAB  F          \
 -ndim  1          -aq2D  States     \
  -out test.fid -verb -ov
