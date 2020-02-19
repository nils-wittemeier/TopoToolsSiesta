#!/bin/sh
mkdir work
cd work;
cp ../Bi.psml .;
cp ../Bi2D_BHex.fdf .;
siesta < Bi2D_BHex.fdf > Bi2D_BHex.out
cd ..
