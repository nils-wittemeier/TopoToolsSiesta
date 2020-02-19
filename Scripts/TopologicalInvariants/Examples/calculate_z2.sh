#!/bin/sh
python3 ../../Z2run.py 'work/Bi2D_BHex.fdf' -o Bi2D_BHex.spin.bands \
    -k "0.5 0. 0., 0. 0. 0., 0.33333 0.33333 0., 0. 0. 0." -n 300 \
    -l M $\\Gamma$ K M \
    --hsx
python3 ../../PlotSpinTexture.py Bi2D_BHex.spin.bands --yrange="-2:2" \
    -o SpinTexture.svg
