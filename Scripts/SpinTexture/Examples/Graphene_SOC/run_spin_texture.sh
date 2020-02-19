#!/bin/sh
python3 ../../SpinTexture.py 'work/in.fdf' -o graphene.spin.bands \
    -k "0.3333400000000000 0.3333200000000000 0., \
        0.3333333333333333 0.3333333333333333 0., \
        0.3333266666666667 0.3333266666666667 0." \
    -n  5 5 \
    -l M K $\\Gamma$\' \
    --hsx
python3 ../../PlotSpinTexture.py graphene.spin.bands \
    --yrange="0.1629:0.16325" \
    -o SpinTexture.svg

# 163325
#    25
# 163075
#    25
# 162825
# 
# 16320
#    25
# 16295
#    25
# 16270
