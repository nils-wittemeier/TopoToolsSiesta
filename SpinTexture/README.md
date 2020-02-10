
# 1. Calculate the spin texture (SpinTexture.py):
Spin texture can be calculated from siesta outputs if either the HSX or
TSHS file was generated. The TSHS-file is sufficient on its on. When using
the HSX-file more output files are required.

Required files:
    TSHS
        or
    HSX + fdf + ion(.nc)

The infile argument can be either the fdf-file or the TSHS-file.

## Usage:
usage: SpinTexture.py [-h] [-o outfile] -k KPOINTS -n NKPOINTS
                      [-l [KLABELS [KLABELS ...]]]
                      infile

Plot two bands strucutres in one.

positional arguments:
  infile                input file to calculate spintexture from (TSHS/fdf-file)

optional arguments:
  -h, --help            show this help message and exit
  -o outfile, --outfile outfile
                        output file
  -k KPOINTS, --kpoints KPOINTS
                        corner points of k-path in reciprocal
                        coordinates;coordinates space-separated,kpoints comma-
                        separated
  -n NKPOINTS, --nkpoints NKPOINTS
                        number of k-points along the path
  -l [KLABELS [KLABELS ...]], --klabels [KLABELS [KLABELS ...]]
                        labels for corner points of k-path

##	Example: 
    python ../SpinTexture.py scf.fdf -o f-hex.spin.bands -k "0.5 0. 0., 0. 0. 0., 0.33333 0.33333 0., 0. 0. 0." -n 300 -l M $\\Gamma$ K M

# 2. Plot spin texture (PlotSpinTexture.py)

Required file:
    spin.bands file genearted by SpinTexture.py

## Usage:
usage: PlotSpinTexture.py [-h] [-o outfile] [-x XRANGE] [-y YRANGE]
                              [--title TITLE] [--xlabel XLABEL] [--ylabel YLABEL]
                              [--with-grid] [--with-legend]
                              [--legend-labels LEGEND_LABELS [LEGEND_LABELS ...]]
                              bandsfile

    Plot the spin texture of a system.

    positional arguments:
      bandsfile             file with spin texture file (spin.bands-file)

    optional arguments:
      -h, --help            show this help message and exit
      -o outfile, --outfile outfile
                            output file
      -x XRANGE, --xrange XRANGE
                            range on x-axis; format "min:max"
      -y YRANGE, --yrange YRANGE
                            range on y-axis; format "min:max"
      --title TITLE         plot tile
      --xlabel XLABEL       label for x-axis
      --ylabel YLABEL       label for y-axis
      --with-grid           enable grid in plot
      --with-legend         enable legend in plot
      --legend-labels LEGEND_LABELS [LEGEND_LABELS ...]
                            labels for each dataset

## Example:
    python ../PlotSpinTexture.py f-hex.spin.bands --yrange="-2:2" -o SpinTexture.svg
