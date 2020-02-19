#!/usr/bin/env python3
"""Bandstrucutre plot for comaprison
===================================

usage: compare_bands.py [-h] [-x XRANGE] [-y YRANGE] [--title TITLE]
                        [--xlabel XLABEL] [--ylabel YLABEL] [--with-grid]
                        [--with-legend]
                        bandsfile1 bandsfile2 [outfile]

Plot two bands strucutres in one.

positional arguments:
  bandsfiles            siesta .bands files to be compared
  outfile               output file

optional arguments:
  -h, --help            show this help message and exit

  -x XRANGE, --xrange XRANGE
                        range on x-axis; format "min:max"
  -y YRANGE, --yrange YRANGE
                        range on y-axis; format "min:max"
  --title TITLE         plot tile
  --xlabel XLABEL       label for y-axis
  --ylabel YLABEL       label for y-axis
  --with-grid           enable grid in plot
  --with-legend         enable legend in plot

"""

import numpy as np
import sisl as s


def parse_args():
    # parse command line arguments
    import argparse

    def valid_range(s):
        # Validity test for ranges parsed from command line
        try:
            return tuple(float(a) for a in s.split(':'))
        except ValueError:
            msg = "Not a valid range: '{0}'.".format(s)
            raise argparse.ArgumentTypeError(msg)

    parser = argparse.ArgumentParser(
            description="Plot the spin texture of a system.")

    parser.add_argument("bandsfile", type=str,
                        help="file with spin texture file (spin.bands-file)")

    parser.add_argument("-o", "--outfile", metavar="outfile", type=str,
                        help="output file", default="")

    parser.add_argument("-x", "--xrange", type=valid_range,
                        help="range on x-axis; format \"min:max\"")

    parser.add_argument("-y", "--yrange", type=valid_range,
                        help="range on y-axis; format \"min:max\"",
                        default=(-3, 3))

    parser.add_argument("--title", type=str, help="plot tile", default='')

    parser.add_argument("--xlabel", type=str, help="label for x-axis",
                        default='')

    parser.add_argument("--ylabel", type=str, help="label for y-axis",
                        default='Energy (eV)')

    parser.add_argument("--with-grid", action="store_true",
                        help="enable grid in plot")

    parser.add_argument("--with-legend", action="store_true",
                        help="enable legend in plot")

    parser.add_argument("--legend-labels", nargs='+', type=str,
                        help="labels for each dataset")

    return parser.parse_args()


def read_SpinBands(filename):
    with open(filename, 'r') as f:
        content = f.readlines()

    kt = []
    kl = []
    ibnd = 0
    ik = 0
    for line in content:
        if line.strip().startswith("#"):
            if "nbnds" in line:
                nk = int(line.split("=")[1].split(",")[0].strip())
                nbnds = int(line.split("=")[1].split(",")[1].strip())
                data = np.zeros((nk, nbnds, 5))
            if "lkmin" in line:
                lkmin = float(line.split("=")[1].split(",")[0].strip())
                lkmax = float(line.split("=")[1].split(",")[1].strip())
            if "tick" in line and "Ticks and labels" not in line:
                kt.append(float(line.split(":")[1].split(",")[0].strip()))
                kl.append(line.split(":")[1].split(",")[1].strip())
        elif line.strip() == "":
            ibnd += 1
            ik = 0
        else:
            data[ik, ibnd] = np.asarray([float(x) for x in line.split()])
            ik += 1
    
    return data[:,0,0], data[...,1:], kt, kl


def cleanfortex(s):
    clean_s = s.replace("_", "\\textunderscore")
    return clean_s


if __name__ == '__main__':
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    args = parse_args()
    save = args.outfile != ""

    strSpin = ['UP', 'DOWN']
    ls = np.asarray([['-'], ['-'], ['-'], ['-'], ['-']])
    mk = np.asarray([[''], [''], [''], [''], ['']])
    lw = np.asarray([[1.0], [1.0], [1.0], [1.0], [1.0]])
    ms = np.asarray([[1.0], [1.0], [1.0], [1.0], [1.0]])
    # colors = np.asarray([['b'], ['orange'], ['r'], ['green']])
    colors = np.asarray([
        ['#332288'],
        ['#88CCEE'],
        ['#44AA99'],
        ['#117733'],
        ['#999933'],
        ['#DDCC77'],
        ['#CC6677'],
        ['#882255'],
        ['#AA4499']
    ])

    if save:
        matplotlib.use('pdf')
        plt.rcParams['text.usetex'] = True
        plt.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern}']
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Computer Modern']
        plt.rcParams['font.size'] = 24

    # Open read input files
    print('Reading file: {}'.format(args.bandsfile))
    lk, bands, xtick, xtick_label = read_SpinBands(args.bandsfile)

    # Initialize Figure
    fig, axes = plt.subplots(1, 3, figsize=(9, 6), sharex=True, sharey=True)
    norm = plt.Normalize(-1, 1)
    lw = 2
    # Plot the data from files
    for icomp, component in enumerate(['$S_x$', '$S_y$', '$S_z$']):
        for ibnd in range(bands.shape[1]):
            if (np.max(bands[:, ibnd, 0]) < args.yrange[0]
                    or np.min(bands[:, ibnd, 0]) > args.yrange[1]):
                continue
            points = np.array([lk, bands[:, ibnd, 0]]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lc = LineCollection(segments, cmap='coolwarm', norm=norm)
            lc.set_array(bands[:, ibnd, 1+icomp])
            lc.set_linewidth(lw)
            line = axes[icomp].add_collection(lc)
        axes[icomp].set_title(component)
        
    # All subplots share the same axis settings, so we can just them once    
    args.xrange = (min(lk), max(lk)) if args.xrange is None else args.xrange
    axes[0].set_xlim(*args.xrange)
    axes[0].set_ylim(*args.yrange) 
    axes[0].set_ylabel('Eigenspectrum [eV]')
        
    # Plot vertical lines at high symmetry points
    lbound = np.where(np.greater_equal(xtick, args.xrange[0]))[0][0]
    ubound = np.where(np.less_equal(xtick, args.xrange[1]))[0][-1] + 1
    axes[0].xaxis.set_ticks(xtick[lbound:ubound])
    axes[0].set_xticklabels(xtick_label[lbound:ubound])
    for axis in axes:
        for tick in xtick[lbound:ubound]:
            axis.plot([tick, tick], args.yrange, 'k', linewidth=0.5)

    fig.colorbar(line, ax=axes.ravel().tolist())

    if save:
        plt.savefig(args.outfile)
    else:
        plt.show()
