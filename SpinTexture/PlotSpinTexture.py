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
                nbnds = int(line.split("=")[1].split(",")[0].strip())
                nk = int(line.split("=")[1].split(",")[1].strip())
                data = np.zeros((nbnds, nk, 5))
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
            data[ibnd, ik] = np.asarray([float(x) for x in line.split()])
            ik += 1
    lk = np.linspace(lkmin, lkmax, nk)

    return (kt, kl), lk, data


def cleanfortex(s):
    clean_s = s.replace("_", "\\textunderscore")
    return clean_s


if __name__ == '__main__':
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import ListedColormap, BoundaryNorm

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

    # Initialize Figure
    fig, axs = plt.subplots(1, 3, figsize=(9, 6), sharex=True, sharey=True)

    # Open read input files
    handles = []
    labels = []
    first = True
    print('Reading file: {}'.format(args.bandsfile))
    xtick, lk, bands = read_SpinBands(args.bandsfile)
    nbnds = bands.shape[0]
    nk = bands.shape[1]

    args.xrange = (min(lk), max(lk)) if args.xrange is None else args.xrange

    norm = plt.Normalize(-1, 1)
    # Plot the data from files
    first = True
    for idir, dir in enumerate(['x', 'y', 'z']):
        for ibnd in range(nbnds):
            if (np.max(bands[ibnd, :, 1]) < args.yrange[0]
                    or np.min(bands[ibnd, :, 1]) > args.yrange[1]):
                continue
            points = np.array([lk, bands[ibnd, :, 1]]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lc = LineCollection(segments, cmap='coolwarm', norm=norm)
            lc.set_array(bands[ibnd, :, 2+idir])
            lc.set_linewidth(2)
            line = axs[idir].add_collection(lc)
            if first:
                fig.colorbar(line, ax=axs[-1])
                first = False

        axs[idir].set_title(f'$S_{dir}$')
        axs[idir].set_xlim(args.xrange[0], args.xrange[1])
        axs[idir].set_ylim(args.yrange[0], args.yrange[1])

        # Plot vertical lines at high symmetry points
        lbound = np.where(np.greater_equal(xtick[0], args.xrange[0]))[0][0]
        ubound = np.where(np.less_equal(xtick[0], args.xrange[1]))[0][-1] + 1
        for tick in xtick[0][lbound:ubound]:
            axs[idir].plot([tick, tick], args.yrange, 'k', linewidth=0.5)

        axs[idir].xaxis.set_ticks(xtick[0][lbound:ubound])
        axs[idir].set_xticklabels(xtick[1][lbound:ubound])

    if save:
        plt.savefig(args.outfile)
    else:
        plt.show()
