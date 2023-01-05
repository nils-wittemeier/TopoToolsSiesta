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

    parser.add_argument("-m", "--mode", type=str,
                        help="plot mode; \"lk\" for linear plot (band structure like), \"quiver\" for arrows on grid/set of k points",
                        default="lk")

    parser.add_argument("-x", "--xrange", type=valid_range,
                        help="range on x-axis; format \"min:max\"")

    parser.add_argument("-y", "--yrange", type=valid_range,
                        help="range on y-axis; format \"min:max\"",
                        default=(-3, 3))

    parser.add_argument("-e", "--erange", type=valid_range,
                        help="energy range; format \"min:max\"; only takes effect is \"mode\"=\"quiver\"",
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

    parser.add_argument("--only", type=str,
                        help="only plot one component (x/y/z)")

    parser.add_argument("--no-spin", action="store_true",
                        help="plot no spin component")

    return parser.parse_args()


def read_SpinBands(filename, mode):
    with open(filename, 'r') as f:
        content = f.readlines()

    kt = []
    kl = []
    ibnd = 0
    ik = 0
    first_after_header=True
    for line in content:
        if line.strip().startswith("#"):
            if "nbnds" in line:
                nk = int(line.split("=")[1].split(",")[0].strip())
                nbnds = int(line.split("=")[1].split(",")[1].strip())
                if mode == "lk":
                    data = np.zeros((nk, nbnds, 5))
                else: # mode == "quiver"
                    data = np.zeros((nk, nbnds, 7))
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
    
    if mode == "lk": 
        return data[:,0,0], data[...,1:], kt, kl
    else: # mode = "quiver"
        return data[:,0,0:3], data[...,3:], kt , kl


def cleanfortex(s):
    clean_s = s.replace("_", "\\textunderscore")
    return clean_s


if __name__ == '__main__':
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    import matplotlib.colors as mcolors

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
    k, bands, xtick, xtick_label = read_SpinBands(args.bandsfile, args.mode)

    if args.mode == 'lk':

        # Initialize Figure
        if args.only is None:
            if args.no_spin:
                fig, ax = plt.subplots(1, 1, figsize=(4, 5), sharex=True, sharey=True)
                axes = np.array([ax], dtype=object)
            else:
                components = [(0, '$S_x$'),
                            (1, '$S_y$'),
                            (2, '$S_z$')]
                fig, axes = plt.subplots(1, 3, figsize=(9, 6), sharex=True, sharey=True)
        else:
            components = {'x': [(0, '$S_x$')],
                        'y': [(1, '$S_y$')],
                        'z': [(2, '$S_z$')]}.get(args.only)
            fig, ax = plt.subplots(1, 1, figsize=(4, 5), sharex=True, sharey=True)
            axes = np.array([ax], dtype=object)

        # Plot the data from files
        selected_bands = bands[:,  (np.max(bands[:, :, 0], 0) > args.yrange[0])
                                 * (np.min(bands[:, :, 0], 0) < args.yrange[1])]
        if args.no_spin:
            for ibnd in range(selected_bands.shape[1]):
                axes[0].plot(lk, bands[:, ibnd,0], color='black', linewidth=2)
        else:
            norm = plt.Normalize(-1, 1)
            lw = 2
            for iaxes, (icomp, component) in enumerate(components):
                for ibnd in range(bands.shape[1]):
                    if (np.max(bands[:, ibnd, 0]) < args.yrange[0]
                            or np.min(bands[:, ibnd, 0]) > args.yrange[1]):
                        continue
                    if mode == 'lk':
                        points = np.array([lk, bands[:, ibnd, 0]]).T.reshape(-1, 1, 2)
                    else:
                        points = np.array([k[:,0], bands[:, ibnd, 0]]).T.reshape(-1, 1, 2)
                    segments = np.concatenate([points[:-1], points[1:]], axis=1)
                    lc = LineCollection(segments, cmap='coolwarm', norm=norm)
                    lc.set_array(bands[:, ibnd, 1+icomp])
                    lc.set_linewidth(lw)
                    line = axes[iaxes].add_collection(lc)
                axes[iaxes].set_title(component)
            fig.colorbar(line, ax=axes.ravel().tolist())
        
        # Plot vertical lines at high symmetry points
        args.xrange = (min(lk), max(lk)) if args.xrange is None else args.xrange
        lbound = np.where(np.greater_equal(xtick, args.xrange[0]))[0][0]
        ubound = np.where(np.less_equal(xtick, args.xrange[1]))[0][-1] + 1
        axes[0].xaxis.set_ticks(xtick[lbound:ubound])
        axes[0].set_xticklabels(xtick_label[lbound:ubound])
        for axis in axes:
            for tick in xtick[lbound:ubound]:
                axis.plot([tick, tick], args.yrange, 'k', linewidth=0.5)
                    
        # All subplots share the same axis settings, so we can just them once    
        axes[0].set_xlim(*args.xrange)
        axes[0].set_ylim(*args.yrange) 
        axes[0].set_ylabel('Eigenspectrum [eV]')
        axes[0].locator_params(axis='y', nbins=5)

    else:
        # Plot the data from files
        selected_bands = bands[:,  (np.max(bands[:, :, 0], 0) > args.erange[0])
                                 * (np.min(bands[:, :, 0], 0) < args.erange[1])]
        nbnds = selected_bands.shape[1]
        nrows = int(np.sqrt(nbnds))
        ncols = (nbnds + nrows - 1) // nrows
        fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 5*nrows),
                                 sharex=True, sharey=True)
        if nrows == 1:
            axes = axes.reshape(1,-1)
        for ibnd in range(selected_bands.shape[1]):
            row = ibnd // ncols
            col = ibnd  % ncols
            axes[row,col].set_title(f'$\langle E \\rangle$={np.mean(selected_bands[:,ibnd,0]):.3f} eV', fontsize=18)
            def color(x):
                y = np.empty_like(x, dtype=object)
                y[np.where(x < 0)] = "tab:blue"
                y[np.where(x == 0)] = "black"
                y[np.where(0 < x)] = "tab:red"
                return y
            axes[row,col].plot(
                k[:,0], k[:,1],
                color='black', linewidth=1)
            axes[row,col].scatter(
                k[:,0], k[:,1],
                s=np.abs(selected_bands[:, ibnd, 3])*100,
                c=color(selected_bands[:,ibnd, 3]),
                label="($S_x$,$S_y$) [arb. units]"
                )#, color='black', linewidth=2)
            axes[row,col].quiver(
                k[:,0], k[:,1],
                selected_bands[:, ibnd, 1],selected_bands[:, ibnd, 2],
                pivot='tail',
                scale=4,
                #  units='xy', scale = 10,
                color="tab:orange"
                )#, color='black', linewidth=2)
            # axes[iaxes].scatter(
            #     np.roll(k,icomp,1)[:,0], selected_bands[:,ibnd,0],
            #     )#, color='black', linewidth=2)
        
        args.xrange = (min(k[:,0]), max(k[:,0])) if args.xrange is None else args.xrange
            
        # All subplots share the same axis settings, so we can just them once
        axes[0,0].set(xlim=args.xrange, ylim=args.yrange)
        for ax in axes.ravel():
            ax.set(aspect=1)
        fig.supxlabel('$k_x$ [1/Å]')
        fig.supylabel('$k_y$ [1/Å]')
        # plt.legend()

    for axis in axes:
        if args.with_grid:
            axis.grid()
    plt.tight_layout()

    if save:
        plt.savefig(args.outfile)
    else:
        plt.show()
