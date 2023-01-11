#!/usr/bin/env python3
import numpy as np
import sisl as s
from sisl.quaternion import Quaternion
from sisl.utils.mathematics import cart2spher, fnorm
import sisl._array as _a


def parse_args():
    import argparse

    parser = argparse.ArgumentParser(
            description="Calculate the spin texture of a system.")

    parser.add_argument(
            "infile", type=str,
            help="input file. json formatted file with calculation parameters")

    return parser.parse_args()

def parse_json(filename):
    import json
    with open(filename, 'r') as json_file:
        json_data = json.load(json_file)
    return json_data


def read_SpinBands(filename):
    with open(filename, 'r') as f:
        content = f.readlines()

    kt = []   # tick markers along the path 
    kl = []   # labels for each tick
    ibnd = 0  # band index
    ik = 0    # k point index

    header=True
    last_line_empty=False
    for line in content:
        if line.strip().startswith("#"):
            # Lines starting with # denote the header section of the file
            if "nbnds" in line:
                # This line contains dimension
                nk = int(line.split("=")[1].split(",")[0].strip())
                nbnds = int(line.split("=")[1].split(",")[1].strip())
            if "lkmin" in line:
                lkmin = float(line.split("=")[1].split(",")[0].strip())
                lkmax = float(line.split("=")[1].split(",")[1].strip())
            if "tick" in line and "Ticks and labels" not in line:
                kt.append(float(line.split(":")[1].split(",")[0].strip()))
                kl.append(line.split(":")[1].split(",")[1].strip())                
            
        elif line.strip() == "":
            # An empty line indicates a new band
            if last_line_empty:
                last_line_empty = False
            else:
                ibnd += 1
                ik = 0
                last_line_empty = True

        else:            
            if header:
                # Just after we are done reading the header we will have to check
                # whether this file in 'bands' format, i.e., the first
                # column is the position along a path OR wether the 
                # first 3 column contain k point coorindates
                header = False
                dim = len(line.split())
                data = np.zeros((nk, nbnds, dim))
            # Read the position along the path or the kpoint and the
            # corresponding energies / moments
            data[ik, ibnd] = np.asarray([float(x) for x in line.split()])
            ik += 1
    
    return data[:,0,0:-4], data[...,-4:], kt, kl


def cleanfortex(s):
    clean_s = s.replace("_", "\\textunderscore")
    return clean_s


if __name__ == '__main__':
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    cl_args = parse_args()
    json_args = parse_json(cl_args.infile)
    save = json_args['outfile'] != ""

    if save:
        matplotlib.use('pdf')
        plt.rcParams['text.usetex'] = True
        plt.rcParams['text.latex.preamble'] = r'\usepackage{lmodern}'
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = 'Computer Modern'
        plt.rcParams['font.size'] = 16

    # Open read input files
    print('Reading file: {}'.format(json_args['infile']))
    k, bands, xtick, xtick_label = read_SpinBands(json_args['infile'])

    if json_args['mode'] == 'bands':
        if k.shape[1] == 3:
            k[1:]-= k[:-1]                 # k[i] = k[i] - k[i-1]
            dk = np.linalg.norm(k, axis=1) # sqrt(kx**2 + ky**2 + kz**2)
            dk[0] = 0 
            lk = np.cumsum(dk)
        else:
            lk = k[:,0]

        # Select bands in the plotting range
        if json_args.get('yrange', None) is not None:
            bands = bands[:,  (np.max(bands[:, :, 0], 0) > json_args['yrange'][0])
                            * (np.min(bands[:, :, 0], 0) < json_args['yrange'][1])]
        

        # Plotting
        if json_args.get('no-spin', False):
            # No spins, only bands
            fig, ax = plt.subplots(1, 1, figsize=(4, 5),
                                   dpi=json_args.get('dpi',None),
                                   sharex=True, sharey=True)
            axes = np.array([ax], dtype=object)

            for ibnd in range(bands.shape[1]):
                axes[0].plot(lk, bands[:, ibnd,0], color='black', linewidth=2)
        else:
            # With spin
            if json_args.get('only', None) is None:
                components = [(0, '$S_x$'),
                            (1, '$S_y$'),
                            (2, '$S_z$')]
                fig, axes = plt.subplots(1, 3, figsize=(9, 6),
                                       dpi=json_args.get('dpi',None),
                                       sharex=True, sharey=True)
            else:
                # Only one component
                components = {'x': [(0, '$S_x$')],
                            'y': [(1, '$S_y$')],
                            'z': [(2, '$S_z$')]}.get(json_args['only'])
                fig, ax = plt.subplots(1, 1, figsize=(4, 5),
                                       dpi=json_args.get('dpi',None),
                                       sharex=True, sharey=True)
                axes = np.array([ax], dtype=object)

            # Create the colored lines
            norm = plt.Normalize(-1, 1)
            lw = 2
            for iaxes, (icomp, component) in enumerate(components):
                for ibnd in range(bands.shape[1]):
                    # Define pairs of x and y values
                    points = np.array([lk, bands[:, ibnd, 0]]).T.reshape(-1, 1, 2)
                    # Define segments from each pair to the next
                    segments = np.concatenate([points[:-1], points[1:]], axis=1)
                    # Create a collection of the segements
                    lc = LineCollection(segments, cmap='coolwarm', norm=norm)
                    # Use the spin moment to define the color
                    lc.set_array(bands[:, ibnd, 1+icomp])
                    # Adjust the linewidth
                    lc.set_linewidth(lw)
                    # Add the collection to the plot
                    line = axes[iaxes].add_collection(lc)
                # Label the subplot
                axes[iaxes].set_title(component)

        if json_args.get('no-spin', False):
            plt.tight_layout()
        elif 'only' in json_args:
            fig.subplots_adjust(right=0.7)
            pos = axes[0].get_position() # get the original position
            height = pos.y1 - pos.y0
            cbar_ax = fig.add_axes([0.75, pos.y0, 0.05, height])
            fig.colorbar(line, cax=cbar_ax)
            if 'ylabel' in json_args:
                fig.subplots_adjust(left=0.2)
        else:
            fig.subplots_adjust(right=0.8)
            pos = axes[0].get_position() # get the original position
            height = pos.y1 - pos.y0
            cbar_ax = fig.add_axes([0.85, pos.y0, 0.025, height])
            fig.colorbar(line, cax=cbar_ax)
            # if 'ylabel' in json_args:
            #     fig.subplots_adjust(left=0.2)

        # Fix the axis labels and ranges
        # All subplots share the same axis settings, so we can just them once    
        if 'xrange' not in json_args:
            json_args['xrange'] = (min(lk), max(lk)) 
        axes[0].set(
                xlim=json_args.get('xrange', None),
                ylim=json_args.get('yrange', None))

        axes[0].set_ylabel(json_args.get('ylabel', None))
        axes[0].set_xlabel(json_args.get('xlabel', None))
        axes[0].locator_params(axis='y', nbins=5)

        # Plot vertical lines at high symmetry points
        if len(xtick) > 0:
            lbound = np.where(np.greater_equal(xtick, json_args['xrange'][0]))[0][0]
            ubound = np.where(np.less_equal(xtick, json_args['xrange'][1]))[0][-1] + 1
            axes[0].xaxis.set_ticks(xtick[lbound:ubound])
            axes[0].set_xticklabels(xtick_label[lbound:ubound])
            for axis in axes:
                for tick in xtick[lbound:ubound]:
                    axis.axvline(x=tick, color='k', linewidth=0.5)

    else:
        # Plot the data from files
        if json_args.get('erange',None):
            selected_bands = np.where( (np.mean(bands[:, :, 0], 0) > json_args['erange'][0])
                                      *(np.mean(bands[:, :, 0], 0) < json_args['erange'][1]))[0]
        else:
            selected_bands = np.arange(bands.shape[1])

        # Generate rotation
        k_n = _a.asarrayd(json_args.get('normal', [0., 0., 1.]))
        _, phi, theta = cart2spher(k_n)
        if theta != 0:
            q = Quaternion(-theta, [0, 1, 0], rad=True) * Quaternion(-phi, [0, 0, 1], rad=True)
        else:
            # No rotation required
            q = Quaternion(0., [0, 1, 0], rad=True)

        # Project path onto plane
        k = q.rotate(k)

        nbnds = selected_bands.shape[0]
        nrows = int(np.sqrt(nbnds))
        ncols = (nbnds + nrows - 1) // nrows

        fig, axes = plt.subplots(nrows, ncols, figsize=(4.5*ncols+2.0, 4.5*nrows+1.5),
                                 sharex=True, sharey=True)
        if nrows == 1:
            axes = axes.reshape(1,-1)
        for iax, ibnd in enumerate(selected_bands):
            row = iax // ncols
            col = iax  % ncols
            axes[row,col].set_title(f'Band: {ibnd} - $\langle E \\rangle$={np.mean(bands[:,ibnd,0]):.3f} eV')

            # Spin moments onto plane
            bands[:, ibnd, 1:] = q.rotate(bands[:, ibnd, 1:])

            def color(x):
                y = np.empty_like(x, dtype=object)
                y[np.where(x < 0)] = "tab:blue"
                y[np.where(x == 0)] = "black"
                y[np.where(0 < x)] = "tab:red"
                return y
            if json_args.get('line', False):
                axes[row,col].plot(
                    k[:,0], k[:,1],
                    color='black',
                    linewidth=1
                )
            else:
                axes[row,col].scatter(
                    k[:,0], k[:,1],
                    color='black',
                )
            axes[row,col].quiver(
                k[:,0], k[:,1],
                bands[:, ibnd, 1],bands[:, ibnd, 2],
                pivot='tail',
                color="tab:orange",
                scale=3,
            )
            
        # All subplots share the same axis settings, so we can just them once
        axes[0,0].set(xlim=json_args.get('xrange', None), ylim=json_args.get('yrange', None))
        for ax in axes.ravel():
            ax.set(aspect=1)
        fig.supxlabel(json_args.get('xlabel', '$k_x$ [1/Å]'))
        fig.supylabel(json_args.get('ylabel', '$k_y$ [1/Å]'))
        plt.tight_layout()

    for axis in axes:
        if json_args.get('grid',False):
            axis.grid()

    if save:
        plt.savefig(json_args['outfile'])
    else:
        plt.show()
