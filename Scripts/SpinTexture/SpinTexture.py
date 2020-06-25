#!/usr/bin/env python3
import sys
import sisl
import numpy as np
np.set_printoptions(threshold=sys.maxsize)


def parse_args():
    import argparse

    def valid_klist(s):
        try:
            klist = np.asarray(
                    [[float(a)
                      for a in point.strip().split(' ')]
                     for point in s.split(',')]
            )
            assert klist.shape[1] == 3
        except ValueError:
            msg = "Not a valid range: '{0}'.".format(s)
            raise argparse.ArgumentTypeError(msg)
        return klist

    parser = argparse.ArgumentParser(
            description="Calculate the spin texture of a system.")

    parser.add_argument(
            "infile", type=str,
            help="input file to calculate spintexture from (fdf-file)")

    parser.add_argument("-o", "--outfile", metavar="outfile", type=str,
                        help="output file",
                        default="out.spin.bands")

    parser.add_argument(
            "-k", "--kpoints", type=valid_klist, required=True,
            help=("corner points of k-path in reciprocal coordinates;"
                  " coordinates space-separated and kpoints comma-separated"))

    parser.add_argument("-n", "--division", type=int, nargs='+', required=True,
                        help="number of k-points along the path")

    parser.add_argument("-l", "--klabels", type=str, nargs='*',
                        help="labels for corner points of k-path",
                        default=None)

    parser.add_argument("--hsx", action='store_true',
                        help="shift fermi level when reading from HSX")

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    sile = sisl.get_sile(args.infile)
    geom = sile.read_geometry()
    H = sile.read_hamiltonian(geometry=geom)

    assert H is not None, "Could not read Hamiltonian."
    
    if len(args.division) == 1: args.division = args.division[0]
    kpath = sisl.BandStructure(
            H,
            point=args.kpoints,
            division=args.division,
            name=args.klabels,
    )

    def wrap(es, parent, k, weight):
        return np.concatenate(
            (es.eig.reshape(-1, 1), es.spin_moment()),
            axis=-1,
        )

    data = kpath.asarray().eigenstate(wrap=wrap, eta=True)
    if args.hsx:
        Ef = sile.read_fermi_level()
        data[:, :, 0] -= Ef

    nk, nbands, _ = data.shape

    lk = kpath.lineark()
    xtick, xtick_label = kpath.lineartick()

    enc = "latin1"
    with open(args.outfile, "wb") as out:
        out.write("# Spintexture data\n".encode(enc))
        out.write("#     nk, nbnds = {0}, {1}\n".format(nk, nbands).encode(enc))
        out.write("#     lkmin, lkmax = {0}, {1}\n".format(min(lk), max(lk)).encode(enc))
        out.write("#\n".encode(enc))
        out.write("# Ticks and labels\n".encode(enc))
        for i in range(len(xtick)):
            out.write("#     tick {0}: {1}, {2}\n".format(i, xtick[i], xtick_label[i]).encode(enc))
        out.write("#\n".encode(enc))
        out.write("#Postion along Path, Eigenvalue, Sx, Sy, Sz\n".encode(enc))

        for i in range(data.shape[1]):
            np.savetxt(
                out,
                np.concatenate((lk.reshape(-1, 1), data[:, i, :]), axis=-1)
                )
            out.write("\n".encode(enc))
