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
            description="Plot two bands strucutres in one.")

    parser.add_argument(
            "infile", type=str,
            help="input file to calculate spintexture from (fdf-file)")

    parser.add_argument("-o", "--outfile", metavar="outfile", type=str,
                        help="output file",
                        default="out.spin.bands")

    parser.add_argument(
            "-k", "--kpoints", type=valid_klist, required=True,
            help=("corner points of k-path in reciprocal coordinates;"
                  " coordinates space-separated and kpoints comma-separated."))

    parser.add_argument("-n", "--nkpoints", type=int, required=True,
                        help="number of k-points along the path")

    parser.add_argument("-l", "--klabels", type=str, nargs='*',
                        help="labels for corner points of k-path",
                        default=None)

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    H = sisl.get_sile(args.infile).read_hamiltonian()
    # TODO labels
    kpath = sisl.BandStructure(
            H,
            point=args.kpoints,
            division=args.nkpoints,
            name=args.klabels
    )

    nk = args.nkpoints
    nbands = len(H)

    lk, kt, kl = kpath.lineark(True)

    def wrap(es, parent, k, weight):
        return np.concatenate((es.eig.reshape(-1, 1), es.spin_moment()),
                              axis=-1)

    data = kpath.asarray().eigenstate(wrap=wrap, eta=True)

    enc = "latin1"
    with open(args.outfile, "wb") as out:
        out.write("# Spintexture data\n".encode(enc))
        out.write(f"#     nbnds, nk = {nbands}, {nk}\n".encode(enc))
        out.write(f"#     lkmin, lkmax = {min(lk)}, {max(lk)}\n".encode(enc))
        out.write("#\n".encode(enc))
        out.write("# Ticks and labels\n".encode(enc))
        for i in range(len(kt)):
            out.write(f"#     tick {i}: {kt[i]}, {kl[i]}\n".encode(enc))
        out.write("#\n".encode(enc))
        out.write("#Postion along Path, Eigenvalue, Sx, Sy, Sz\n".encode(enc))

        for i in range(data.shape[1]):
            np.savetxt(
                outfile,
                np.concatenate((lk.reshape(-1, 1), data[:, i, :]), axis=-1)
                )
            out.write("\n".encode(enc))
