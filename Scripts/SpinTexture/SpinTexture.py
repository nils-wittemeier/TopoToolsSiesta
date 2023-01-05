#!/usr/bin/env python3
import sys
import sisl
import numpy as np
np.set_printoptions(threshold=sys.maxsize)


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

def get_BZ_obj(parent, mode, **kwargs):
    if mode == 'BandStructure':
        return sisl.BandStructure(parent, **kwargs)
    elif mode == 'MonkhorstPack':
        return sisl.MonkhorstPack(parent, **kwargs)
    elif mode == 'BrillouinZone':
        return sisl.BrillouinZone(parent, **kwargs)
    elif mode == 'Circle':
        return sisl.BandStructure.param_circle(parent, **kwargs)
    else:
        assert False, f'Unknown kpath mode "{mode}". Allowed modes are: "BandStructure", "MonkhorstPack", "BrillouinZone", "Circle"'


if __name__ == '__main__':
    cl_args = parse_args()
    json_args = parse_json(cl_args.infile)

    sile = sisl.get_sile(json_args['fdf'])
    geom = sile.read_geometry()
    H = sile.read_hamiltonian(geometry=geom)

    assert H is not None, "Could not read Hamiltonian."
    
    kpath = get_BZ_obj(H, json_args['kpath']['mode'], **json_args['kpath']['kwargs'])


    def wrap(es, parent, k, weight):
        return np.concatenate(
            (es.eig.reshape(-1, 1), es.spin_moment()),
            axis=-1,
        )

    data = kpath.apply.array.renew(zip=False).eigenstate(wrap=wrap, eta=True)

    nk, nbands, _ = data.shape

    try: 
        lk = kpath.lineark()
        xtick, xtick_label = kpath.lineartick()

        enc = "latin1"
        with open(json_args['out'], "wb") as out:
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
    except AttributeError:
        enc = "latin1"
        with open(json_args['out'], "wb") as out:
            out.write("# Spintexture data\n".encode(enc))
            out.write("#     nk, nbnds = {0}, {1}\n".format(nk, nbands).encode(enc))
            out.write("#\n".encode(enc))
            
            out.write("#k point [kx, ky, kz], Eigenvalue, Sx, Sy, Sz\n".encode(enc))

            for i in range(data.shape[1]):
                np.savetxt(
                    out,
                    np.concatenate((kpath.tocartesian(kpath.k), data[:, i, :]), axis=-1)
                    )
                out.write("\n".encode(enc))

