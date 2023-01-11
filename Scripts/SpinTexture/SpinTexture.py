#!/usr/bin/env python3
import numpy as np
from numbers import Integral
import sisl
from sisl.quaternion import Quaternion
from sisl.utils.mathematics import cart2spher, fnorm
import sisl._array as _a
from scipy.optimize import root_scalar
import sys
from tqdm.auto import trange

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
        assert False, f'Unknown kpath mode \"{mode}\". Allowed modes are: "BandStructure", "MonkhorstPack", "BrillouinZone", "Circle"'

class KPointLoop(sisl.BrillouinZone):
    def __init__(self, parent, N_or_dk, kR, normal, origin, loop=False):
        super().__init__(parent)

        if isinstance(N_or_dk, Integral):
            N = N_or_dk
        else:
            # Calculate the required number of points
            N = int(kR ** 2 * np.pi / N_or_dk + 0.5)
            if N < 4:
                N = 4
                info('BrillouinZone.param_circle increased the number of circle points to 4.')

        normal = _a.asarrayd(normal)
        origin = _a.asarrayd(origin)
        k_n = self.tocartesian(normal)
        self._k_o = self.tocartesian(origin)

        # Generate a preset list of k-points on the unit-circle
        if loop:
            radians = _a.aranged(N) / (N-1) * 2 * np.pi
        else:
            radians = _a.aranged(N) / N * 2 * np.pi
        k = _a.emptyd([N, 3])
        k[:, 0] = np.cos(radians)
        k[:, 1] = np.sin(radians)
        k[:, 2] = 0.

        # Now generate the rotation
        _, theta, phi = cart2spher(k_n)
        if theta != 0:
            pv = _a.arrayd([k_n[0], k_n[1], 0])
            pv /= fnorm(pv)
            q = Quaternion(phi, pv, rad=True) * Quaternion(theta, [0, 0, 1], rad=True)
        else:
            q = Quaternion(0., [0, 0, k_n[2] / abs(k_n[2])], rad=True)

        # Calculate k-points
        k = q.rotate(k)
        self._k0 = k / fnorm(k).reshape(-1, 1)
        self._k = self.toreduced(kR*self._k0 + self._k_o)

        # The sum of weights is equal to the BZ area
        W = np.pi * kR ** 2
        self._w = np.repeat([W / N], N)

    def get_radii(self):
        k = self.tocartesian(self._k) - self._k_o
        return fnorm(k)

    def update_radius(self, ik, radius):
        self._k[ik] = self.toreduced(radius * self._k0[ik] + self._k_o)

    def update_radii(self, radii):
        assert len(self) == len(radii)
        self._k = self.toreduced(radii.reshape((-1,1)) * self._k0 + self._k_o)


if __name__ == '__main__':
    cl_args = parse_args()
    json_args = parse_json(cl_args.infile)

    sile = sisl.get_sile(json_args['fdf'])
    geom = sile.read_geometry()
    H = sile.read_hamiltonian(geometry=geom)

    def wrap(es, parent, k, weight):
        return np.concatenate(
            (es.eig.reshape(-1, 1), es.spin_moment()),
            axis=-1,
        )

    assert H is not None, "Could not read Hamiltonian."
    
    if json_args['kpath']['mode'] == 'ConstantEnergy':
        # We have to optimize the path to match approximate a constant energy path
        #
        E0 = json_args['kpath']['energy']
        kmin = json_args['kpath']['kRmin']
        kmax = json_args['kpath']['kRmax']
        etol = json_args['kpath']['etol']

        json_args['kpath']['kwargs']['kR'] = 0.5*(kmax+kmin)
        kpath = KPointLoop(H, **json_args['kpath']['kwargs'])
        radii = kpath.get_radii()

        def get_E(r, ik):
            kpath.update_radius(ik, r)
            Ek = H.eigh(kpath.k[ik])
            idx = np.argmin((Ek-E0)**2)
            # print(kpath.k[ik], idx, Ek[idx], Ek[idx]-E0)
            return Ek[idx]-E0

        for ik in trange(len(kpath)):
            try:
                result = root_scalar(
                        f=get_E, 
                        args=(ik,),
                        bracket=(kmin, kmax),
                )
                radii[ik] = result.root
                # print(ik, radii[ik])
            except ValueError:
                kpath.update_radius(ik, kmin)
                Ekmin = H.eigh(kpath.k[ik])
                idxmin = np.argmin((Ekmin-E0)**2)
                kpath.update_radius(ik, kmax)
                Ekmax = H.eigh(kpath.k[ik])
                idxmax = np.argmin((Ekmax-E0)**2)
                raise ValueError(
                        "    E must lie between E(kRmin) and E(kRmax)\n"
                        +f"    Problem occured for ik={ik}:\n"
                        +f"      For kRmin I find E={Ekmin[idxmin]} (ibnd={idxmin})\n"
                        +f"      For kRmax I find E={Ekmax[idxmax]} (ibnd={idxmax})"
                )

        kpath.update_radii(radii)
        data = kpath.apply.array.renew(zip=False).eigenstate(wrap=wrap, eta=True)
        index_array = np.where(
                np.logical_and(E0-etol < data[0,:,0], data[0,:,0] < E0+etol)
                )
        data = data[:,index_array[0],:]
        nk, nbands, _ = data.shape

    else:
        kpath = get_BZ_obj(H, json_args['kpath']['mode'], **json_args['kpath']['kwargs'])
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
                out.write("\n\n".encode(enc))
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
                out.write("\n\n".encode(enc))

