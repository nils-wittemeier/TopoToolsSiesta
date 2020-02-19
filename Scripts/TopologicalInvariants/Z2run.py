#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess

import numpy as np
import matplotlib.pyplot as plt
import sisl
import z2pack


def parse_args():
    import argparse
    import json

    class LoadFromFile (argparse.Action):
        def __call__ (self, parser, namespace, values, option_string = None):
            with values as config_file:
                config = json.load(config_file)
                for key in config:
                    setattr(namespace, key, config[key])

    parser = argparse.ArgumentParser(
            description="Calculate Z2 invariant.")

    parser.add_argument(
            "-c", "--config", type=open, action=LoadFromFile,
            help=("json file with options (command line arguments passed"
                  " before this flag will be overwritten, arguments passed"
                  " after this will overwrite options from the file)"))

    parser.add_argument(
            "-i", "--infile", type=str,
            help="input file to calculate spintexture from (fdf-file)")

    parser.add_argument(
            "-s", "--surface", type=str,
            help=("Expression describing the surface on which the WCC / Wilson"
                  " loops should be calculated. The argument should be the"
                  " right side of a lambda expression which parametrizes the"
                  " surface :math:`\\mathbf{k}(t_1, t_2)`, in reduced"
                  " coordinates. It should take two arguments (``float``) and"
                  " return a nested list of ``float`` describing the points in"
                  " k-space. Note that the surface must be closed at least"
                  " along the :math:`t_2` - direction, that is"
                  " :math:`\\mathbf{k}(t_1, 0) = \\mathbf{k}(t_1, 1) +"
                  " \\mathbf{G}`, where :math:`\\mathbf{G}` is an inverse"
                  " lattice vector.  Example: ' [t1 / 2, t2, 0]'"))

    parser.add_argument(
            "--pos_tol", type=float,
            help=("The maximum movement of a WCC for the iteration w.r.t. the"
                  " number of k-points in a single string to converge. The"
                  " iteration can be turned off by setting ``pos_tol=None``."),
            default=1e-2)

    parser.add_argument(
            "--gap-tol", type=float,
            help=("Determines the smallest distance between a gap and its"
                  " neighbouring WCC for the gap check to be satisfied. The"
                  " distance must be larger than ``gap_tol`` times the size of"
                  " the gap. This check is performed only for the largest gap"
                  " in each string of WCC. The check can be turned off by"
                  " setting ``gap_tol=None``."),
            default=0.3)

    parser.add_argument(
            "--move-tol", type=float,
            help=("Determines the largest possible movement between WCC of"
                  " neighbouring strings for the move check to be satisfied."
                  " The movement can be no larger than ``move_tol`` times the"
                  " size of the largest gap between two WCC (from the two"
                  " neighbouring strings, the smaller value is chosen).  The"
                  " check can be turned off by setting ``move_tol=None``."),
            default=0.3)

    parser.add_argument(
            "--num_lines", type=int,
            help="Initial number of strings.",
            default=11)

    parser.add_argument(
            "--min_neighbour_dist", type=float,
            help=("Minimum distance between two strings (no new strings will"
                  " be added, even if the gap check or move check fails)."),
            default=0.01)

    parser.add_argument(
            "--save_file", type=str,
            help="Path to a file where the result should be stored",
            default='./res.json')

    parser.add_argument(
            "--load", type=bool,
            help=("Determines whether the initial result is loaded from"
                  " ``save_file``."),
            default=False)

    parser.add_argument(
            "--iterator", type=int, nargs=3,
            help=("Defines iterator for the number of points in a k-point"
                  "  string. Specify in order: minimum, maximum, step"),
            default=[8, 27, 2])


    parser.add_argument(
            "--bands", type=int, nargs='*',
            help=("Specifies either the number of occupied bands (if it is an"
                  " integer) or which bands should be taken into consideration"
                  " (if it is a list of space separated indices). If no value"
                  " is given, half the given bands are considered."))

    parser.add_argument(
            "--hermitian_tol", type=float,
            help=("Maximum absolute value in the difference between the"
                  " Hamiltonian and its hermitian conjugate. By default the"
                  " test is deactivated entirely."),
            default=None)

    return parser.parse_args()
    


# The Z2Pack tight-binding interface only allows for orthonormal basis sets
# So we have to use the hm (Hamiltonian matrix) interface, to work around
# this limitation.
args = parse_args()
if args.hermitian_tol is None:
    print("Warning: Check whether Hamiltonian is hermitian is deactivated. To"
          " activate the test please specify \'hermitian_tol\' option.")

sile = sisl.get_sile(args.infile)
geom = sile.read_geometry()
H = sile.read_hamiltonian(geometry=geom)

pos = H.fxyz[list(map(H.geom.o2a, np.arange(H.no)))]
if H.spin.spins >= 4:
    # In non-collinear systems we need to double the amount of orbitals
    pos = np.repeat(pos, 2, axis=0)

# print("hermitian_tol", type(args.hermitian_tol))
# print("bands", type(args.bands))
# print("surface", type(args.surface))
# print("save_file", type(args.save_file))
# print("num_lines", type(args.num_lines))
# print("iterator", type(args.iterator))
# print("min_neighbour_dist", type(args.min_neighbour_dist))
# print("load", type(args.load))
# assert False
system = z2pack.hm.System(
        hamilton=lambda k: H.Hk(k=k, dtype=np.complex64, format='array'),
        hermitian_tol=args.hermitian_tol,
        overlap=lambda k: H.Sk(k=k, dtype=np.complex64, format='array'),
        dim=3,
        pos=pos,
        bands=args.bands,
)

# Run the WCC calculations
result = z2pack.surface.run(
    system=system,
    surface=lambda t1, t2: eval(args.surface),
    save_file=args.save_file,
    num_lines=args.num_lines,
    iterator=range(*args.iterator),
    min_neighbour_dist=args.min_neighbour_dist,
    load=args.load,
)

fig, ax = plt.subplots(1, 1, figsize=(9,5))
z2pack.plot.wcc(result, axis=ax)
plt.savefig('WCC.svg', bbox_inches='tight')

print(
    'Z2 topological invariant: {0}'.format(
        z2pack.invariant.z2(result)
    )
)

# Run the WCC calculations
result = z2pack.surface.run(
    system=system,
    surface=z2pack.shape.Sphere(center=(0.,0.,0), radius=0.001),
    save_file=args.save_file,
    num_lines=args.num_lines*3,
    iterator=range(*args.iterator),
    min_neighbour_dist=args.min_neighbour_dist,
    load=args.load,
)

print(
    'Chern Number: {0}'.format(
        z2pack.invariant.chern(result)
    )
)

