#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess

import numpy as np
import matplotlib.pyplot as plt
import sisl
import z2pack

sile = sisl.get_sile('work/Bi2D_BHex.fdf')
geom = sile.read_geometry()
H = sile.read_hamiltonian(geometry=geom)

pos = H.fxyz[list(map(H.geom.o2a, np.arange(H.no)))]
if H.spin.spins >= 4:
    # In non-collinear systems we need to double the amount of orbitals
    pos = np.repeat(pos, 2, axis=0)

# Treat the highest occupied bands separately from the rest
system_full = z2pack.hm.System(
        hamilton=lambda k: H.Hk(k=k, dtype=np.complex64, format='array'),
        hermitian_tol=1e-4,
        overlap=lambda k: H.Sk(k=k, dtype=np.complex64, format='array'),
        dim=3,
        pos=pos,
        #bands=30,
        bands=range(24,30),
)

system_hocc = z2pack.hm.System(
        hamilton=lambda k: H.Hk(k=k, dtype=np.complex64, format='array'),
        hermitian_tol=1e-4,
        overlap=lambda k: H.Sk(k=k, dtype=np.complex64, format='array'),
        dim=3,
        pos=pos,
        bands=[28,29],
)

system_occ = z2pack.hm.System(
        hamilton=lambda k: H.Hk(k=k, dtype=np.complex64, format='array'),
        hermitian_tol=1e-4,
        overlap=lambda k: H.Sk(k=k, dtype=np.complex64, format='array'),
        dim=3,
        pos=pos,
        bands=range(0,28),
)

# Run the WCC calculations
result_full = z2pack.surface.run(
    system=system_full,
    surface=lambda t1, t2: [t1/2, t2, 0],
    save_file='./res_full.json',
    iterator=range(8,40,2),
    pos_tol=1e-2,
    gap_tol=2e-2,
    move_tol=0.3,
    load=False,
)

print(
    'Z2 topological invariant of the full system: {0}'.format(
#        (z2pack.invariant.z2(result) + z2pack.invariant.z2(result2)) % 2
        z2pack.invariant.z2(result_full)
    )
)

fig, ax = plt.subplots(1, 2, figsize=(9,5), sharex=True, sharey=True)
z2pack.plot.wcc(result_full, axis=ax[0])
ax[0].set_xlabel(r'$\bar{x} [a_x]$')
ax[0].set_xlabel(r'$k_y [\pi/a_y]$')
ax[1].set_xlabel(r'$k_y [\pi/a_y]$')
handles, labels = ax[0].get_legend_handles_labels()
handles, labels = ax[1].get_legend_handles_labels()
plt.savefig('WCC.svg')

assert False

result_hocc = z2pack.surface.run(
    system=system_hocc,
    surface=lambda t1, t2: [t1/2, t2, 0],
    save_file='./res_hocc.json',
    load=True,
)

result_lower = z2pack.surface.run(
    system=system_lower,
    surface=lambda t1, t2: [t1/2, t2, 0],
    save_file='./res_lower.json',
    load=True,
)

fig, ax = plt.subplots(1, 2, figsize=(9,5), sharex=True, sharey=True)
z2pack.plot.wcc(result_full, axis=ax[0])
z2pack.plot.wcc(result_lower, axis=ax[1])
ax[0].set_xlabel(r'$\bar{x} [a_x]$')
ax[0].set_xlabel(r'$k_y [\pi/a_y]$')
ax[1].set_xlabel(r'$k_y [\pi/a_y]$')
handles, labels = ax[0].get_legend_handles_labels()
handles, labels = ax[1].get_legend_handles_labels()
plt.savefig('WCC.svg')

print(
    'Z2 topological invariant for bands 1-28: {0}'.format(
        z2pack.invariant.z2(result_lower)
    )
)

print(
    'Z2 topological invariant for bands 29, 30: {0}'.format(
        z2pack.invariant.z2(result_hocc)
    )
)

print(
    'Z2 topological invariant of the full system: {0}'.format(
#        (z2pack.invariant.z2(result) + z2pack.invariant.z2(result2)) % 2
        z2pack.invariant.z2(result_full)
    )
)
