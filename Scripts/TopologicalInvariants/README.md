# Topological invariants:
The script accepts various options either via a .json file or via command line
arguments.

To calculate topological invariants from siesta outputs either the HSX or
TSHS-file is required. The TSHS-file is sufficient on its own. When using the
HSX-file more output files are required.

Required files:
 - TSHS or
 - HSX + fdf + ion(.nc)

The 'infile' argument can be either the fdf-file or the TSHS-file.

## Example:
    python ../Z2run.py b-hex.TSHS "[t1/2,t2,0]" --hermitian-tol=1e-4 --min-neighbour-dist=1e-10

## Usage:
    usage: Z2run.py [-h] [-i INFILE] [-s SURFACE] [--pos_tol POS_TOL]
                    [--gap-tol GAP_TOL] [--move-tol MOVE_TOL]
                    [--num_lines NUM_LINES]
                    [--min_neighbour_dist MIN_NEIGHBOUR_DIST]
                    [--save_file SAVE_FILE] [--load LOAD]
                    [--iterator ITERATOR ITERATOR ITERATOR]
                    [--bands [BANDS [BANDS ...]]] [--hermitian_tol HERMITIAN_TOL]
                    config

    Calculate Z2 invariant.

    positional arguments:
      config                json file with options (command line arguments passed
                            before this flag will be overwritten, arguments passed
                            after this will overwrite options from the file)

    optional arguments:
      -h, --help            show this help message and exit
      -i INFILE, --infile INFILE
                            input file to calculate spintexture from (fdf-file)
      -s SURFACE, --surface SURFACE
                            Expression describing the surface on which the WCC /
                            Wilson loops should be calculated. The argument should
                            be the right side of a lambda expression which
                            parametrizes the surface :math:`\mathbf{k}(t_1, t_2)`,
                            in reduced coordinates. It should take two arguments
                            (``float``) and return a nested list of ``float``
                            describing the points in k-space. Note that the
                            surface must be closed at least along the :math:`t_2`
                            - direction, that is :math:`\mathbf{k}(t_1, 0) =
                            \mathbf{k}(t_1, 1) + \mathbf{G}`, where
                            :math:`\mathbf{G}` is an inverse lattice vector.
                            Example: ' [t1 / 2, t2, 0]'
      --pos_tol POS_TOL     The maximum movement of a WCC for the iteration w.r.t.
                            the number of k-points in a single string to converge.
                            The iteration can be turned off by setting
                            ``pos_tol=None``.
      --gap-tol GAP_TOL     Determines the smallest distance between a gap and its
                            neighbouring WCC for the gap check to be satisfied.
                            The distance must be larger than ``gap_tol`` times the
                            size of the gap. This check is performed only for the
                            largest gap in each string of WCC. The check can be
                            turned off by setting ``gap_tol=None``.
      --move-tol MOVE_TOL   Determines the largest possible movement between WCC
                            of neighbouring strings for the move check to be
                            satisfied. The movement can be no larger than
                            ``move_tol`` times the size of the largest gap between
                            two WCC (from the two neighbouring strings, the
                            smaller value is chosen). The check can be turned off
                            by setting ``move_tol=None``.
      --num_lines NUM_LINES
                            Initial number of strings.
      --min_neighbour_dist MIN_NEIGHBOUR_DIST
                            Minimum distance between two strings (no new strings
                            will be added, even if the gap check or move check
                            fails).
      --save_file SAVE_FILE
                            Path to a file where the result should be stored
      --load LOAD           Determines whether the initial result is loaded from
                            ``save_file``.
      --iterator ITERATOR ITERATOR ITERATOR
                            Defines iterator for the number of points in a k-point
                            string. Specify in order: minimum, maximum, step
      --bands [BANDS [BANDS ...]]
                            Specifies either the number of occupied bands (if it
                            is an integer) or which bands should be taken into
                            consideration (if it is a list of space separated
                            indices). If no value is given, half the given bands
                            are considered.
      --hermitian_tol HERMITIAN_TOL
                            Maximum absolute value in the difference between the
                            Hamiltonian and its hermitian conjugate. By default
                            the test is deactivated entirely.
