Relaxed Scans
*************

Pysisyphus currently supports 1d-relaxed scans, where one internal
coordinate is kept frozen, while the remaining ones are relaxed. A sample
input is given below.

.. code:: yaml

    geom:
     type: redund
     fn: lib:h2o2_hf_321g_opt.xyz
    calc:
     type: pyscf
     basis: 321g
     pal: 4
    scan:
     type: BOND
     indices: [2, 3]
     # start is optional; if not specified the initial value of the primitive
     # internal coordinate is used.
     start: 3.0
     end: 5.0
     # steps is optional and can be omitted, but then step_size must be given
     steps: 5
     # step_size: 0.4
     # Enable symmetric scan in both directions, starting from the initial coordinates
     # symmetric: False
     # 
     # By specifying an 'opt' block, the default optimization values, e.g.
     # convergence thresholds can be altered.
     # opt:
      # thresh: gau

If no `start` value is specified in the `scan:` block,
then the initial value is used. Specification of `steps` is mandatory, as it determines
the number of optimizations to carry out.
If not directly given (`step_size`), a suitable step size is derived from
an `end` value. Either `end` or `step_size` must be given.

The values `step_size`, `start`, `end` *have to be provided in atomic units or radians*.
When the appropriate YAML constructor is used, inputs can still be in Ångström or degrees.
The constructor `!angstrom` converts input in Ångström to Bohr and `!deg` converts input
from degrees to radians. See below for an example.

.. code:: yaml

    [...]
    scan:
     type: DIHEDRAL
     indices: [0, 1, 2, 3]
     steps: 5
     step_size: !deg 10.
    [...]

In the example above, a relaxed scan around a torsion is carried out. In total,
5 steps of 10 degree each will be taken. By using the YAML constructor `!deg`, the input
can be given in degrees, instead of radians. Internally, the degrees input is then
converted to the appropriate unit, in this case radians.

Naive relaxed scan around the equilibrium value of a primitive internal may be affected
by hysteresis, if the system under study is sufficiently complicated.

Relaxed scans are always carried out in redundant internal coordinates. Please
see :ref:`Optimization of Minima` for the available optimization options. Supported
primitive internals (`type`) are found in: :ref:`Types of Primitive Coordinates`.

All geometries, including the initial one are written to `relaxed_scan.trj` in the
current working directory.
