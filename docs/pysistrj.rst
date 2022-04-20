Working with Geometries
***********************
**TBD**

Pleas see `pysistrj --help` for a list of possible invocations.

.. code:: bash


    usage: Utility to transform .xyz and .trj files. [-h]
                                                     (--between BETWEEN | --align | --split | --reverse | --cleantrj | --spline | --first FIRST | --every EVERY | --center | --centerm | --translate TRANSLATE TRANSLATE TRANSLATE | --append | --join | --match | --std | --shake | --internals | --get GET | --origin | --fragsort)
                                                     [--scale SCALE] [--seed SEED]
                                                     [--idpp | --lst | --redund]
                                                     [--noipalign] [--bohr]
                                                     [--noxyz]
                                                     [--atoms ATOMS [ATOMS ...]]
                                                     [--add_prims ADD_PRIMS]
                                                     fns [fns ...]

    positional arguments:
      fns                   Filenames of .xyz and/or .trj files (xyz and trj can
                            be mixed).

    optional arguments:
      -h, --help            show this help message and exit
      --between BETWEEN     Interpolate additional images.
      --align               Align geometries onto the first geometry.
      --split               Split a supplied geometries in multiple .xyz files.
      --reverse             Reverse a .trj file.
      --cleantrj            Keep only the first four columns of xyz/trj files.
      --spline              Evenly redistribute geometries along a splined path.
      --first FIRST         Copy the first N geometries to a new .trj file.
      --every EVERY         Create new .trj with every N-th geometry. Always
                            includes the first and last point.
      --center              Move the molecules centroid into the origin.
      --centerm             Move the molecules center of mass into the origin.
      --translate TRANSLATE TRANSLATE TRANSLATE
                            Translate the molecule by the given vector given in
                            Ångström.
      --append              Combine the given .xyz files into one .xyz file.
      --join                Combine the given .xyz/.trj files into one .trj file.
      --match               Resort the second .xyz file so the atom order matches
                            the first .xyz file. Uses the hungarian method.
      --std                 Move supplied geometry to its standard orientation.
      --shake               Shake (randomly displace) coordiantes.
      --internals           Print automatically generated internals.
      --get GET             Get n-th geometry. Expects 0-based index input.
      --origin              Translate geometry, so that min(X/Y/Z) == 0.
      --fragsort            Resort atoms by fragments.
      --idpp                Interpolate using Image Dependent Pair Potential.
      --lst                 Interpolate by linear synchronous transit.
      --redund              Interpolate in internal coordinates.
      --noipalign           Don't align geometries when interpolating.
      --bohr                Input geometries are in Bohr instead of Angstrom.
      --noxyz               Disable dumping of single .xyz files.
      --atoms ATOMS [ATOMS ...]
                            Used with --internals. Only print primitives including
                            the given atoms.
      --add_prims ADD_PRIMS
                            Used with --internals. Define additional primitives.
                            Expects a string representation of a nested list that
                            can be parsed as YAML e.g.
                            [[10,30],[1,2,3],[4,5,6,7]].

      --scale SCALE         Scales the displacement in --shake.
      --seed SEED           Initialize the RNG for reproducible results.
