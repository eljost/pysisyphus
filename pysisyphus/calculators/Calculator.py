#!/usr/bin/env python3

import logging
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile

from natsort import natsorted

from pysisyphus.config import Config
from pysisyphus.constants import BOHR2ANG


class Calculator:

    conf_key = None

    def __init__(self, calc_number=0, charge=0, mult=1,
                 base_name="calculator", pal=1,
                 last_calc_cycle=None, clean_after=True, out_dir="./"):
        """Base-class of all calculators.

        Meant to be extended.

        Parameters
        ----------
        calc_number : int
            Identifier of the Calculator. Used in distinguishing it from
            other Calculators, e.g. in ChainOfStates calculations. Also
            used in the creation of filenames.
        charge : int
            Molecular charge.
        mult : int
            Molecular multiplicity (1 = singlet, 2 = doublet, ...)
        base_name : str
            Generated filenames will start with this string.
        pal : int
            Positive integer that gives the number of physical cores to
            use on 1 node.
        last_calc_cycle : int
            Internal variable used in restarts.
        clean_after : bool
            Delete the temporary directory after calculations.
        out_dir : str
            Path that is prepended to generated filenames.
        """
        self.calc_number = calc_number
        self.charge = int(charge)
        self.mult = int(mult)
        self.base_name = base_name
        self.pal = int(pal)
        self.out_dir = Path(out_dir).resolve()

        assert pal > 0, "pal must be a non-negative integer!"

        # Extensions of the files to keep after running a calculation.
        # Usually overridden in derived classes.
        self.to_keep = ()
        # How many calculations were already run
        self.calc_counter = 0
        # Handle restarts
        if last_calc_cycle:
            self.calc_counter = int(last_calc_cycle)+1
            self.reattach(int(last_calc_cycle))
            self.log(f"Set {self.calc_counter} for this calculation")
        self.clean_after = clean_after

        self.inp_fn = "calc.inp"
        self.out_fn = "calc.out"
        # When this is set the run() method will use this path
        # instead of creating a new one.
        # Currently this is only used with the Turbomole calculator.
        self.path_already_prepared = None
        self.last_run_path = None

    def get_cmd(self, key):
        assert self.conf_key, \
            "Tried loading a cmd from the config file but no conf_key " \
            "was specified in the Calculator!"

        try:
            return Config[self.conf_key][key]
        except KeyError:
            print(f"Failed to load key '{key}' from section '{self.conf_key}'. "
                   "Exiting!")
            sys.exit()

    @property
    def name(self):
        return f"{self.base_name}_{self.calc_number:03d}"

    def log(self, message=""):
        """Write a log message.

        Wraps the logger variable.

        Parameters
        ----------
        message : str
            Message to be logged.
        """

        logger = logging.getLogger("calculator")
        logger.debug(f"{self.name}_cyc_{self.calc_counter:03d}, {message}")

    def get_energy(self, atoms, coords):
        """Meant to be extended."""
        raise Exception("Not implemented!")

    def get_forces(self, atoms, coords):
        """Meant to be extended."""
        raise Exception("Not implemented!")

    def get_hessian(self, atoms, coords):
        """Meant to be extended."""
        raise Exception("Not implemented!")

    def make_fn(self, name, counter=None, return_str=False):
        """Make a full filename.

        Return a full filename including the calculator name and the
        current counter given a suffix.

        Parameters
        ----------
        name: str
            Suffix of the filename.
        counter : int, optional
            If not given use the current calc_counter.
        return_str : int, optional
            Return a string instead of a Path when True.

        Returns
        -------
        fn : str
            Filename.
        """            

        if not counter:
            counter = self.calc_counter
        fn = self.out_dir / f"{self.name}.{counter:03d}.{name}"
        if return_str:
            fn = str(fn)
        return fn

    def prepare_path(self, use_in_run=False):
        """Get a temporary directory handle.

        Create a temporary directory that can later be used in a calculation.

        Parameters
        ----------
            use_in_run : bool, option
                Sets the internal variable ``self.path_already_prepared`` that
                is later read by ``self.run()``. No new temporary directory will
                be created in ``self.run()``.

        Returns
        -------
            path: Path
                Prepared directory.
        """
        
        prefix = f"{self.name}_{self.calc_counter:03d}_"
        path = Path(tempfile.mkdtemp(prefix=prefix))
        if use_in_run:
            self.path_already_prepared = path
        return path

    def prepare(self, inp):
        """Prepare a temporary directory and write input.

        Similar to prepare_path, but the input is also written into
        the prepared directory.

        Paramters
        ---------
        inp : str
            Input to be written into the file ``self.inp_fn`` in
            the prepared directory.

        Returns
        -------
            path: Path
                Prepared directory.
        """

        if not self.path_already_prepared:
            path = self.prepare_path()
        else:
            path = self.path_already_prepared

        # Calculators like Turbomole got no input.
        if inp:
            inp_path = path / self.inp_fn
            with open(inp_path, "w") as handle:
                handle.write(inp)

        return path

    def prepare_input(self, atoms, coords, calc_type):
        """Meant to be extended."""
        raise Exception("Not implemented!")

    def print_out_fn(self, path):
        """Print calculation output.

        Prints the output of a calculator after a calculation.

        Parameters
        ----------
        path : Path
            Temporary directory of the calculation.
        """
        with open(path / self.out_fn) as handle:
            text = handle.read()
        print(text)

    def prepare_turbo_coords(self, atoms, coords):
        """Get a Turbomole coords string.

        Parameters
        ----------
        atoms : iterable
            Atom descriptors (element symbols).
        coords: np.array, 1d
            1D-array holding coordinates in Bohr.

        Returns
        -------
        coords: str
            String holding coordinates in Turbomole coords format.
        """
        fmt = "{:<20.014f}"
        coord_str = "$coord\n"
        for atom, coord in zip(atoms, coords.reshape(-1, 3)):
            coord_line = (fmt+fmt+fmt).format(*coord) + atom.lower() + "\n"
            coord_str += coord_line
        coord_str += "$end"
        return coord_str

    def prepare_coords(self, atoms, coords):
        """Get 3d coords in Angstrom.

        Reshape internal 1d coords to 3d and convert to Angstrom.

        Parameters
        ----------
        atoms : iterable
            Atom descriptors (element symbols).
        coords: np.array, 1d
            1D-array holding coordinates in Bohr.

        Returns
        -------
        coords: np.array, 3d
            3D-array holding coordinates in Angstrom.
        """
        coords = coords.reshape(-1, 3) * BOHR2ANG
        coords = "\n".join(
                ["{} {:10.08f} {:10.08f} {:10.08f}".format(a, *c) for a, c in zip(atoms, coords)]
        )
        return coords

    def prepare_xyz_string(self, atoms, coords):
        """Returns a xyz string in Angstrom.

        Parameters
        ----------
        atoms : iterable
            Atom descriptors (element symbols).
        coords: np.array, 1d
            1D-array holding coordinates in Bohr.

        Returns
        -------
        xyz_str: string
            Coordinates in .xyz format.
        """

        return f"{len(atoms)}\n\n{self.prepare_coords(atoms, coords)}"

    def run(self, inp, calc, add_args=None, env=None, shell=False,
            hold=False, keep=True, cmd=None, inc_counter=True,
            run_after=True, parser_kwargs=None, symlink=True):
        """Run a calculation.

        The bread-and-butter method to actually run an external quantum
        chemistry code.

        Parameters
        ----------
        inp : str
            Input for the external program that is written to the temp-dir.
        calc : str, hashable
            Key (and more or less type of calculation) to select the right
            parsing function from ``self.parser_funcs``.
        add_args : iterable, optional
            Additional arguments that will be appended to the program call.
        env : Environment, optional
            A potentially modified environment for the subprocess call.
        shell : bool, optional
            Use a shell to execute the program call. Need for Turbomole were
            we chain program calls like dscf; escf.
        hold : bool, optional
            Wether to remove the temporary directory after the calculation.
        keep : bool, optional
            Wether to backup files as specified in ``self.to_keep()``. Usually
            you want this.
        cmd : str, optional
            Overwrites ``self.base_cmd``.
        inc_counter : bool, optional
            Wether to increment the counter after a calculation.

        Returns
        -------
        results : dict
            Dictionary holding all applicable results of the calculations
            like the energy, a forces vector and/or excited state energies
            from TDDFT.
        """

        path = self.prepare(inp)
        self.log(f"Running in {path} on {platform.node()}")
        if cmd:
            args = [cmd, self.inp_fn]
        else:
            args = [self.base_cmd, self.inp_fn]
        if add_args:
            args.extend(add_args)
        if not env:
            env = os.environ.copy()
        with open(path / self.out_fn, "w") as handle:
            if symlink:
                # We can't use resolve here as a previous symlink may already
                # exist. Calling resolve would translate this to the original
                # out file in some tempdir (that is already deleted ...).
                # sym_fn = Path("cur_out").resolve()
                sym_fn = self.out_dir / "cur_out"
                try:
                    os.remove(sym_fn)
                except FileNotFoundError:
                    pass
                os.symlink(path / self.out_fn, sym_fn)
                self.log(f"Created symlink in '{sym_fn}'")
            result = subprocess.Popen(args, cwd=path,
                                      stdout=handle, stderr=subprocess.PIPE,
                                      env=env, shell=shell)
            result.wait()
        try:
            if run_after:
                self.run_after(path)
            parser_kwargs = {} if parser_kwargs is None else parser_kwargs
            results = self.parser_funcs[calc](path, **parser_kwargs)
            if keep:
                self.keep(path)
        except Exception as err:
            print("Crashed input:")
            print(inp)
            backup_dir = Path(os.getcwd()) / f"crashed_{self.name}"
            if backup_dir.exists():
                shutil.rmtree(backup_dir)
            shutil.copytree(path, backup_dir)
            raise err
        finally:
            if (not hold) and self.clean_after:
                self.clean(path)
                if inc_counter:
                    self.calc_counter += 1

        self.path_already_prepared = None
        self.last_run_path = path
        return results

    def run_after(self, path):
        """Meant to be extended.

        This method is called after a calculation was done, but before
        entering ``self.keep()`` and ``self.clean()``. Can be used to call
        tools like formchk or ricctools.
        """

    def prepare_pattern(self, raw_pat):
        """Prepare globs.

        Transforms an entry of ``self.to_keep`` into a glob and a key
        suitable for the use in ``self.keep()``.

        Parameters
        ----------
        raw_pat : str
            Entry of ``self.to_keep``

        Returns
        -------
        pattern : str
            Glob that can be used in Path.glob()
        multi : bool
            Flag if glob may match multiple files.
        key : str
            A key to be used in the ``kept_fns`` dict.
        """
        key_given = None
        if ":" in raw_pat:
            key_given, raw_pat = raw_pat.split(":")
        # Indicates if multiple files are expected
        multi = "*" in raw_pat
        # Drop '*' as it just indicates if we expect multiple matches
        raw_pat = raw_pat.replace("*", "")
        # Interpret it as prefix and drop the two underscores
        if raw_pat.startswith("__"):
            pattern = f"{raw_pat[2:]}*"
            pattern_key = f"{raw_pat[2:]}s"
        # Use raw_pat as suffix
        else:
            pattern = f"*{raw_pat}"
            pattern_key = f"{raw_pat}"
        if key_given:
            pattern_key = key_given
        pattern_key = pattern_key.lower()
        return pattern, multi, pattern_key

    def keep(self, path):
        """Backup calculation results.

        Parameters
        ----------
        path : Path
            Temporary directory of the calculation.

        Returns
        -------
        kept_fns : dict
            Dictonary holding the filenames that were backed up. The keys
            correspond to the type of file.
        """

        kept_fns = dict()
        for raw_pattern in self.to_keep:
            pattern, multi, key = self.prepare_pattern(raw_pattern)
            globbed = natsorted(path.glob(pattern))
            if not multi:
                assert(len(globbed) <= 1), f"Expected at most one file " \
                 f"matching {pattern} in {path}. Found {len(globbed)} " \
                 f"files instead ({', '.join([g.name for g in globbed])})!"
            else:
                kept_fns[key] = list()
            for tmp_fn in globbed:
                base = tmp_fn.name
                new_fn = self.make_fn(base)
                shutil.copy(tmp_fn, new_fn)
                if multi:
                    kept_fns[key].append(new_fn)
                else:
                    kept_fns[key] = new_fn
        return kept_fns

    def clean(self, path):
        """Delete the temporary directory.

        Parameters
        ----------
        path : Path
            Directory to delete.
        """
        shutil.rmtree(path)
        self.log(f"Cleaned {path}")

    def reattach(self, last_calc_cycle):
        """Meant to be extended.
        
        When restarting the calculator set all attributes to restore the
        previous state."""
        self.calc_counter = last_calc_cycle
