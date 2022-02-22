import logging
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile

from natsort import natsorted

from pysisyphus.config import get_cmd, OUT_DIR_DEFAULT
from pysisyphus.constants import BOHR2ANG
from pysisyphus import logger
from pysisyphus import helpers_pure


class Calculator:

    conf_key = None

    def __init__(
        self,
        calc_number=0,
        charge=0,
        mult=1,
        base_name="calculator",
        pal=1,
        mem=1000,
        check_mem=True,
        retry_calc=1,
        last_calc_cycle=None,
        clean_after=True,
        out_dir=OUT_DIR_DEFAULT,
    ):
        """Base-class of all calculators.

        Meant to be extended.

        Parameters
        ----------
        calc_number : int, default=0
            Identifier of the Calculator. Used in distinguishing it from
            other Calculators, e.g. in ChainOfStates calculations. Also
            used in the creation of filenames.
        charge : int, default=0
            Molecular charge.
        mult : int, default=1
            Molecular multiplicity (1 = singlet, 2 = doublet, ...)
        base_name : str, default=calculator
            Generated filenames will start with this string.
        pal : int, default=1
            Positive integer that gives the number of physical cores to
            use on 1 node.
        mem : int, default=1000
            Mememory per core in MB. The total amount of memory is given as
            mem*pal.
        check_mem : bool, default=True
            Whether to adjust the requested memory if too much is requested.
        retry_calc : int, default=0
            Number of additional retries when calculation failed.
        last_calc_cycle : int
            Internal variable used in restarts.
        clean_after : bool
            Delete temporary directory were calculations were executed
            after a calculation.
        out_dir : str
            Path that is prepended to generated filenames.
        """

        self.logger = logging.getLogger("calculator")

        self.calc_number = calc_number
        self.charge = int(charge)
        self.mult = int(mult)
        self.base_name = base_name
        self.pal = int(pal)
        assert self.pal > 0, "pal must be a non-negative integer!"
        if check_mem:
            mem = helpers_pure.check_mem(int(mem), pal, logger=self.logger)
        self.mem = mem
        # Disasble retries if check_termination method is not implemented
        self.retry_calc = int(retry_calc) if hasattr(self, "check_termination") else 0
        assert self.retry_calc >= 0
        try:
            self.out_dir = Path(out_dir).resolve()
        except TypeError:
            self.out_dir = Path(OUT_DIR_DEFAULT).resolve()
        self.out_dir.mkdir(parents=True, exist_ok=True)

        # Extensions of the files to keep after running a calculation.
        # Usually overridden in derived classes.
        self.to_keep = ()
        # How many calculations were already run
        self.calc_counter = 0
        # Handle restarts
        if last_calc_cycle:
            self.calc_counter = int(last_calc_cycle) + 1
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
        self.backup_dir = None

    def get_cmd(self, key="cmd"):
        assert self.conf_key, "'conf_key'-attribute is missing for this calculator!"

        try:
            return get_cmd(section=self.conf_key, key=key, use_defaults=True)
        except KeyError:
            logger.debug(f"Failed to load key '{key}' from section '{self.conf_key}'!")

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

        self.logger.debug(f"{self.name}, cycle {self.calc_counter:03d}: {message}")

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

        if counter is None:
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
            coord_line = (fmt + fmt + fmt).format(*coord) + atom.lower() + "\n"
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
            [
                "{} {:10.08f} {:10.08f} {:10.08f}".format(a, *c)
                for a, c in zip(atoms, coords)
            ]
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

    def run(
        self,
        inp,
        calc,
        add_args=None,
        env=None,
        shell=False,
        hold=False,
        keep=True,
        cmd=None,
        inc_counter=True,
        run_after=True,
        parser_kwargs=None,
        symlink=True,
    ):
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
        cmd : str or iterable, optional
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

        self.backup_dir = None
        path = self.prepare(inp)
        self.log(f"Running in {path} on {platform.node()}")
        if cmd is None:
            cmd = self.base_cmd

        if isinstance(cmd, str):
            cmd = [cmd]

        args = cmd + [self.inp_fn]
        if add_args:
            args.extend(add_args)
        if not env:
            env = os.environ.copy()
        tmp_out_fn = path / self.out_fn
        with open(tmp_out_fn, "w") as handle:
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

                try:
                    os.symlink(tmp_out_fn, sym_fn)
                    self.log(f"Created symlink in '{sym_fn}'")
                # This may happen if we use a dask scheduler
                except FileExistsError:
                    self.log("Symlink already exists. Skipping generation.")

            # Do at least one cycle. When retries are disabled retry_calc == 0
            # and range(0+1) will result in one cycle
            added_retry_args = False
            for retry in range(self.retry_calc + 1):
                result = subprocess.Popen(
                    args,
                    cwd=path,
                    stdout=handle,
                    stderr=subprocess.PIPE,
                    env=env,
                    shell=shell,
                )
                result.wait()
                try:
                    normal_termination = False
                    # Calling check_termination may result in an exception and
                    # normal_termination will stay at False
                    normal_termination = self.check_termination(tmp_out_fn)
                # Method check_termination may not be implemented, so we will always
                # do only one try.
                except AttributeError:
                    normal_termination = True
                # The out file may not be present
                except FileNotFoundError:
                    self.log(
                        f"Could not find out-file {str(tmp_out_fn)} for termination status check!"
                    )

                if normal_termination:
                    break
                else:
                    print("Abnormal termination! Retrying calculation.")
                    shutil.copy(tmp_out_fn, str(tmp_out_fn) + f".fail_{retry:02d}")
                    try:
                        self.clean_tmp(path)
                    except AttributeError:
                        self.log(f"'self.clean_path()' not implemented!")
                    # Clear tmp_out_fn
                    handle.seek(0)
                    handle.truncate()
                    self.log("Detected abnormal termination! Retrying calculation.")

                    if not added_retry_args:
                        try:
                            args += self.get_retry_args()
                            added_retry_args = True
                        except AttributeError:
                            self.log(f"'self.get_retry_args()' not implemented!")

        # Parse results for desired quantities
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
            self.backup_dir = backup_dir
            if backup_dir.exists():
                shutil.rmtree(backup_dir)
            shutil.copytree(path, backup_dir)
            print(
                f"Copied contents of\n\t'{path}'\nto\n\t'{backup_dir}'.\n"
                "Consider checking the log files there.\n"
            )
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
                assert len(globbed) <= 1, (
                    f"Expected at most one file "
                    f"matching {pattern} in {path}. Found {len(globbed)} "
                    f"files instead ({', '.join([g.name for g in globbed])})!"
                )
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

    def get_restart_info(self):
        """Return a dict containing chkfiles.

        Returns
        -------
        restart_info : dict
            Dictionary holding the calculator state. Used for restoring calculaters
            in restarted calculations.
        """
        try:
            # Convert possible Paths to str
            chkfiles = {k: str(v) for k, v in self.get_chkfiles().items()}
        except AttributeError:
            chkfiles = dict()

        restart_info = {
            "base_name": self.base_name,
            "calc_number": self.calc_number,
            "calc_counter": self.calc_counter,
            "chkfiles": chkfiles,
        }

        return restart_info

    def verify_chkfiles(self, chkfiles):
        """Checks if given chkfiles exist and return them as Paths

        Parameters
        ----------
        chkfiles : dict
            Dictionary holding the chkfiles. The keys correspond to the attribute
            names, the values are strs holding the (potentially full) filename (path).

        Returns
        -------
        paths : dict
            Dictionary of Paths.
        """
        paths = {}
        for key, chkfile in chkfiles.items():
            chkfile = Path(chkfile)
            # If the chkfile exists at the given path we use it as it is.
            if not chkfile.exists():
                self.log(
                    f"Given chkfile '{chkfile}' could not be found! Dropping "
                    "absolute part and trying only its name."
                )
                # Check if relative path exists. This may happen if the calculation
                # has been moved to a different folder.
                name = Path(chkfile.name)
                if name.exists():
                    chkfile = name
                else:
                    self.log(f"'{name}' could not be found! Skipping this chkfile.")
                    continue
            paths[key] = chkfile
        return paths

    def set_restart_info(self, restart_info):
        """Sets restart information (chkfiles etc.) on the calculator.

        Parameters
        -------
        restart_info : dict
            Dictionary holding the calculator state. Used for restoring calculaters
            in restarted calculations.
        """
        try:
            chkfiles = self.verify_chkfiles(restart_info.pop("chkfiles"))
            self.set_chkfiles(chkfiles)
        except KeyError:
            self.log("No chkfiles preset in restart_info")
        except AttributeError:
            self.log(
                "Found chkfiles on restart_info, but 'set_chkfiles' is not "
                "implemented for Calculator."
            )

        self.log("Setting restart_info")
        for key, value in restart_info.items():
            setattr(self, key, value)
            self.log(f"\t{key}: {value}")

    def print_capabilities(self):
        print(
            f"    Can retry?: {hasattr(self, 'check_termination')}\n"
            f"Can track ES??: {hasattr(self, 'prepare_overlap_data')}\n"
        )
