from pathlib import Path
import zipfile

from pysisyphus.typing import PathLike


def zip_directory(inp_dir: PathLike, zip_fn: PathLike):
    directory = Path(inp_dir).expanduser().resolve()
    with zipfile.ZipFile(zip_fn, mode="w", compression=zipfile.ZIP_DEFLATED) as archive:
        for file_path in directory.rglob("*"):
            archive.write(file_path, arcname=file_path.relative_to(directory))
