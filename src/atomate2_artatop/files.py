"""Module defining functions for manipulating artatop files."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from atomate2.common.files import copy_files, get_zfile, gunzip_files
from atomate2.utils.file_client import FileClient, auto_fileclient
from atomate2.utils.path import strip_hostname

if TYPE_CHECKING:
    from pathlib import Path


ARTATOP_OUTPUT_FILES = [
    "re_lin",
    "re_nlin",
    "re_art",
    "result.art_IND",
    "result.art_TOT",
    "result.re",
    "energy-band.dat",
]

ARTATOP_INPUT_FILES = [
    "input_lin",
    "input_nlin",
    "input_art",
]

VASP_OUTPUT_FILES = [
    "INCAR",
    "POSCAR",
    "OUTCAR",
    "POTCAR",
    "vasprun.xml",
    "PROCAR",
    "CONTCAR",
    "OPTIC",
    "WAVEDER",
    "CHG",
    "CHGCAR",
    "DOSCAR",
    "EIGENVAL",
    "IBZKPT",
    "OSZICAR",
    "WAVECAR",
    "XDATCAR",
]

ARTATOP_OUTPUT_FOLDERS = ["out_lin", "out_nonlin"]

logger = logging.getLogger(__name__)


@auto_fileclient
def copy_artatop_files(
    src_dir: Path | str,
    dest_dir: Path | str,
    src_host: str | None = None,
    file_client: FileClient = None,
    force: bool = False,
) -> None:
    """
    Copy ARTATOP files from the source directory to the destination directory.

    Parameters
    ----------
    src_dir : Path or str
        The source directory containing the input files.
    dest_dir : Path or str
        The destination directory where the files will be copied.
    src_host : str or None
        The source hostname for remote filesystems.
    file_client : FileClient
        A file client for performing file operations.
    force : bool
        Whether to force gunzipping even if files already exist.
    """
    src_dir = strip_hostname(src_dir)
    dest_dir = strip_hostname(dest_dir)

    logger.info(f"Copying ARTATOP inputs from {src_dir} to {dest_dir}")
    directory_listing = file_client.listdir(src_dir, host=src_host)

    files = []
    for file in VASP_OUTPUT_FILES:
        found_file = get_zfile(directory_listing, file, allow_missing=True)
        if found_file is not None:
            files.append(found_file)

    # Copy required files
    copy_files(
        src_dir,
        dest_dir,
        src_host=src_host,
        include_files=files,
        file_client=file_client,
    )

    # Decompress .gz files if necessary
    gunzip_files(
        include_files=files,
        allow_missing=True,
        file_client=file_client,
        force=force,
    )

    logger.info("Finished copying inputs")
