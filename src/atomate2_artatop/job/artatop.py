"""
Module defines job makers for ARTATOP workflows.

It includes:
- get_artatop_jobs: Creates ARTATOP jobs for optical response analysis.
"""
from __future__ import annotations
import os
import shutil
from pathlib import Path
import logging

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from jobflow import Maker, Response, job, Flow
from atomate2_artatop.sets.core import InputFileHandler
from atomate2.common.files import copy_files
from pymatgen.io.vasp import Vasprun
from monty.os.path import zpath
from atomate2.utils.path import strip_hostname

from atomate2.vasp.sets.core import NonSCFSetGenerator
from atomate2_artatop.jobs import ARTATOPMaker
from atomate2_artatop.files import copy_artatop_files
from atomate2_artatop.schemas import (
    ArtatopTaskDocument,
    ArtatopInputModel,
    ArtatopOutputModel,
)

if TYPE_CHECKING:
    from pathlib import Path
    from pymatgen.core import Structure

logger = logging.getLogger(__name__)


@job
def get_artatop_jobs(
    optics_job_output,
    hse_job_output,
    job_paths: dict | None = None,
    artatop_maker: ARTATOPMaker | None = None,
    additional_metadata: dict | None = None,
) -> Response:
    """
    Create a list of ARTATOP jobs (LIN, NLIN, ART) dynamically.

    Parameters
    ----------
    optics_job_output : TaskDoc
        The output reference from OpticsMaker, which is a TaskDoc object.
    hse_job_output : TaskDoc
        The output reference from HSE calculation.
    job_paths : dict or None
        Dictionary of job paths for reference.
    artatop_maker : ARTATOPMaker or None
        Maker for the ARTATOP jobs.
    additional_metadata : dict or None
        Additional metadata for the ARTATOP task document.

    Returns
    -------
    Response
        A response containing the ARTATOP jobs and output directories.
    """
    # Both static calculations are done now, so files exist

    pbe_path = strip_hostname(optics_job_output.dir_name)
    hse_path = strip_hostname(hse_job_output.dir_name)

    pbe_file = zpath(os.path.join(pbe_path, "vasprun.xml.gz"))
    hse_file = zpath(os.path.join(hse_path, "vasprun.xml.gz"))

    pbe_vasprun = Vasprun(pbe_file, parse_dos=False)
    hse_vasprun = Vasprun(hse_file, parse_dos=False)

    gap_pbe = pbe_vasprun.eigenvalue_band_properties[0]
    gap_hse = hse_vasprun.eigenvalue_band_properties[0]
    scissor = round(gap_hse - gap_pbe, 3)

    outputs = {
        "optics_dir": optics_job_output.dir_name,
        "artatop_dirs": [],
        "artatop_task_documents": [],
    }

    artatop_maker = artatop_maker or ARTATOPMaker()
    artatop_maker.scissor = scissor

    # Add vasprun paths to be passed later to from_directory()
    artatop_maker.task_document_kwargs.update(
        {
            "pbe_vasprun_file": str(pbe_file),
            "hse_vasprun_file": str(hse_file),
            "job_paths": job_paths,
        }
    )
    artatop_maker.task_document_kwargs["additional_metadata"] = additional_metadata

    # Run ARTATOP in the new directory
    artatop_job = artatop_maker.make(art_input_dir=optics_job_output.dir_name)

    # Store job details
    outputs["artatop_dirs"].append(str(Path(".")))
    outputs["artatop_task_documents"].append(artatop_job.output)

    return Response(replace=Flow(artatop_job, output=outputs))
