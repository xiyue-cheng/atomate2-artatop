"""Functions to run ARTATOP."""

from __future__ import annotations

import logging
import shlex
import subprocess
from os.path import expandvars
from typing import TYPE_CHECKING, Any

from custodian import Custodian
from custodian.artatop.handlers import ArtatopFilesValidator
from custodian.artatop.jobs import ArtatopJob
from jobflow.utils import ValueEnum
from atomate2.settings import Atomate2Settings

AT2_SETTINGS = Atomate2Settings()
from atomate2_artatop.settings import SETTINGS

if TYPE_CHECKING:
    from collections.abc import Sequence

    from custodian.custodian import Validator


_DEFAULT_VALIDATORS = [ArtatopFilesValidator()]
_DEFAULT_HANDLERS = ()

logger = logging.getLogger(__name__)


class JobType(ValueEnum):
    """
    Type of ARTATOP job.

    - ``DIRECT``: Run Artatop without custodian.
    - ``NORMAL``: Normal custodian :obj:`.ArtatopJob`.
    """

    DIRECT = "direct"
    NORMAL = "normal"


def run_artatop(
    job_type: JobType | str = JobType.NORMAL,
    artatop_cmd: str = SETTINGS.ARTATOP_CMD,
    max_errors: int = SETTINGS.ARTATOP_CUSTODIAN_MAX_ERRORS,
    scratch_dir: str = AT2_SETTINGS.CUSTODIAN_SCRATCH_DIR,
    validators: Sequence[Validator] = _DEFAULT_VALIDATORS,
    artatop_job_kwargs: dict[str, Any] = None,
    custodian_kwargs: dict[str, Any] = None,
) -> None:
    """
    Run ARTATOP with Custodian support.
    Examples
    --------
    >>> run_artatop(job_type="normal", artatop_cmd="artatop < input > output")

    Parameters
    ----------
    job_type : str or .JobType
        The job type ('direct' or 'normal').
    artatop_cmd : str
        Command to execute ARTATOP. Environment variables in the path will be expanded.
    max_errors : int
        Maximum number of errors allowed by Custodian.
    scratch_dir : str
        Directory for scratch files.
    validators : list of Validator
        Validators for Custodian.
    artatop_job_kwargs : dict
        Keyword arguments for ArtatopJob.
    custodian_kwargs : dict
        Keyword arguments for Custodian.
    """
    artatop_job_kwargs = artatop_job_kwargs or {}
    custodian_kwargs = custodian_kwargs or {}

    artatop_cmd = expandvars(artatop_cmd)
    split_artatop_cmd = shlex.split(artatop_cmd)

    if job_type == JobType.DIRECT:
        logger.info(f"Running command: {artatop_cmd}")
        return_code = subprocess.call(artatop_cmd, shell=True)  # noqa: S602
        logger.info(f"Command finished with return code {return_code}")
        return

    if job_type == JobType.NORMAL:
        jobs = [ArtatopJob(**artatop_job_kwargs)]
    else:
        raise ValueError(f"Unsupported job type: {job_type}")

    handlers: list = []

    # Validate that max_errors is an integer
    if not isinstance(max_errors, int):
        raise ValueError(f"`max_errors` should be an integer, got {type(max_errors)}")

    c = Custodian(
        handlers,
        jobs,
        validators=validators,
        max_errors=max_errors,
        scratch_dir=str(scratch_dir),
        **custodian_kwargs,
    )

    logger.info("Running ARTATOP with Custodian.")
    c.run()
