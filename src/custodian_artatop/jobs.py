import logging
import os
import shlex
import shutil
import subprocess
from monty.io import zopen
from monty.shutil import compress_file
from custodian.custodian import Job

ARTATOPINPUT_FILES = ["input_lin", "input_nlin", "input_art"]
ARTATOP_OUTPUT_FILES = ["re_lin", "re_nlin", "re_art"]
FW_FILES = ["custodian.json", "FW.json", "FW_submit.script"]

logger = logging.getLogger(__name__)

class ArtatopJob(Job):
    """Runs the ARTATOP Job."""

    def __init__(
        self,
        artatop_cmd: str,
        output_file: str = "artatopout",
        stderr_file: str = "std_err_artatop.txt",
        gzipped: bool = True,
        add_files_to_gzip=None,
        backup: bool = True,
    ) -> None:
        """
        Args:
            artatop_cmd: Command to run ARTATOP.
            output_file: Standard output file name.
            stderr_file: Standard error file name.
            gzipped: If True, ARTATOP files will be gzipped.
            add_files_to_gzip: List of additional files to gzip.
            backup: If True, input files will be backed up.
        """
        self.artatop_cmd = artatop_cmd
        self.output_file = output_file
        self.stderr_file = stderr_file
        self.gzipped = gzipped
        self.add_files_to_gzip = add_files_to_gzip or []
        self.backup = backup

    def setup(self, directory="./") -> None:
        """Back up input files if required."""
        if self.backup:
            for file in ARTATOPINPUT_FILES:
                src = os.path.join(directory, file)
                if os.path.isfile(src):
                    backup_file = f"{src}.orig"
                    shutil.copy(src, backup_file)
                    logger.info(f"Backed up {file} to {backup_file}")

    def run(self, directory="./"):
        """Run the ARTATOP command."""
        cmd = self.artatop_cmd if isinstance(self.artatop_cmd, str) else shlex.join(self.artatop_cmd)
        logger.info(f"Running ARTATOP command: {cmd}")

        with zopen(os.path.join(directory, self.output_file), "wt") as f_std, \
                zopen(os.path.join(directory, self.stderr_file), "wt") as f_err:
            result = subprocess.run(cmd, stdout=f_std, stderr=f_err, shell=True)
            if result.returncode != 0:
                raise RuntimeError(f"ARTATOP failed with return code {result.returncode}")

    def postprocess(self, directory="./") -> None:
        """Compress output files if required."""
        if self.gzipped:
            for file in ARTATOP_OUTPUT_FILES:
                file_path = os.path.join(directory, file)
                if os.path.isfile(file_path):
                    compress_file(file_path, compression="gz")
            for file in self.add_files_to_gzip:
                compress_file(os.path.join(directory, file), compression="gz")
            logger.info("Compressed output files.")
