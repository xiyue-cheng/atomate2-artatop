import os
import logging
from custodian.custodian import Validator

logger = logging.getLogger(__name__)

class ArtatopFilesValidator(Validator):
    """
    Validates the existence of required ARTATOP output files.
    """

    def __init__(self, required_files=None) -> None:
        """
        Args:
            required_files: List of required ARTATOP output files.
        """
        self.required_files = required_files or []
        
    def check(self, directory: str = "./") -> bool:
        missing_files = [f for f in self.required_files if not os.path.isfile(os.path.join(directory, f))]
        if missing_files:
            logger.error(f"Missing ARTATOP output files: {missing_files}")
            return True  # Error
        return False  # No error
        


