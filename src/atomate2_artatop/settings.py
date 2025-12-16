"""Settings for the atomate2-artatop extension."""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Literal, Union

from monty.serialization import loadfn
from pydantic import Field, model_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

_DEFAULT_CONFIG_FILE_PATH = "~/.artatop.yaml"
_ENV_PREFIX = "artatop_"


class ArtatopSettings(BaseSettings):
    """
    Settings for running ARTATOP outside of the main atomate2 package.

    Settings can be provided via a YAML file (defaults to ~/.atomate2_artatop.yaml)
    or environment variables prefixed with ``ARTATOP_`` (e.g. ``ARTATOP_ARTATOP_CMD``).
    """

    CONFIG_FILE: str = Field(
        _DEFAULT_CONFIG_FILE_PATH, description="File to load alternative defaults from."
    )

    ARTATOP_CMD: str = Field(
        default="artatop", description="Command to run the ARTATOP binary."
    )
    ARTATOP_CUSTODIAN_MAX_ERRORS: int = Field(
        5, description="Maximum number of errors to correct before custodian gives up."
    )
    ARTATOP_ZIP_FILES: Union[bool, Literal["atomate"]] = Field(
        "atomate",
        description=(
            "If True compress all files, 'atomate' compresses a curated subset, "
            "and False skips compression."
        ),
    )
    ARTATOP_STORE_ADDITIONAL_JSON: bool = Field(
        default=True,
        description=(
            "Ingest any additional JSON data present when parsing ARTATOP directories."
        ),
    )

    model_config = SettingsConfigDict(env_prefix=_ENV_PREFIX)

    @model_validator(mode="before")
    @classmethod
    def load_default_settings(cls, values: dict[str, Any]) -> dict[str, Any]:
        """Load settings from YAML config or environment variables."""
        config_file_path = values.get(key := "CONFIG_FILE", _DEFAULT_CONFIG_FILE_PATH)
        env_var_name = f"{_ENV_PREFIX.upper()}{key}"
        config_file_path = Path(config_file_path).expanduser()

        new_values: dict[str, Any] = {}
        if config_file_path.exists():
            if config_file_path.stat().st_size == 0:
                warnings.warn(
                    f"Using {env_var_name} at {config_file_path} but it's empty",
                    stacklevel=2,
                )
            else:
                try:
                    new_values.update(loadfn(config_file_path))
                except ValueError:
                    raise SyntaxError(
                        f"{env_var_name} at {config_file_path} is unparsable"
                    ) from None
        elif config_file_path != Path(_DEFAULT_CONFIG_FILE_PATH).expanduser():
            warnings.warn(
                f"{env_var_name} at {config_file_path} does not exist", stacklevel=2
            )

        return {**new_values, **values}


SETTINGS = ArtatopSettings()

__all__ = ("SETTINGS", "ArtatopSettings")
