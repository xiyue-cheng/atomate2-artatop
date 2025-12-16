"""Module defining ARTATOP document schemas."""

import gzip
import json
from pathlib import Path
from typing import Any, Optional, Union, List, Tuple

from pydantic import BaseModel, Field
from atomate2_artatop import SETTINGS, __version__
from emmet.core.structure import StructureMetadata
from monty.json import MontyDecoder, jsanitize
import datetime as _dt


def _json_default(o):
    """Fallback JSON serializer for objects not handled by the encoder.

    Currently converts datetime/date/time objects to ISO strings. Any other
    unknown object will be converted to str(o) as a last resort so json.dump
    doesn't raise a TypeError.
    """
    try:
        if isinstance(o, (_dt.datetime, _dt.date, _dt.time)):
            return o.isoformat()
    except Exception:
        pass
    return str(o)


from pymatgen.io.vasp import Vasprun
from pymatgen.core import Structure


from atomate2.utils.datetime import datetime_str

try:
    import ijson
except ImportError:
    ijson = None


class LinearOpticalResponse(BaseModel):
    """Represents the linear optical response data."""

    energy: float = Field(..., description="Energy value in eV.")
    real_part: float = Field(..., description="Real part of the dielectric function.")
    imaginary_part: float = Field(
        ..., description="Imaginary part of the dielectric function."
    )
    absorption_coefficient: Optional[float] = Field(
        None, description="Absorption coefficient at the given energy."
    )
    refractive_index: Optional[float] = Field(
        None, description="Refractive index at the given energy."
    )
    extinction_coefficient: Optional[float] = Field(
        None, description="Extinction coefficient at the given energy."
    )


class NonlinearOpticalResponse(BaseModel):
    """Represents the nonlinear optical response data."""

    energy: float = Field(..., description="Energy value in eV.")
    total_real_part: float = Field(
        ..., description="Total real part of the second harmonic generation response."
    )
    total_imaginary_part: float = Field(
        ...,
        description="Total imaginary part of the second harmonic generation response.",
    )
    contributions: Optional[dict[str, Any]] = Field(
        None, description="Detailed contributions to the nonlinear optical response."
    )


class AtomicContributions(BaseModel):
    atom: str
    orbital_contributions: Optional[dict] = None
    valence_contributions: Optional[dict] = None
    conduction_contributions: Optional[dict] = None
    total_contribution: Optional[float] = None


class DTensorValues(BaseModel):
    energy: float
    components: dict[str, float]


class DeffValues(BaseModel):
    energy: float
    deff: float


class BirefringenceValues(BaseModel):
    energy: float
    delta_n: float


class ARTSummaryEntry(BaseModel):
    atom_type: str
    num_atoms: int
    ind: float
    total: float
    vb: float
    cb: float
    vb_s: float = 0.0
    vb_p: float = 0.0
    vb_d: float = 0.0
    vb_f: float = 0.0
    cb_s: float = 0.0
    cb_p: float = 0.0
    cb_d: float = 0.0
    cb_f: float = 0.0
    tot_s: float = 0.0
    tot_p: float = 0.0
    tot_d: float = 0.0
    tot_f: float = 0.0


class DPmVEntry(BaseModel):
    nbands: float
    value: float
    label: str


class ArtatopInputModel(BaseModel):
    """Definition of input settings for the ARTATOP computation."""

    calc_type: str = Field(
        ..., description="Type of calculation (e.g., linear, nonlinear, ART)."
    )
    input_file: str = Field(..., description="Path to the ARTATOP input file.")
    output_dir: str = Field(..., description="Directory for storing output files.")
    custom_components: Optional[dict[str, Any]] = Field(
        None, description="Custom components used for ART calculations."
    )
    task: str = Field(..., description="Specific ARTATOP task or mode to execute.")

    @classmethod
    def from_file(cls, filename: str) -> "ArtatopInputModel":
        # Just use filename to guess the type
        fname = Path(filename).name.lower()

        if "lin" in fname:
            calc_type = "linear"
        elif "nlin" in fname or "nonlin" in fname:
            calc_type = "nonlinear"
        elif "art" in fname:
            calc_type = "art"
        else:
            raise ValueError(f"Could not detect calc_type from filename: {filename}")

        # Default task and other metadata
        task = "default_task"  # You can later extract this from file if needed
        input_file = str(filename)
        output_dir = str(Path(filename).parent)
        custom_components = None  # Set if needed in the future

        return cls(
            calc_type=calc_type,
            input_file=input_file,
            output_dir=output_dir,
            task=task,
            custom_components=custom_components,
        )


class ArtatopOutputModel(BaseModel):
    """Definition of output results from the ARTATOP computation."""

    dir_name: str = Field(..., description="Directory containing ARTATOP outputs.")

    relaxed_structure: Optional[Structure] = None

    # UV region
    linear_response_uv: Optional[list[LinearOpticalResponse]] = None
    d_tensor_uv: Optional[list[DTensorValues]] = None
    deff_values_uv: Optional[list[DeffValues]] = None
    birefringence_uv: Optional[list[BirefringenceValues]] = None

    # IR region
    linear_response_ir: Optional[list[LinearOpticalResponse]] = None
    d_tensor_ir: Optional[list[DTensorValues]] = None
    deff_values_ir: Optional[list[DeffValues]] = None
    birefringence_ir: Optional[list[BirefringenceValues]] = None

    atomic_contributions: Optional[List[AtomicContributions]] = None

    art_summary: Optional[list[ARTSummaryEntry]] = None
    d_pmV: Optional[list[DPmVEntry]] = None
    dshg_pmV: Optional[list[DPmVEntry]] = None

    d_energy: Optional[List[Tuple[float, float]]] = None  # list of (energy, d-PmV)
    dshg_energy: Optional[
        List[Tuple[float, float]]
    ] = None  # list of (energy, dshg-PmV)

    chemical_formula: Optional[str] = None
    space_group: Optional[str] = None
    point_group: Optional[str] = None
    n_atoms: Optional[int] = None

    ori_a: Optional[float] = None
    ori_b: Optional[float] = None
    ori_c: Optional[float] = None
    ori_alpha: Optional[float] = None
    ori_beta: Optional[float] = None
    ori_gamma: Optional[float] = None
    original_structure: Optional[Structure] = None

    relax_a: Optional[float] = None
    relax_b: Optional[float] = None
    relax_c: Optional[float] = None
    relax_alpha: Optional[float] = None
    relax_beta: Optional[float] = None
    relax_gamma: Optional[float] = None
    relax_volume: Optional[float] = None

    v_over_eg_exp_per_atom: Optional[float] = None
    v_over_eg_hse_per_atom: Optional[float] = None
    bandgap_pbe: Optional[float] = None
    bandgap_exp: Optional[float] = None
    bandgap_hse: Optional[float] = None
    scissor_exp: Optional[float] = None
    scissor_hse: Optional[float] = None

    nbands_static: Optional[int] = None
    nbands_optics: Optional[int] = None
    total_energy_static: Optional[float] = None
    enthalpy: Optional[float] = Field(None, description="Enthalpy per atom (eV/atom)")

    kpoints_relax1: Optional[list[int]] = None
    kpoints_relax2: Optional[list[int]] = None
    kpoints_static: Optional[list[int]] = None
    kpoints_optics: Optional[list[int]] = None
    kpoints_hse06: Optional[list[int]] = None

    ediffg_relax2: Optional[float] = None
    aexx_hse: Optional[float] = None

    art_highest_component: Optional[str] = None
    art_highest_value: Optional[float] = None

    @classmethod
    def from_directory(
        cls,
        dir_name: str,
        input_file: str,
        job_paths: Optional[dict] = None,
        pbe_vasprun_file: Optional[Path] = None,
        hse_vasprun_file: Optional[Path] = None,
        additional_metadata: dict = None,
    ) -> "ArtatopOutputModel":
        from atomate2_artatop.artatop_parser import (
            parse_artatop_outputs,
            get_artatop_functional_data,
        )

        return parse_artatop_outputs(
            dir_name=dir_name,
            input_file=input_file,
            job_paths=job_paths,
            pbe_vasprun_file=pbe_vasprun_file,
            hse_vasprun_file=hse_vasprun_file,
            additional_metadata=additional_metadata,
        )

    class Config:
        arbitrary_types_allowed = True
        json_encoders = {Structure: lambda v: v.as_dict()}


class ArtatopTaskDocument(StructureMetadata, extra="allow"):
    """Main schema for an ARTATOP task document."""

    # Ensure builder_meta exists so downstream code can set/update it
    builder_meta: Optional[dict[str, Any]] = Field(
        default=None, description="Builder metadata such as source and version."
    )

    dir_name: Union[str, Path] = Field(
        ..., description="Directory containing ARTATOP outputs."
    )

    input_data: ArtatopInputModel = Field(
        ..., description="Input parameters for the ARTATOP computation."
    )
    output_data: ArtatopOutputModel = Field(
        ..., description="Parsed results from the ARTATOP computation."
    )

    additional_metadata: dict[str, Any] = Field(
        default_factory=dict,
        description="Additional optional metadata related to the task.",
    )
    last_updated: str = Field(
        default_factory=datetime_str,
        description="Timestamp when this task document was last updated.",
    )

    @classmethod
    def from_directory(
        cls,
        dir_name: str,
        input_file: str,
        pbe_vasprun_file: Optional[Path] = None,
        hse_vasprun_file: Optional[Path] = None,
        store_additional_json: bool = SETTINGS.ARTATOP_STORE_ADDITIONAL_JSON,
        additional_metadata: dict = None,
        job_paths: Optional[dict] = None,
    ) -> "ArtatopTaskDocument":
        """Parse ARTATOP inputs and outputs, then construct the task document."""
        from pymatgen.io.vasp import Vasprun

        # --- Load ARTATOP input model ---
        input_data = ArtatopInputModel.from_file(input_file)

        # --- Load ARTATOP parsed outputs ---
        output_data = ArtatopOutputModel.from_directory(
            dir_name=dir_name,
            input_file=input_file,
            job_paths=job_paths,
            pbe_vasprun_file=pbe_vasprun_file,
            hse_vasprun_file=hse_vasprun_file,
            additional_metadata=additional_metadata,
        )

        if additional_metadata is None:
            additional_metadata = {}

        structure = output_data.relaxed_structure

        base_doc = cls(
            structure=structure,
            dir_name=str(dir_name),
            input_data=input_data,
            output_data=output_data,
            additional_metadata=additional_metadata,
        )

        if base_doc.builder_meta is None:
            base_doc.builder_meta = {}

        builder_dict = (
            base_doc.builder_meta.dict()
            if hasattr(base_doc.builder_meta, "dict")
            else base_doc.builder_meta
        )

        base_doc.builder_meta = {
            **builder_dict,
            "source": "artatop",
            "version": __version__,
        }

        # --- Optionally save JSON summary ---
        if store_additional_json:
            art_json_path = Path(dir_name) / "artatop_summary.json.gz"
            with gzip.open(art_json_path, "wt", encoding="UTF-8") as file:
                file.write("[")
                blocks = [
                    {"output_data": output_data},
                    {"additional_metadata": additional_metadata},
                    {"builder_meta": base_doc.builder_meta},
                ]
                for i, block in enumerate(blocks):
                    json.dump(
                        jsanitize(block, strict=False, allow_bson=True),
                        file,
                        default=_json_default,
                    )
                    if i < len(blocks) - 1:
                        file.write(",")
                file.write("]")

        # --- Return the task document ---
        return base_doc

    class Config:
        arbitrary_types_allowed = True
        json_encoders = {Structure: lambda v: v.as_dict()}


def read_saved_json(
    filename: str, pymatgen_objs: bool = True, query: str = "relaxed_structure"
) -> dict[str, Any]:
    """
    Read the data from .json.gz files corresponding to query.

    Uses ijson to parse specific keys (memory efficient)

    Parameters
    ----------
    filename: str
        Path to the JSON file (gzipped).
    pymatgen_objs: bool
        If True, convert any pymatgen-compatible dicts (e.g., Structure) to actual objects.
    query: str or None
        Field name to query from the JSON file. If None, load all fields.

    Returns
    -------
    dict
        Dictionary of the requested field(s) from the JSON file.
    """
    with gzip.open(filename, "rb") as file:
        artatop_data = {
            field: data
            for obj in ijson.items(file, "item", use_float=True)
            for field, data in obj.items()
            if query is None or query == field
        }

        if not artatop_data:
            raise ValueError(
                f"No data associated with 'query={query}' found in the JSON file. "
                "Please check your query string."
            )

    if pymatgen_objs:
        decoder = MontyDecoder()
        for key, value in artatop_data.items():
            artatop_data[key] = decoder.process_decoded(value)

    return artatop_data


def save_to_json(self, filename: str) -> None:
    """Save the task document as a compressed JSON file."""
    with gzip.open(filename, "wt", encoding="UTF-8") as file:
        json.dump(self.dict(), file, indent=4, default=_json_default)
