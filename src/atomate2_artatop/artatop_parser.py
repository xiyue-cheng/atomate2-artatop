"""
ARTATOP output parser for optical properties and atomic contributions.
"""
import logging
import numpy as np
import pathlib
from pathlib import Path
import warnings
from typing import Union, List, Optional
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar, Vasprun
from pymatgen.io.vasp.inputs import Incar
from atomate2.utils.path import strip_hostname
from atomate2_artatop.schemas import (
    LinearOpticalResponse,
    AtomicContributions,
    DTensorValues,
    DeffValues,
    BirefringenceValues,
    ARTSummaryEntry,
    DPmVEntry,
    ArtatopOutputModel,
)

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from atomate2_artatop.outputs import (
    write_artatop_summary,
    write_artatop_summary_result_line,
    parse_fixed_lines,
    write_artatop_analysis_file,
)

logger = logging.getLogger(__name__)


def read_value_from_line(path: Path, line_number: int, col: int = 2) -> float:
    """
    Read a specific value from a file by line number and column.

    Parameters
    ----------
    path : Path
        Path to the file to read.
    line_number : int
        Line number to read (1-indexed).
    col : int
        Column index to extract (0-indexed).

    Returns
    -------
    float
        Extracted numerical value, or 0.0 if line doesn't exist.
    """
    with open(path) as f:
        lines = f.readlines()
        if line_number <= len(lines):
            return float(lines[line_number - 1].split()[col])
    return 0.0


def parse_d_tensor_and_deff(
    nonlin_dir: Path, line_idx_0: int, line_idx_1: int, energies: list[float]
) -> tuple[list[DTensorValues], list[DeffValues]]:
    """
    Parse nonlinear optical tensors and effective coefficients.

    Reads d-tensor components from ARTATOP output files and computes
    effective nonlinear coefficients using Kleinman symmetry relations.

    Parameters
    ----------
    nonlin_dir : Path
        Directory containing nonlinear output files (nonlin_*.dat).
    line_idx_0 : int
        Line index for first energy point.
    line_idx_1 : int
        Line index for second energy point.
    energies : list[float]
        Energy values corresponding to the line indices.

    Returns
    -------
    tuple[list[DTensorValues], list[DeffValues]]
        D-tensor components and effective coefficients for both energies.
    """
    components = [
        "xxx",
        "xyy",
        "xzz",
        "xyz",
        "xxz",
        "xxy",
        "yxx",
        "yyy",
        "yzz",
        "yyz",
        "yxz",
        "yxy",
        "zxx",
        "zyy",
        "zzz",
        "zyz",
        "zxz",
        "zxy",
    ]

    def compute_d_tensor(line_idx: int) -> dict[str, float]:
        d = {}
        for comp in components:
            path = nonlin_dir / f"nonlin_{comp}.dat"
            val = read_value_from_line(path, line_idx)
            d[f"d_{comp}"] = (val / 2.0) * (4 * np.pi / 3.0) * 10
        return d

    def compute_deff(d):

        d11, d22, d33 = d["d_xxx"], d["d_yyy"], d["d_zzz"]
        d12, d13, d21, d23, d31, d32 = (
            d["d_xyy"],
            d["d_xzz"],
            d["d_yxx"],
            d["d_yzz"],
            d["d_zxx"],
            d["d_zyy"],
        )
        d16, d15, d26, d24, d35, d34 = (
            d["d_xxy"],
            d["d_xxz"],
            d["d_yxy"],
            d["d_yyz"],
            d["d_zxz"],
            d["d_zyz"],
        )
        d14, d25, d36 = d["d_xyz"], d["d_yxz"], d["d_zxy"]
        dtemp = (
            19 / 105 * (d11 ** 2 + d22 ** 2 + d33 ** 2)
            + 13
            / 105
            * (d11 * d12 + d11 * d13 + d22 * d21 + d22 * d23 + d33 * d31 + d33 * d32)
            + 44
            / 105
            * (d16 ** 2 + d15 ** 2 + d26 ** 2 + d24 ** 2 + d35 ** 2 + d34 ** 2)
            + 13
            / 105
            * (d16 * d23 + d15 * d32 + d26 * d13 + d24 * d31 + d35 * d12 + d34 * d21)
            + 5 / 7 * ((d14 + d25 + d36) / 3) ** 2
        )
        return np.sqrt(abs(dtemp))

    d0 = compute_d_tensor(line_idx_0)
    d1 = compute_d_tensor(line_idx_1)
    return (
        [
            DTensorValues(energy=energies[0], components=d0),
            DTensorValues(energy=energies[1], components=d1),
        ],
        [
            DeffValues(energy=energies[0], deff=compute_deff(d0)),
            DeffValues(energy=energies[1], deff=compute_deff(d1)),
        ],
    )


def parse_linear_optics(
    lin_dir: Path, line_idx_0: int, line_idx_1: int, energies: list[float]
) -> tuple[list[LinearOpticalResponse], list[BirefringenceValues]]:
    """
    Parse linear optical properties from ARTATOP output.

    Extracts dielectric function components and computes refractive indices
    and birefringence for specified energy points.

    Parameters
    ----------
    lin_dir : Path
        Directory containing linear output files (lin_*.dat).
    line_idx_0 : int
        Line index for first energy point.
    line_idx_1 : int
        Line index for second energy point.
    energies : list[float]
        Energy values corresponding to the line indices.

    Returns
    -------
    tuple[list[LinearOpticalResponse], list[BirefringenceValues]]
        Linear optical responses and birefringence values for both energies.
    """
    comps = ["xx", "yy", "zz"]

    def extract_eps(path: Path, idx: int) -> tuple[float, float]:
        with open(path) as f:
            lines = f.readlines()
            if idx <= len(lines):
                parts = lines[idx - 1].split()
                return float(parts[2]), float(parts[1])  # real, imag
        return 0.0, 0.0

    def compute_n(eps_real, eps_imag):
        return np.sqrt((np.sqrt(eps_real ** 2 + eps_imag ** 2) + eps_real) / 2)

    def compute_biref(nvals):
        return max(nvals) - min(nvals)

    results = []
    biref_list = []
    for idx, e in zip([line_idx_0, line_idx_1], energies):
        nvals = []
        for comp in comps:
            eps_r, eps_i = extract_eps(lin_dir / f"lin_{comp}.dat", idx)
            n = compute_n(eps_r, eps_i)
            results.append(
                LinearOpticalResponse(
                    energy=e,
                    real_part=eps_r,
                    imaginary_part=eps_i,
                    refractive_index=n,
                    absorption_coefficient=None,
                )
            )
            nvals.append(n)
        biref_list.append(BirefringenceValues(energy=e, delta_n=compute_biref(nvals)))
    return results, biref_list


def parse_optical_response(dir_name: str) -> ArtatopOutputModel:
    """
    Parse complete optical response from ARTATOP calculation directory.

    Processes both UV and IR regions for linear and nonlinear optical properties.

    Parameters
    ----------
    dir_name : str
        Directory containing ARTATOP output files.

    Returns
    -------
    ArtatopOutputModel
        Complete parsed optical response data.
    """
    dir_path = Path(dir_name)
    if not dir_path.exists():
        return ArtatopOutputModel(dir_name=str(dir_path))

    uv_energies = [0.0, 1.167]
    ir_energies = [0.0, 0.65]

    d_tensor_uv, deff_uv = parse_d_tensor_and_deff(
        dir_path / "out_nonlin", 10, 110, uv_energies
    )
    d_tensor_ir, deff_ir = parse_d_tensor_and_deff(
        dir_path / "out_nonlin", 10, 66, ir_energies
    )

    lin_resp_uv, biref_uv = parse_linear_optics(
        dir_path / "out_lin", 8, 108, uv_energies
    )
    lin_resp_ir, biref_ir = parse_linear_optics(
        dir_path / "out_lin", 8, 64, ir_energies
    )

    return ArtatopOutputModel(
        dir_name=str(dir_path),
        linear_response_uv=lin_resp_uv,
        birefringence_uv=biref_uv,
        d_tensor_uv=d_tensor_uv,
        deff_values_uv=deff_uv,
        linear_response_ir=lin_resp_ir,
        birefringence_ir=biref_ir,
        d_tensor_ir=d_tensor_ir,
        deff_values_ir=deff_ir,
    )


def parse_full_spectrum_optics(dir_path: Path) -> None:
    """
    Generate full spectral data files from ARTATOP linear outputs.

    Creates REAL.dat, IMAG.dat, refractive.dat, and birefringence.dat
    with complete spectral data for dielectric function analysis.

    Parameters
    ----------
    dir_path : Path
        Directory containing ARTATOP output files.

    Returns
    -------
    None
        Writes spectral data files to the directory.
    """
    lin_dir = dir_path / "out_lin"
    output_dir = dir_path

    def read_data(file: Path):
        with open(file, "r") as f:
            lines = f.readlines()[7:2578]
            data = []
            for line in lines:
                if line.strip().startswith("#") or not line.strip():
                    continue
                try:
                    data.append(list(map(float, line.split())))
                except ValueError:
                    continue
            return data

    esp_xx = read_data(lin_dir / "lin_xx.dat")
    esp_yy = read_data(lin_dir / "lin_yy.dat")
    esp_zz = read_data(lin_dir / "lin_zz.dat")

    imag_path = output_dir / "IMAG.dat"
    real_path = output_dir / "REAL.dat"
    refr_path = output_dir / "refractive.dat"
    biref_path = output_dir / "birefringence.dat"

    with open(imag_path, "w") as imag_file, open(real_path, "w") as real_file:
        for x, y, z in zip(esp_xx, esp_yy, esp_zz):
            energy = x[0]
            imag_file.write(f"{energy:.5f} {x[1]:.5f} {y[1]:.5f} {z[1]:.5f}\n")
            real_file.write(f"{energy:.5f} {x[2]:.5f} {y[2]:.5f} {z[2]:.5f}\n")

    imag = np.loadtxt(imag_path)
    real = np.loadtxt(real_path)

    with open(refr_path, "w") as refr_file, open(biref_path, "w") as bire_file:
        for i in range(len(imag)):
            energy = imag[i][0]
            im_xx, im_yy, im_zz = imag[i][1:]
            re_xx, re_yy, re_zz = real[i][1:]

            refra_xx = np.sqrt((np.sqrt(im_xx ** 2 + re_xx ** 2) + re_xx) / 2)
            refra_yy = np.sqrt((np.sqrt(im_yy ** 2 + re_yy ** 2) + re_yy) / 2)
            refra_zz = np.sqrt((np.sqrt(im_zz ** 2 + re_zz ** 2) + re_zz) / 2)
            refra = (refra_xx + refra_yy + refra_zz) / 3

            refr_file.write(
                f"{energy:.5f} {refra_xx:.5f} {refra_yy:.5f} {refra_zz:.5f} {refra:.5f}\n"
            )

            delta = max(refra_xx, refra_yy, refra_zz) - min(
                refra_xx, refra_yy, refra_zz
            )
            bire_file.write(f"{energy:.5f} {delta:.5f}\n")


def write_result_re(output: ArtatopOutputModel, filename: str = "result.re") -> None:
    """
    Write formatted optical response summary to result.re file.

    Creates a human-readable summary of optical properties in the
    standard ARTATOP result.re format.

    Parameters
    ----------
    output : ArtatopOutputModel
        Parsed optical response data.
    filename : str
        Output filename for the summary.

    Returns
    -------
    None
        Writes formatted summary to file.
    """

    def write_block(f, label, d_tensor, deff_values, lin_resp, biref):
        f.write(f"\n{label}\n")

        # write deff values (grouped by energy)
        for deff in deff_values:
            f.write(f"deff at omega = {deff.energy:.3f} eV\n{deff.deff:.3f}\n\n")

        # write d-tensor (grouped by energy)
        for dset in d_tensor:
            f.write(f"d at omega = {dset.energy:.3f} eV\n")
            d = dset.components
            rows = [
                [
                    d.get("d_xxx", 0),
                    d.get("d_xyy", 0),
                    d.get("d_xzz", 0),
                    d.get("d_xyz", 0),
                    d.get("d_xxz", 0),
                    d.get("d_xxy", 0),
                ],
                [
                    d.get("d_yxx", 0),
                    d.get("d_yyy", 0),
                    d.get("d_yzz", 0),
                    d.get("d_yyz", 0),
                    d.get("d_yxz", 0),
                    d.get("d_yxy", 0),
                ],
                [
                    d.get("d_zxx", 0),
                    d.get("d_zyy", 0),
                    d.get("d_zzz", 0),
                    d.get("d_zyz", 0),
                    d.get("d_zxz", 0),
                    d.get("d_zxy", 0),
                ],
            ]
            for row in rows:
                f.write(" ".join(f"{val:.3f}" for val in row) + "\n")
            f.write("\n")

        # write linear optic values
        energies = sorted(set(l.energy for l in lin_resp))
        for energy in energies:
            lin_vals = [l for l in lin_resp if l.energy == energy]
            biref_val = next((b for b in biref if b.energy == energy), None)

            if lin_vals:
                f.write(f"linear optic epsilon at {energy:.3f} eV\n")
                f.write(" ".join(f"{l.real_part:.3f}" for l in lin_vals) + "\n\n")

                f.write(f"refractive n at {energy:.3f} eV\n")
                f.write(
                    " ".join(f"{l.refractive_index:.3f}" for l in lin_vals) + "\n\n"
                )

            if biref_val:
                f.write(
                    f"birefringence at {energy:.3f} eV\n{biref_val.delta_n:.3f}\n\n"
                )

    # Only open the file once
    with open(filename, "w") as f:
        write_block(
            f,
            "UV Region",
            output.d_tensor_uv,
            output.deff_values_uv,
            output.linear_response_uv,
            output.birefringence_uv,
        )
        write_block(
            f,
            "IR Region",
            output.d_tensor_ir,
            output.deff_values_ir,
            output.linear_response_ir,
            output.birefringence_ir,
        )


def detect_orbital_type(all_lines: list[str], atoms_count: int) -> int:
    """
    Detect orbital basis type from ARTATOP output.

    Parameters
    ----------
    all_lines : list[str]
        Lines from ARTATOP orbital output file.
    atoms_count : int
        Number of atoms in the structure.

    Returns
    -------
    int
        Orbital type: 3 (spd), 4 (spdf), or 9 (lorbit=11). Returns 0 if unknown.

    """
    # Skip header (9 lines) and empty lines
    orbital_lines = [line for line in all_lines[9:] if line.strip()]
    orbitals_per_atom = len(orbital_lines) / atoms_count

    if orbitals_per_atom == 3:
        return 3  # spd
    elif orbitals_per_atom == 4:
        return 4  # spdf
    elif orbitals_per_atom == 9:
        return 9  # lorbit=11
    else:
        logger.info(
            f"Unknown orbital configuration: {orbitals_per_atom} orbitals per atom"
        )
        return 0


def parse_orbital_atomic_contributions(
    structure: Structure, out_dir: Path
) -> list[AtomicContributions]:
    """
    Parse atomic orbital contributions to nonlinear optical response.

    Processes valence, conduction, and total orbital contributions
    from ARTATOP atomic resolution projection files.

    Parameters
    ----------
    structure : Structure
        Crystal structure for atom labeling.
    out_dir : Path
        Directory containing ARTATOP atomic projection files.

    Returns
    -------
    list[AtomicContributions]
        Atomic contributions with orbital resolution.
    """
    val_path = next(out_dir.glob("arp_shg_val_*.txt"))
    con_path = next(out_dir.glob("arp_shg_con_*.txt"))
    all_path = out_dir / "arp_nonlin.txt"

    def read_clean_lines(path: Path, skip_lines: int = 0) -> list[str]:
        with path.open() as f:
            return [
                line
                for line in f.readlines()[skip_lines:]
                if line.strip()
                and not line.lstrip().startswith("#")
                and not line.lstrip().lower().startswith("tot")
            ]

    val_lines = read_clean_lines(val_path, skip_lines=1)
    con_lines = read_clean_lines(con_path, skip_lines=1)

    atomic_contribs = []

    # Read flattened arp_nonlin (skip header 9 lines)
    all_lines = all_path.read_text().splitlines()[9:]
    clean_lines = [
        ln
        for ln in all_lines
        if ln.strip()
        and not ln.lstrip().startswith("#")
        and not ln.lstrip().lower().startswith("tot")
    ]

    total_sum = max(
        sum(float(line.split()[2]) for line in clean_lines) if clean_lines else 0.0,
        1e-10,
    )

    n_atoms = len(structure)
    if n_atoms == 0:
        logger.info("Structure has no sites; skipping orbital parsing")
        return atomic_contribs

    idx = 0
    base_labels = ["s", "p", "d"]

    for i, site in enumerate(structure.sites):
        # read 3 lines per atom sequentially (spd)
        values = []
        for j in range(3):
            if idx + j < len(clean_lines):
                try:
                    values.append(float(clean_lines[idx + j].split()[2]))
                except Exception:
                    values.append(0.0)
            else:
                values.append(0.0)

        # read corresponding val and con lines defensively
        v_vals = []
        c_vals = []
        for j in range(3):
            try:
                v_vals.append(float(val_lines[idx + j].split()[2]))
            except Exception:
                v_vals.append(0.0)
            try:
                c_vals.append(float(con_lines[idx + j].split()[2]))
            except Exception:
                c_vals.append(0.0)

        idx += 3

        orb_acc = {lbl: val for lbl, val in zip(base_labels, values)}
        val_acc = {lbl: val for lbl, val in zip(base_labels, v_vals)}
        con_acc = {lbl: val for lbl, val in zip(base_labels, c_vals)}

        orb_pct = {k: 100.0 * v / total_sum for k, v in orb_acc.items()}
        val_pct = {k: 100.0 * v / total_sum for k, v in val_acc.items()}
        con_pct = {k: 100.0 * v / total_sum for k, v in con_acc.items()}

        atomic_contribs.append(
            AtomicContributions(
                atom=f"{i+1}-{site.species_string}",
                orbital_contributions=orb_pct,
                valence_contributions=val_pct,
                conduction_contributions=con_pct,
                total_contribution=sum(orb_pct.values()),
            )
        )
    logger.info(f"Parsed orbital contributions for {len(atomic_contribs)} atoms")
    return atomic_contribs


def write_result_art_IND(
    summary_entries: list[ARTSummaryEntry],
    filename: Union[str, Path] = "result.art_IND",
) -> None:
    """
    Write atomic resolution tensor (ART) individual contributions to file.

    Creates result.art_IND file with per-atom-type contributions normalized
    by number of atoms for easy comparison across different structures.

    Parameters
    ----------
    summary_entries : list[ARTSummaryEntry]
        ART summary entries with atomic contributions.
    filename : str or Path
        Output filename for individual contributions.

    Returns
    -------
    None
        Writes normalized ART contributions to file.
    """
    filename = Path(filename)
    with open(filename, "w") as f:
        f.write(
            "Type Natom IND TOT VB CB VB_s VB_p VB_d CB_s CB_p CB_d TOT_s TOT_p TOT_d\n"
        )
        for e in summary_entries:
            nat = e.num_atoms
            f.write(
                f"{e.atom_type} {nat} "
                f"{e.ind:.4f} {e.total:.4f} "
                f"{e.vb / nat:.4f} {e.cb / nat:.4f} "
                f"{e.vb_s / nat:.4f} {e.vb_p / nat:.4f} {e.vb_d / nat:.4f} "
                f"{e.cb_s / nat:.4f} {e.cb_p / nat:.4f} {e.cb_d / nat:.4f} "
                f"{e.tot_s / nat:.4f} {e.tot_p / nat:.4f} {e.tot_d / nat:.4f}\n"
            )


def parse_art_summary_and_dpmv(
    poscar_path: Path, out_dir: Path
) -> tuple[list[ARTSummaryEntry], list[DPmVEntry], list[DPmVEntry]]:
    """
    Parse ART summary and generate d-PmV data from atomic resolution projections.

    Processes ARTATOP atomic projection files to create:
    - ART summary with valence/conduction contributions
    - d-PmV and dshg-PmV data for band convergence analysis
    - Energy-dependent d-tensor data files

    Parameters
    ----------
    poscar_path : Path
        Path to POSCAR file for structure information.
    out_dir : Path
        Directory containing ARTATOP atomic projection files.

    Returns
    -------
    tuple[list[ARTSummaryEntry], list[DPmVEntry], list[DPmVEntry], list[Tuple[float, float]], list[Tuple[float, float]]]
        - ART summary entries
        - d-PmV data for band convergence
        - dshg-PmV data for band convergence
        - d-energy data as (energy, value) pairs
        - dshg-energy data as (energy, value) pairs
    """

    structure = Structure.from_file(poscar_path)
    dir_path = poscar_path.parent

    def copy_trim(source_pattern: str, dest_name: str) -> Path:
        """
        Copy file with header removal.

        Parameters
        ----------
        source_pattern : str
            Glob pattern for source files.
        dest_name : str
            Destination filename.

        Returns
        -------
        Path
            Path to created file.
        """
        src = sorted(out_dir.glob(source_pattern))[0]
        dest = out_dir / dest_name
        with src.open() as f_in, dest.open("w") as f_out:
            f_out.writelines(f_in.readlines()[1:])  # skip header
        return dest

    # Create arp-*.dat files
    arp_val_path = copy_trim("arp_shg_val_*.txt", "arp-val.dat")
    arp_con_path = copy_trim("arp_shg_con_*.txt", "arp-con.dat")
    arp_nshg_path = copy_trim("arp_nshg_*.txt", "arp-nshg.dat")
    arp_dshg_path = copy_trim("arp_dshg_*.txt", "arp-dshg.dat")

    # Generate d-PmV.dat from arp-nshg.dat
    d_pmv = []
    with arp_nshg_path.open() as f:
        for line in f:
            parts = line.split()
            nbands = int(float(parts[0]))
            val = float(parts[2]) * (4 / 3 * np.pi * 10 / 2)
            d_pmv.append(DPmVEntry(nbands=nbands, value=val, label="prf-PmV"))

    with open(out_dir / "d-PmV.dat", "w") as f:
        for entry in d_pmv:
            f.write(f"{entry.nbands} {entry.value:.4f}\n")

    # Generate dshg-PmV.dat from arp-dshg.dat
    dshg_pmv = []
    with arp_dshg_path.open() as f:
        for line in f:
            parts = line.split()
            nbands = int(float(parts[0]))
            val = float(parts[2]) * (4 / 3 * np.pi * 10 / 2)
            dshg_pmv.append(DPmVEntry(nbands=nbands, value=val, label="dprf-PmV"))

    with open(out_dir / "dshg-PmV.dat", "w") as f:
        for entry in dshg_pmv:
            f.write(f"{entry.nbands} {entry.value:.6f}\n")

    # Parse summary from original TXT files
    val_lines = arp_val_path.read_text().splitlines()
    con_lines = arp_con_path.read_text().splitlines()
    all_lines = (out_dir / "arp_nonlin.txt").read_text().splitlines()[9:]

    clean_lines = [
        ln
        for ln in all_lines
        if ln.strip()
        and not ln.lstrip().startswith("#")
        and not ln.lstrip().startswith("tot")
    ]

    total_sum = sum(float(line.split()[2]) for line in clean_lines)

    # Check if total_sum is zero to avoid division by zero
    if total_sum == 0:
        logger.info(
            "total_sum is zero, which would cause division by zero. "
            "There are no valid lines in arp_nonlin.txt or all values are zero."
        )
        return []

    atom_order = structure.symbol_set
    amt_dict = structure.composition.get_el_amt_dict()

    atom_type_list = atom_order
    atom_counts = [int(amt_dict[el]) for el in atom_order]

    # Filter val_lines and con_lines the same way as clean_lines
    clean_val_lines = [
        ln
        for ln in val_lines
        if ln.strip()
        and not ln.lstrip().startswith("#")
        and not ln.lstrip().startswith("tot")
    ]

    clean_con_lines = [
        ln
        for ln in con_lines
        if ln.strip()
        and not ln.lstrip().startswith("#")
        and not ln.lstrip().startswith("tot")
    ]

    # CREATE PER-ATOM FILES IN MEMORY
    # Ensure total_atoms is an int for use in range() and list sizing
    total_atoms = int(sum(atom_counts))
    orbitals_per_atom = 3

    # Create per-atom data structures
    atom_total_data = [[] for _ in range(total_atoms)]
    atom_valence_data = [[] for _ in range(total_atoms)]
    atom_conduction_data = [[] for _ in range(total_atoms)]

    # DISTRIBUTE LINES TO PER-ATOM FILES
    # Total data (arp_nonlin.txt)
    total_idx = 0
    for atom_idx in range(total_atoms):
        for orbital in range(orbitals_per_atom):
            if total_idx < len(clean_lines):
                atom_total_data[atom_idx].append(clean_lines[total_idx])
                total_idx += 1
            else:
                atom_total_data[atom_idx].append("0 0 0.0")  # Pad with zeros if missing

    # Valence data (arp-val.dat)
    valence_idx = 0
    for atom_idx in range(total_atoms):
        for orbital in range(orbitals_per_atom):
            if valence_idx < len(clean_val_lines):
                atom_valence_data[atom_idx].append(clean_val_lines[valence_idx])
                valence_idx += 1
            else:
                atom_valence_data[atom_idx].append(
                    "0 0 0.0"
                )  # Pad with zeros if missing

    # Conduction data (arp-con.dat)
    conduction_idx = 0
    for atom_idx in range(total_atoms):
        for orbital in range(orbitals_per_atom):
            if conduction_idx < len(clean_con_lines):
                atom_conduction_data[atom_idx].append(clean_con_lines[conduction_idx])
                conduction_idx += 1
            else:
                atom_conduction_data[atom_idx].append(
                    "0 0 0.0"
                )  # Pad with zeros if missing

    # PROCESS FROM PER-ATOM FILES (like shell script reading temp_$k files)
    summary_entries = []
    global_atom_idx = 0

    for atype, count in zip(atom_type_list, atom_counts):
        occ_s = occ_p = occ_d = vb_s = vb_p = vb_d = cb_s = cb_p = cb_d = 0

        for atom_offset in range(int(count)):
            # Get this atom's data from our in-memory "files"
            atom_idx = global_atom_idx + atom_offset

            # Read from TOTAL "file" for this atom (like temp_$k)
            total_lines = atom_total_data[atom_idx]
            s = float(total_lines[0].split()[2]) if len(total_lines) > 0 else 0.0
            p = float(total_lines[1].split()[2]) if len(total_lines) > 1 else 0.0
            d = float(total_lines[2].split()[2]) if len(total_lines) > 2 else 0.0

            # Read from VALENCE "file" for this atom (like temp_vb_$k)
            valence_lines = atom_valence_data[atom_idx]
            sv = float(valence_lines[0].split()[2]) if len(valence_lines) > 0 else 0.0
            pv = float(valence_lines[1].split()[2]) if len(valence_lines) > 1 else 0.0
            dv = float(valence_lines[2].split()[2]) if len(valence_lines) > 2 else 0.0

            # Read from CONDUCTION "file" for this atom (like temp_cb_$k)
            conduction_lines = atom_conduction_data[atom_idx]
            sc = (
                float(conduction_lines[0].split()[2])
                if len(conduction_lines) > 0
                else 0.0
            )
            pc = (
                float(conduction_lines[1].split()[2])
                if len(conduction_lines) > 1
                else 0.0
            )
            dc = (
                float(conduction_lines[2].split()[2])
                if len(conduction_lines) > 2
                else 0.0
            )

            # Calculate percentages (same as before)
            occ_s += 100 * s / total_sum
            occ_p += 100 * p / total_sum
            occ_d += 100 * d / total_sum
            vb_s += 100 * sv / total_sum
            vb_p += 100 * pv / total_sum
            vb_d += 100 * dv / total_sum
            cb_s += 100 * sc / total_sum
            cb_p += 100 * pc / total_sum
            cb_d += 100 * dc / total_sum

        global_atom_idx += int(count)

        tot = occ_s + occ_p + occ_d
        vb = vb_s + vb_p + vb_d
        cb = cb_s + cb_p + cb_d
        ind = tot / int(count)

        entry = ARTSummaryEntry(
            atom_type=str(atype),
            num_atoms=int(count),
            ind=ind,
            total=tot,
            vb=vb,
            cb=cb,
            vb_s=vb_s,
            vb_p=vb_p,
            vb_d=vb_d,
            cb_s=cb_s,
            cb_p=cb_p,
            cb_d=cb_d,
            tot_s=occ_s,
            tot_p=occ_p,
            tot_d=occ_d,
        )
        summary_entries.append(entry)

    # Write result.art_TOT
    with open(dir_path / "result.art_TOT", "w") as f:
        f.write(
            "Type Natom IND TOT VB CB VB_s VB_p VB_d CB_s CB_p CB_d TOT_s TOT_p TOT_d\n"
        )
        for e in summary_entries:
            f.write(
                f"{e.atom_type} {e.num_atoms} "
                f"{e.ind:.4f} {e.total:.4f} {e.vb:.4f} {e.cb:.4f} "
                f"{e.vb_s:.4f} {e.vb_p:.4f} {e.vb_d:.4f} "
                f"{e.cb_s:.4f} {e.cb_p:.4f} {e.cb_d:.4f} "
                f"{e.tot_s:.4f} {e.tot_p:.4f} {e.tot_d:.4f}\n"
            )

    # Clean up temp files (keep only arp-dshg.dat like shell script)
    for temp in [arp_val_path, arp_con_path, arp_nshg_path]:
        try:
            temp.unlink()
        except Exception as e:
            logger.info(f"Could not remove temporary file {temp}: {e}")

    # --- Generate d_energy.dat and dshg_energy.dat (energy-band.dat must exist) ---
    energy_band = dir_path / "energy-band.dat"
    dpmv_file = out_dir / "d-PmV.dat"
    dshgpmv_file = out_dir / "dshg-PmV.dat"

    if energy_band.exists() and dpmv_file.exists() and dshgpmv_file.exists():

        def paste_and_extract(
            energy_file: Path, data_file: Path, output_file: Path
        ) -> None:
            """
            Combine energy and data files column-wise.

            Parameters
            ----------
            energy_file : Path
                File with energy data.
            data_file : Path
                File with d/dshg data.
            output_file : Path
                Output combined file.
            """
            with open(energy_file) as f1, open(data_file) as f2, open(
                output_file, "w"
            ) as fout:
                for l1, l2 in zip(f1, f2):
                    p1 = l1.split()
                    p2 = l2.split()
                    if len(p1) >= 2 and len(p2) >= 2:
                        # column 2 from energy-band.dat and column 2 from data file
                        fout.write(f"{p1[1]} {p2[1]}\n")

        # Generate d_energy.dat and dshg_energy.dat in main dir
        if energy_band.exists() and dpmv_file.exists() and dshgpmv_file.exists():
            d_energy_out = dir_path / "d_energy.dat"
            dshg_energy_out = dir_path / "dshg_energy.dat"

            paste_and_extract(energy_band, dpmv_file, d_energy_out)
            paste_and_extract(energy_band, dshgpmv_file, dshg_energy_out)

        d_energy = []
        dshg_energy = []

        if d_energy_out.exists():
            with d_energy_out.open() as f:
                for line in f:
                    parts = line.split()
                    if len(parts) >= 2:
                        d_energy.append((float(parts[0]), float(parts[1])))

        if dshg_energy_out.exists():
            with dshg_energy_out.open() as f:
                for line in f:
                    parts = line.split()
                    if len(parts) >= 2:
                        dshg_energy.append((float(parts[0]), float(parts[1])))
    else:
        logger.info(
            "Warning: One or more input files for d_energy.dat generation are missing."
        )

    return summary_entries, d_pmv, dshg_pmv, d_energy, dshg_energy


def parse_procar_and_generate_band_data(
    procar_path: Path, outcar_path: Path, output_file: Path = Path("energy-band.dat")
) -> tuple[Path, int, int, int, float]:
    """
    Parse PROCAR and OUTCAR to generate energy-band data for ART analysis.

    Extracts band energies relative to Fermi level and creates energy-band.dat
    file used for d-energy and dshg-energy analysis.

    Parameters
    ----------
    procar_path : Path
        Path to PROCAR file for band structure data.
    outcar_path : Path
        Path to OUTCAR for Fermi energy.
    output_file : Path
        Output file for energy-band data.

    Returns
    -------
    tuple[Path, int, int, int, float]
        Output file path, number of k-points, bands, ions, and Fermi energy.
    """
    with open(procar_path) as f:
        lines = f.readlines()

    # Extract metadata from line 2
    line2 = lines[1].split()
    nk = int(line2[3])
    nband = int(line2[7])
    nion = int(line2[11])

    # Extract E-fermi
    with open(outcar_path) as f:
        efermi_lines = [line for line in f if "E-fermi" in line]
        efermi = float(efermi_lines[-1].split()[2]) if efermi_lines else 0.0

    # Collect all 'band' lines
    temp1 = [line for line in lines if "band" in line]

    with open(output_file, "w") as fout:
        for i in range(1, nband + 1):
            temp2 = [
                line
                for line in temp1
                if f"band" in line and line.strip().split()[1] == str(i)
            ]

            e_temp = []
            for j in range(nk):
                try:
                    e_j = float(temp2[j].split()[4])
                    e_bandj = e_j - efermi
                    e_temp.append((j + 1, e_bandj))
                except Exception:
                    continue

            if not e_temp:
                continue

            # Extract min, max, and final band energy
            e_vals = [x[1] for x in e_temp]
            e_max = max(e_vals)
            e_min = min(e_vals)
            e_bandj_final = e_temp[-1][1]

            # Decide VB/CB representative energy
            e_prf = e_max if e_bandj_final < 0 else e_min

            # Write with full float precision (8+ digits)
            fout.write(f"{i} {e_prf:.8f} {e_max:.8f} {e_min:.8f}\n")

    return output_file, nk, nband, nion, efermi


def get_artatop_functional_data(
    vasprun_pbe: Vasprun,
    vasprun_hse: Vasprun,
    initial_structure: Structure,
    final_structure: Structure,
    eg_exp: float = 0.0,
) -> dict:
    """
    Extract functional data from VASP calculations for ARTATOP analysis.

    Computes band gaps, scissor corrections, structural parameters, and
    volume-to-bandgap ratios from PBE and HSE calculations.

    Parameters
    ----------
    vasprun_pbe : Vasprun
        PBE calculation results for baseline properties.
    vasprun_hse : Vasprun
        HSE calculation results for accurate band gaps.
    initial_structure : Structure
        Initial structure before relaxation.
    final_structure : Structure
        Final relaxed structure.
    eg_exp : float
        Experimental band gap value (if available).

    Returns
    -------
    dict
        Dictionary containing structural, electronic, and functional properties
        for ARTATOP analysis including band gaps, scissor corrections, and
        volume metrics.
    """
    relax = final_structure.lattice
    natoms = len(final_structure)

    # Bandgap info
    bg_pbe = vasprun_pbe.eigenvalue_band_properties[0]
    bg_hse = vasprun_hse.eigenvalue_band_properties[0]

    eg_pbe = bg_pbe
    eg_hse = bg_hse
    scissor_hse = round(eg_hse - eg_pbe, 3)
    scissor_exp = round(max(0.0, eg_exp - eg_pbe), 3)

    # Volume per atom and V/Eg
    volume_per_atom = relax.volume / natoms
    v_over_eg_hse = round(volume_per_atom / eg_hse, 3) if eg_hse else 0.0
    v_over_eg_exp = round(volume_per_atom / eg_exp, 3) if eg_exp else 0.0

    # Parameters
    aexx = vasprun_hse.parameters.get("AEXX", None)
    nbands = vasprun_hse.parameters.get("NBANDS", None)
    fermi_energy = vasprun_hse.efermi

    return {
        "chemical_formula": final_structure.composition.reduced_formula,
        "space_group": SpacegroupAnalyzer(final_structure).get_space_group_symbol(),
        "point_group": SpacegroupAnalyzer(final_structure).get_point_group_symbol(),
        "n_atoms": natoms,
        "relax_a": round(relax.a, 3),
        "relax_b": round(relax.b, 3),
        "relax_c": round(relax.c, 3),
        "relax_alpha": round(relax.alpha, 3),
        "relax_beta": round(relax.beta, 3),
        "relax_gamma": round(relax.gamma, 3),
        "relax_volume": round(relax.volume, 3),
        "v_over_eg_exp_per_atom": v_over_eg_exp,
        "v_over_eg_hse_per_atom": v_over_eg_hse,
        "bandgap_pbe": round(eg_pbe, 3),
        "bandgap_exp": round(eg_exp, 3),
        "bandgap_hse": round(eg_hse, 3),
        "scissor_exp": scissor_exp,
        "scissor_hse": scissor_hse,
    }


def parse_artatop_outputs(
    dir_name: str,
    input_file: str,
    job_paths: Optional[dict] = None,
    pbe_vasprun_file: Optional[Path] = None,
    hse_vasprun_file: Optional[Path] = None,
    additional_metadata: dict = None,
    energies: list[float] = [0.0, 0.65, 1.167],
) -> ArtatopOutputModel:
    """
    Main function to parse complete ARTATOP outputs and generate analysis files.

    Orchestrates the parsing of all ARTATOP output files including:
    - Optical properties (linear and nonlinear)
    - Atomic resolution tensor (ART) analysis
    - Band structure and k-point information
    - Structural and electronic properties

    Parameters
    ----------
    dir_name : str
        Directory containing ARTATOP output files.
    input_file : str
        ARTATOP input file path.
    job_paths : dict, optional
        Dictionary mapping job types to their directory paths.
    pbe_vasprun_file : Path, optional
        Path to PBE vasprun.xml file.
    hse_vasprun_file : Path, optional
        Path to HSE vasprun.xml file.
    additional_metadata : dict, optional
        Additional metadata including experimental band gap.
    energies : list[float]
        Energy points for analysis (eV).

    Returns
    -------
    ArtatopOutputModel
        Complete parsed ARTATOP output data model.
    """
    dir_path = Path(strip_hostname(dir_name))

    if job_paths is None:
        logger.info("job_paths must be provided to parse_artatop_outputs.")
        return ArtatopOutputModel(dir_name=str(dir_path))

    # Clean all job_paths
    job_paths = {k: strip_hostname(v) for k, v in job_paths.items()}

    # Step 1: Run the full optical response parser
    parsed_data = parse_optical_response(str(dir_path))

    # Step 2: Fix POSCAR if needed and get structure
    poscar_path = dir_path / "POSCAR"
    lines = poscar_path.read_text().splitlines()
    if lines and any(char.isdigit() for char in lines[0]):
        lines[0] = "Generated by pymatgen"
        poscar_path.write_text("\n".join(lines) + "\n")

    relaxed_structure = Structure.from_file(poscar_path)

    # Step 3: Parse band data if files exist
    procar_path = dir_path / "PROCAR"
    outcar_path = dir_path / "OUTCAR"
    if procar_path.exists() and outcar_path.exists():
        parse_procar_and_generate_band_data(procar_path, outcar_path)

    # Step 4: Parse ART data if directory exists
    out_dir = dir_path / "out_nonlin"
    art_summary, d_pmv, dshg_pmv, d_energy, dshg_energy = parse_art_summary_and_dpmv(
        poscar_path, out_dir
    )
    write_result_art_IND(art_summary, filename=dir_path / "result.art_IND")

    # Step 5: Parse atomic orbital contributions
    atomic_contribs = parse_orbital_atomic_contributions(relaxed_structure, out_dir)

    # Step 6: Get functional data from VASP calculations
    if pbe_vasprun_file is not None and hse_vasprun_file is not None:
        pbe_path = Path(pbe_vasprun_file)
        hse_path = Path(hse_vasprun_file)
    else:
        pbe_path = Path(job_paths["static"]) / "vasprun.xml.gz"
        hse_path = Path(job_paths["hse06"]) / "vasprun.xml.gz"

    if not pbe_path.exists():
        logger.info(f"PBE vasprun.xml not found: {pbe_path}")
        return ArtatopOutputModel(dir_name=str(dir_path))
    if not hse_path.exists():
        logger.info(f"HSE vasprun.xml not found: {hse_path}")
        return ArtatopOutputModel(dir_name=str(dir_path))

    vasprun_pbe = Vasprun(str(pbe_path), parse_dos=False)
    vasprun_hse = Vasprun(str(hse_path), parse_dos=False)

    # --- Use PBE relaxed structure ---
    final_structure = relaxed_structure

    # --- Prepare additional metadata ---
    if additional_metadata is None:
        additional_metadata = {}

    # --- Handle eg_exp ---
    eg_exp = additional_metadata.get("eg_exp", 0.0)

    # --- Get functional data from vaspruns ---
    functional_data = get_artatop_functional_data(
        vasprun_pbe=vasprun_pbe,
        vasprun_hse=vasprun_hse,
        initial_structure=vasprun_pbe.initial_structure,
        final_structure=final_structure,
        eg_exp=eg_exp,
    )
    # Define working subdirectories
    relax1_dir = Path(job_paths["relax1"])
    relax2_dir = Path(job_paths["relax2"])
    static_dir = Path(job_paths["static"])
    optics_dir = Path(job_paths["optics"])
    hse06_dir = Path(job_paths["hse06"])

    original_structure = None
    ori_a = ori_b = ori_c = ori_alpha = ori_beta = ori_gamma = None
    ediffg_relax2 = None
    aexx_hse = None
    import gzip
    import xml.etree.ElementTree as ET

    def load_kpoints_from_vasprun(vasprun_path: pathlib.Path) -> Optional[list[int]]:
        """
        Extract the automatic k-point grid (divisions) from a vasprun.xml or .gz file.

        Parameters
        ----------
            vasprun_path : pathlib.Path
            Path to vasprun.xml or vasprun.xml.gz

        Returns
        -------
            Optional[list[int]]
            K-point grid like [12, 12, 12], or None if not found.
        """
        if not vasprun_path.exists():
            logger.info(f"Vasprun file not found: {vasprun_path}")
            return None

        try:
            # Load XML content
            if vasprun_path.suffix == ".gz":
                with gzip.open(vasprun_path, "rt") as f:
                    tree = ET.parse(f)

            else:
                tree = ET.parse(vasprun_path)
            root = tree.getroot()
            for gen in root.iter("generation"):
                if gen.attrib.get("param") in ("Gamma", "kpts"):
                    for v in gen.findall("v"):
                        if v.attrib.get("name") == "divisions":
                            divisions = [int(x) for x in v.text.strip().split()]
                            if len(divisions) >= 3:
                                return divisions[:3]
            logger.info("No <generation> divisions found in vasprun.xml")
            return None

        except Exception as e:
            logger.info(f"Manual parsing of vasprun.xml failed: {e}")
            return None

    if relax1_dir.exists():
        try:
            unrelaxed_poscar = Poscar.from_file(relax1_dir / "POSCAR.gz")
            unrelaxed_structure = unrelaxed_poscar.structure
            lat = unrelaxed_structure.lattice
            original_structure = unrelaxed_structure
            ori_a = lat.a
            ori_b = lat.b
            ori_c = lat.c
            ori_alpha = lat.alpha
            ori_beta = lat.beta
            ori_gamma = lat.gamma
        except Exception as e:
            logger.info(f"Could not read unrelaxed structure: {e}")

    if relax2_dir.exists():
        try:
            ediffg_relax2 = Incar.from_file(relax2_dir / "INCAR.gz").get("EDIFFG", None)
        except Exception as e:
            logger.info(f"Could not read EDIFFG from relax2: {e}")

    if hse06_dir.exists():
        try:
            vasprun_path = hse06_dir / "vasprun.xml.gz"
            vasprun = Vasprun(vasprun_path)
            aexx_hse = vasprun.parameters.get("AEXX", 0.25)
        except Exception as e:
            logger.info(f"Could not read AEXX from HSE calculation: {e}")

    try:
        kpoints_relax1 = load_kpoints_from_vasprun(relax1_dir / "vasprun.xml.gz")
        kpoints_relax2 = load_kpoints_from_vasprun(relax2_dir / "vasprun.xml.gz")
        kpoints_static = load_kpoints_from_vasprun(static_dir / "vasprun.xml.gz")
        kpoints_optics = load_kpoints_from_vasprun(optics_dir / "vasprun.xml.gz")
        kpoints_hse06 = load_kpoints_from_vasprun(hse06_dir / "vasprun.xml.gz")

        # NBANDS and energy
        nbands_static = Vasprun(static_dir / "vasprun.xml.gz").parameters.get("NBANDS")
        nbands_optics = Vasprun(optics_dir / "vasprun.xml.gz").parameters.get("NBANDS")
        total_energy_static = Vasprun(static_dir / "vasprun.xml.gz").final_energy
    except Exception as e:
        logger.info(f"Failed to extract some summary fields: {e}")
        kpoints_relax1 = (
            kpoints_relax2
        ) = kpoints_static = kpoints_optics = kpoints_hse06 = None
        nbands_static = nbands_optics = total_energy_static = None

    # Step 6: Populate the final model
    output_model = parsed_data
    output_model.original_structure = original_structure
    output_model.relaxed_structure = final_structure

    output_model.art_summary = art_summary
    output_model.d_pmV = d_pmv
    output_model.dshg_pmV = dshg_pmv
    output_model.atomic_contributions = atomic_contribs
    output_model.d_energy = d_energy
    output_model.dshg_energy = dshg_energy

    # Add summary and input-extracted properties
    output_model.chemical_formula = functional_data["chemical_formula"]
    output_model.space_group = functional_data["space_group"]
    output_model.point_group = functional_data["point_group"]

    output_model.ori_a = ori_a
    output_model.ori_b = ori_b
    output_model.ori_c = ori_c
    output_model.ori_alpha = ori_alpha
    output_model.ori_beta = ori_beta
    output_model.ori_gamma = ori_gamma

    output_model.n_atoms = functional_data["n_atoms"]
    output_model.relax_a = functional_data["relax_a"]
    output_model.relax_b = functional_data["relax_b"]
    output_model.relax_c = functional_data["relax_c"]
    output_model.relax_alpha = functional_data["relax_alpha"]
    output_model.relax_beta = functional_data["relax_beta"]
    output_model.relax_gamma = functional_data["relax_gamma"]
    output_model.relax_volume = functional_data["relax_volume"]
    output_model.v_over_eg_exp_per_atom = functional_data["v_over_eg_exp_per_atom"]
    output_model.v_over_eg_hse_per_atom = functional_data["v_over_eg_hse_per_atom"]
    output_model.bandgap_pbe = functional_data["bandgap_pbe"]
    output_model.bandgap_exp = functional_data["bandgap_exp"]
    output_model.bandgap_hse = functional_data["bandgap_hse"]
    output_model.scissor_exp = functional_data["scissor_exp"]
    output_model.scissor_hse = functional_data["scissor_hse"]

    def clean_kpoints(k):
        if not k:
            return None
        if isinstance(k[0], list):  # case like [[5,5,5]]
            return k[0]
        if isinstance(k, list) and len(k) == 3:  # case like [5,5,5]
            return k
        return None  # fallback

    # Assign cleaned kpoints to the model
    output_model.kpoints_relax1 = clean_kpoints(kpoints_relax1)
    output_model.kpoints_relax2 = clean_kpoints(kpoints_relax2)
    output_model.kpoints_static = clean_kpoints(kpoints_static)
    output_model.kpoints_optics = clean_kpoints(kpoints_optics)
    output_model.kpoints_hse06 = clean_kpoints(kpoints_hse06)
    output_model.nbands_static = nbands_static
    output_model.nbands_optics = nbands_optics
    output_model.total_energy_static = total_energy_static

    # --- Enthalpy (eV/atom) ---
    natoms = len(final_structure)
    output_model.enthalpy = total_energy_static / natoms

    output_model.ediffg_relax2 = ediffg_relax2
    output_model.aexx_hse = aexx_hse

    cif_name = additional_metadata.get("cif_name", "unknown.cif")
    # Write output files
    write_artatop_summary(output_model, cif_filename=cif_name)
    write_artatop_summary_result_line(output_model, cif_filename=cif_name)
    # Determine highest ART component
    d_tensor_0 = next(
        (x for x in output_model.d_tensor_uv if abs(x.energy - 0.0) < 1e-3), None
    )
    if d_tensor_0:
        highest_label = max(
            d_tensor_0.components, key=lambda k: abs(d_tensor_0.components[k])
        )
        highest_value = abs(d_tensor_0.components[highest_label])
        output_model.art_highest_component = highest_label
        output_model.art_highest_value = highest_value

    write_artatop_analysis_file(
        output_model, cif_filename=cif_name, filename_dir=dir_name
    )

    # Write simplified OPTICS from result.re
    result_re_path = Path(dir_name) / "result.re"
    if result_re_path.exists():
        parse_fixed_lines(result_re_path, cif_name)

    return output_model
