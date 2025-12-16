from atomate2_artatop.schemas import ArtatopOutputModel
from pathlib import Path
from pymatgen.io.vasp import Vasprun


def parse_fixed_lines(input_file: str, cif_file_path: str) -> None:
    """
    Parse ARTATOP output file and extract optical properties.

    Extracts dielectric constants, refractive indices, birefringence,
    effective nonlinear coefficients, and d-tensor components.

    Parameters
    ----------
    input_file : str
        Path to ARTATOP output file to parse.
    cif_file_path : str
        Path to CIF file for material name extraction.

    Returns
    -------
    None
        Writes output to OPTIC-{material}.dat file.
    """
    output_file = f"OPTIC-{Path(cif_file_path).stem}.dat"
    material_name = Path(cif_file_path).stem

    with open(input_file, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    # Data containers
    e_xx = e_yy = e_zz = None
    n_data = {}
    biref_data = {}
    deff_data = {}
    d_tensor_data = {}

    idx = 0
    while idx < len(lines):
        line = lines[idx]

        if "linear optic epsilon at 0.000 eV" in line:
            parts = lines[idx + 1].split()
            e_xx, e_yy, e_zz = parts[0], parts[1], parts[2]
            idx += 2

        elif "refractive n at 0.000 eV" in line:
            n_data["0.000"] = lines[idx + 1].split()
            idx += 2

        elif "birefringence at 0.000 eV" in line:
            biref_data["0.000"] = lines[idx + 1].strip()
            idx += 2

        elif "linear optic epsilon at 0.650 eV" in line:
            idx += 2  # skip (not needed)

        elif "refractive n at 0.650 eV" in line:
            n_data["0.650"] = lines[idx + 1].split()
            idx += 2

        elif "birefringence at 0.650 eV" in line:
            biref_data["0.650"] = lines[idx + 1].strip()
            idx += 2

        elif "linear optic epsilon at 1.167 eV" in line:
            idx += 2  # skip (not needed)

        elif "refractive n at 1.167 eV" in line:
            n_data["1.167"] = lines[idx + 1].split()
            idx += 2

        elif "birefringence at 1.167 eV" in line:
            biref_data["1.167"] = lines[idx + 1].strip()
            idx += 2

        elif "deff at omega = 0.000 eV" in line:
            deff_data["0.000"] = lines[idx + 1].strip()
            idx += 2

        elif "deff at omega = 0.650 eV" in line:
            deff_data["0.650"] = lines[idx + 1].strip()
            idx += 2

        elif "deff at omega = 1.167 eV" in line:
            deff_data["1.167"] = lines[idx + 1].strip()
            idx += 2

        elif "d at omega = 0.000 eV" in line:
            d_tensor_data["0.000"] = [
                lines[idx + 1].split(),
                lines[idx + 2].split(),
                lines[idx + 3].split(),
            ]
            idx += 4

        elif "d at omega = 0.650 eV" in line:
            d_tensor_data["0.650"] = [
                lines[idx + 1].split(),
                lines[idx + 2].split(),
                lines[idx + 3].split(),
            ]
            idx += 4

        elif "d at omega = 1.167 eV" in line:
            d_tensor_data["1.167"] = [
                lines[idx + 1].split(),
                lines[idx + 2].split(),
                lines[idx + 3].split(),
            ]
            idx += 4

        else:
            idx += 1

    with open(output_file, "w") as f:
        f.write(f"{material_name}\n")
        f.write("#optic properties for optic_EgHSE\n")
        f.write("dielectric functions e_xx e_yy e_zz\n")
        f.write(f"0eV      {e_xx} {e_yy} {e_zz}\n")

        f.write("refractive n and birefringence\n")
        for e in ["0.000", "0.650", "1.167"]:
            n = n_data.get(e, ["0.0", "0.0", "0.0"])
            biref = biref_data.get(e, "0.0")
            energy_label = (
                e.replace("0.000", "0eV")
                .replace("0.650", "0.65eV")
                .replace("1.167", "1.167eV")
            )
            f.write(f"{energy_label}   {n[0]} {n[1]} {n[2]} {biref}\n")

        f.write(
            f"deff at 0eV 0.65eV 1.167eV {deff_data.get('0.000', '0.0')} {deff_data.get('0.650', '0.0')} {deff_data.get('1.167', '0.0')}\n"
        )

        d_energy_mapping = {"0.000": "0 eV", "0.650": "0.65 eV", "1.167": "1.167 eV"}

        for e in ["0.000", "0.650", "1.167"]:
            if e in d_tensor_data:
                f.write(f"d at omega = {d_energy_mapping[e]}\n")
                for row in d_tensor_data[e]:
                    f.write(" ".join(row) + "\n")


def write_artatop_summary(
    output_model, cif_filename: str, filename_dir: str = "./"
) -> None:
    """
    Write comprehensive ARTATOP summary to OUTPUT-{material}.dat.

    Includes structural parameters, band gaps, optical properties,
    and calculation parameters.

    Parameters
    ----------
    output_model : ArtatopOutputModel
        Parsed ARTATOP output data.
    cif_filename : str
        CIF filename for material identification.
    filename_dir : str
        Output directory for summary file.

    Returns
    -------
    None
        Writes summary to OUTPUT-{material}.dat.
    """
    filename = Path(filename_dir) / f"OUTPUT-{Path(cif_filename).stem}.dat"

    with open(filename, "w") as f:
        # --- Cell Parameters ---
        f.write("#Cell parameters\n")
        f.write(f"Chemical fomula                {Path(cif_filename).stem}\n")

        spacegroup = output_model.space_group or "-"
        pointgroup = output_model.point_group or "-"
        natoms = output_model.n_atoms or 0

        f.write(f"Space group                    {spacegroup:<30}\n")
        f.write(f"Point group                    {pointgroup:<30}\n")
        f.write(f"Natom                          {natoms:<30}\n")

        ori = output_model
        f.write(f"ori a (A)                      {ori.ori_a or 0.0:.3f}\n")
        f.write(f"ori b (A)                      {ori.ori_b or 0.0:.3f}\n")
        f.write(f"ori c (A)                      {ori.ori_c or 0.0:.3f}\n")
        f.write(f"ori alpha (degree)             {ori.ori_alpha or 0.0:.3f}\n")
        f.write(f"ori beta (degree)              {ori.ori_beta or 0.0:.3f}\n")
        f.write(f"ori gama (degree)              {ori.ori_gamma or 0.0:.3f}\n")

        f.write(f"relax a (A)                    {output_model.relax_a or 0.0:.3f}\n")
        f.write(f"relax b (A)                    {output_model.relax_b or 0.0:.3f}\n")
        f.write(f"relax c (A)                    {output_model.relax_c or 0.0:.3f}\n")
        f.write(
            f"relax alpha (degree)           {output_model.relax_alpha or 0.0:.3f}\n"
        )
        f.write(
            f"relax beta (degree)            {output_model.relax_beta or 0.0:.3f}\n"
        )
        f.write(
            f"relax gama (degree)            {output_model.relax_gamma or 0.0:.3f}\n"
        )
        f.write(
            f"relax Volume (A^3)             {output_model.relax_volume or 0.0:.3f}\n"
        )

        f.write(
            f"V/Eg-EXP/Natom (A^3/eV)        {output_model.v_over_eg_exp_per_atom or 0.0:.3f}\n"
        )
        f.write(
            f"V/Eg-HSE/Natom (A^3/eV)        {output_model.v_over_eg_hse_per_atom or 0.0:.3f}\n"
        )

        f.write(
            f"Eg-PBE (eV)                    {output_model.bandgap_pbe or 0.0:.3f}\n"
        )
        f.write(
            f"Eg-EXP (eV)                    {output_model.bandgap_exp or 0.0:.3f}\n"
        )
        aexx_val = (
            float(output_model.aexx_hse)
            if isinstance(output_model.aexx_hse, float)
            else 0.25
        )
        f.write(
            f"Eg-HSE (eV) AEXX={aexx_val:.4f}        {output_model.bandgap_hse or 0.0:.3f}\n"
        )
        f.write(f"AEXX use for HSE               {aexx_val:.4f}\n")
        f.write(
            f"Scissor-EXP (eV)               {output_model.scissor_exp or 0.0:.3f}\n"
        )
        f.write(
            f"Scissor-HSE (eV)               {output_model.scissor_hse or 0.0:.3f}\n"
        )

        f.write(" \n\n")

        f.write("#LO: optic_EgHSE\n")
        for energy in [0.0, 0.65, 1.167]:
            lin_list = [
                l
                for l in (output_model.linear_response_uv or [])
                if abs(l.energy - energy) < 1e-3
            ]
            if not lin_list:
                lin_list = [
                    l
                    for l in (output_model.linear_response_ir or [])
                    if abs(l.energy - energy) < 1e-3
                ]

            if lin_list:
                if abs(energy) < 1e-3:
                    f.write(
                        f"e_xx at 0 eV                   {lin_list[0].real_part:.3f}\n"
                    )
                    f.write(
                        f"e_yy at 0 eV                   {lin_list[1].real_part:.3f}\n"
                    )
                    f.write(
                        f"e_zz at 0 eV                   {lin_list[2].real_part:.3f}\n"
                    )
                else:
                    f.write(
                        f"n_xx at {energy:.3f} eV               {lin_list[0].refractive_index:.3f}\n"
                    )
                    f.write(
                        f"n_yy at {energy:.3f} eV               {lin_list[1].refractive_index:.3f}\n"
                    )
                    f.write(
                        f"n_zz at {energy:.3f} eV               {lin_list[2].refractive_index:.3f}\n"
                    )

            biref_list = (output_model.birefringence_uv or []) + (
                output_model.birefringence_ir or []
            )
            biref_val = next(
                (b for b in biref_list if abs(b.energy - energy) < 1e-3), None
            )
            if biref_val:
                f.write(
                    f"delta_n at {energy:.3f} eV             {biref_val.delta_n:.3f}\n"
                )

        f.write("\n\n")

        f.write("#NLO: optic_EgHSE\n")

        for energy in [0.0, 0.65, 1.167]:
            all_deff = (output_model.deff_values_uv or []) + (
                output_model.deff_values_ir or []
            )
            deff_val = next(
                (d for d in all_deff if abs(d.energy - energy) < 1e-3), None
            )
            if deff_val:
                f.write(f"d_eff at {energy:.3f} eV              {deff_val.deff:.3f}\n")

        # Write d_11 to d_36 only for 0 eV
        for d in output_model.d_tensor_uv or []:
            if abs(d.energy) < 1e-3:
                component_keys = [
                    "d_xxx",
                    "d_xyy",
                    "d_xzz",
                    "d_xyz",
                    "d_xxz",
                    "d_xxy",
                    "d_yxx",
                    "d_yyy",
                    "d_yzz",
                    "d_yyz",
                    "d_yxz",
                    "d_yxy",
                    "d_zxx",
                    "d_zyy",
                    "d_zzz",
                    "d_zyz",
                    "d_zxz",
                    "d_zxy",
                ]
                voigt_labels = [
                    "d_11",
                    "d_12",
                    "d_13",
                    "d_14",
                    "d_15",
                    "d_16",
                    "d_21",
                    "d_22",
                    "d_23",
                    "d_24",
                    "d_25",
                    "d_26",
                    "d_31",
                    "d_32",
                    "d_33",
                    "d_34",
                    "d_35",
                    "d_36",
                ]
                for key, label in zip(component_keys, voigt_labels):
                    val = d.components.get(key, 0.0)
                    f.write(f"{label} at 0 eV                   {val:.3f}\n")

        f.write("#Elastic constants in GPa\n")
        f.write("none\n")
        f.write("#Elastic properties in GPa\n")
        f.write("none\n")

        f.write("#Calculation parameters\n")

        k1 = output_model.kpoints_relax1 or [[1, 1, 1]]
        k4 = output_model.kpoints_relax2 or [[1, 1, 1]]
        k2 = output_model.kpoints_static or [[1, 1, 1]]
        k3 = output_model.kpoints_optics or [[1, 1, 1]]

        def kformat(klist: list[int]) -> str:
            return f"{klist[0]} x {klist[1]} x {klist[2]}"

        f.write(f"K-mesh for relax               {kformat(k1)}\n")
        f.write(f"K-mesh for static              {kformat(k2)}\n")
        f.write(f"K-mesh for optic               {kformat(k3)}\n")
        f.write(f"K-mesh for HSE                 {kformat(k4)}\n")

        aexx_val = float(output_model.aexx_hse) if output_model.aexx_hse else 0.25
        f.write(f"AEXX use for HSE               {aexx_val:.4f}\n")

        ediffg_value = (
            output_model.ediffg_relax2
            if output_model.ediffg_relax2 is not None
            else -0.000
        )
        f.write(f"EDIFFG                         {ediffg_value:.3f}\n")

        # Flags
        f.write("relax convergence problem                                      \n")
        f.write("HSE convergence problem                                        \n")

        total_energy = output_model.total_energy_static or 0.0
        f.write(f"Total energy from static (eV)  {total_energy:.3f}\n")
        f.write(
            f"NBANDS for static              {output_model.nbands_static if output_model.nbands_static else 0}\n"
        )
        f.write(
            f"NBANDS for optic               {output_model.nbands_optics if output_model.nbands_optics else 0}\n"
        )


def write_artatop_summary_result_line(
    output_model, cif_filename: str, filename_dir: str = "./"
) -> None:
    """
    Write condensed result line to RESULT-{material}.dat.

    Creates a single-line summary with key optical and structural properties
    for high-throughput analysis.

    Parameters
    ----------
    output_model : ArtatopOutputModel
        Parsed ARTATOP output data.
    cif_filename : str
        CIF filename for material identification.
    filename_dir : str
        Output directory for result file.

    Returns
    -------
    None
        Writes condensed results to RESULT-{material}.dat.
    """
    filename = Path(filename_dir) / f"RESULT-{Path(cif_filename).stem}.dat"
    material_id = Path(cif_filename).stem

    cell = output_model
    natoms = output_model.n_atoms or 0
    sg = output_model.space_group or "-"
    point_type = getattr(output_model, "point_group_type", "NCS")  # fallback

    total_energy = output_model.total_energy_static or 0.0
    enthalpy = total_energy / natoms if natoms else 0.0

    eg_pbe = output_model.bandgap_pbe or 0.0
    eg_hse = output_model.bandgap_hse or 0.0

    deff_uv_list = output_model.deff_values_uv or []
    deff_ir_list = output_model.deff_values_ir or []

    def get_deff(e):
        val = next(
            (d.deff for d in (deff_uv_list + deff_ir_list) if abs(d.energy - e) < 1e-3),
            0.0,
        )
        return val

    deff_0 = get_deff(0.0)
    deff_ir = get_deff(0.65)
    deff_uv = get_deff(1.167)

    # delta_n
    biref_list = (output_model.birefringence_uv or []) + (
        output_model.birefringence_ir or []
    )

    def get_dn(e):
        return next((b.delta_n for b in biref_list if abs(b.energy - e) < 1e-3), 0.0)

    dn_ir = get_dn(0.65)
    dn_uv = get_dn(1.167)

    # Write line
    with open(filename, "w") as f:
        f.write(
            "System                   Natom   SG_relax       NCS/CS    enthalpy  Eg_PBE    Eg_type   Type    deff_0    deff_IR   deff_UV   dn_IR     dn_UV\n"
        )
        f.write(
            "                                                          (eV/atom)  (eV)      (eV)              (pm/V)    (pm/V)    (pm/V)              \n"
        )

        f.write(
            f"{material_id:<24}{natoms:<8}{sg:<15}{point_type:<10}"
            f"{enthalpy:9.3f}{eg_pbe:9.2f}{eg_hse:10.2f}{'HSE':>8}"
            f"{deff_0:10.3f}{deff_ir:10.3f}{deff_uv:10.3f}{dn_ir:10.3f}{dn_uv:10.3f}\n"
        )


from pathlib import Path


def write_artatop_analysis_file(
    output_model, cif_filename: str, filename_dir: str = "./"
) -> None:
    """
    Write atomic resolution analysis to ART-{material}.dat.

    Generates detailed atomic contribution analysis for the largest
    SHG tensor component.

    Parameters
    ----------
    output_model : ArtatopOutputModel
        Parsed ARTATOP output data.
    cif_filename : str
        CIF filename for material identification.
    filename_dir : str
        Output directory for analysis file.

    Returns
    -------
    None
        Writes atomic analysis to ART-{material}.dat.
    """

    cif_stem = Path(cif_filename).stem
    filename = Path(filename_dir) / f"ART-{cif_stem}.dat"

    summary_entries = output_model.art_summary or []

    highest_label = output_model.art_highest_component or "d33"
    highest_value = (
        output_model.art_highest_value
        if output_model.art_highest_value is not None
        else 0.0
    )

    component_to_voigt = {
        "d_xxx": "d_11",
        "d_xyy": "d_12",
        "d_xzz": "d_13",
        "d_xyz": "d_14",
        "d_xxz": "d_15",
        "d_xxy": "d_16",
        "d_yxx": "d_21",
        "d_yyy": "d_22",
        "d_yzz": "d_23",
        "d_yyz": "d_24",
        "d_yxz": "d_25",
        "d_yxy": "d_26",
        "d_zxx": "d_31",
        "d_zyy": "d_32",
        "d_zzz": "d_33",
        "d_zyz": "d_34",
        "d_zxz": "d_35",
        "d_zxy": "d_36",
    }
    voigt_label = component_to_voigt.get(highest_label, highest_label)

    with open(filename, "w") as f:
        f.write(f"{cif_stem} {voigt_label} {abs(highest_value):.3f}\n")
        f.write(
            f"ART analysis results for the largest component {voigt_label} for optic_EgHSE\n"
        )
        f.write(
            "Type    Natom   IND     TOT     VB      CB      VB_s    VB_p    VB_d    CB_s    CB_p    CB_d    TOT_s   TOT_p   TOT_d   Atao\n"
        )

        for e in summary_entries:
            atao = (e.vb_s + e.vb_p + e.vb_d + e.cb_s + e.cb_p + e.cb_d) / e.num_atoms
            f.write(
                f"{e.atom_type:<8}{e.num_atoms:<8}"
                f"{e.ind:<8.2f}{e.total:<8.2f}{e.vb / e.num_atoms:<8.2f}{e.cb / e.num_atoms:<8.2f}"
                f"{e.vb_s / e.num_atoms:<8.2f}{e.vb_p / e.num_atoms:<8.2f}{e.vb_d / e.num_atoms:<8.2f}"
                f"{e.cb_s / e.num_atoms:<8.2f}{e.cb_p / e.num_atoms:<8.2f}{e.cb_d / e.num_atoms:<8.2f}"
                f"{e.tot_s / e.num_atoms:<8.2f}{e.tot_p / e.num_atoms:<8.2f}{e.tot_d / e.num_atoms:<8.2f}"
                f"{atao:<8.2f}\n"
            )
