# atomate2-artatop
Atomate2 extension integrating the Atomic Response Theory Toolkit (ARTATOP) for optical properties analysis, including full atomic and orbital-resolved SHG contributions powered by atomate2, jobflow, VASP, and ARTATOP. Built and tested with Python 3.10+.

Badges (optional)
-----------------
[![python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/downloads/)

Overview
--------
Taking advantage of Atomate2’s systematic approach, we have integrated ARTATOP (Atomic Response Theory Toolkit for Optical Properties) into Atomate2 to automate and streamline the calculation of nonlinear optical properties of materials. ARTATOP enables the analysis of each atom’s contribution to the second-harmonic generation (SHG) response using advanced theoretical approaches. Workflow execution is managed remotely via jobflow-remote for robust and scalable job handling.

Features
--------
The ARTATOP workflow is developed within the Atomate2 framework to provide fully automated, high-throughput calculations of both linear and nonlinear optical properties in crystalline materials. Its key purpose is to ensure reproducible and efficient computation of optical tensors—including dielectric constants, refractive indices, birefringence, and second-harmonic generation (SHG) coefficients—based on VASP and atomic response theory.

Installation
------------
```bash
pip install -r requirements.txt
pip install .
```

Post-install configuration
--------------------------
Set your executables in `artatop.yaml`:
```yaml
ARTATOP_CMD: <<ARTATOP_CMD>>
```
Please add the following line to your `.bashrc` or `.bash_profile` to define
the environment variable pointing to the ARTATOP configuration file:

```bash
export ARTATOP_CONFIG_FILE="/path/to/artatop.yaml"

```

Quick start (local run)
-----------------------
```python
from jobflow import run_locally
from pymatgen.core import Structure
from atomate2_artatop.flows.artatop import ArtatopWorkflowMaker

structure = Structure(
    lattice=[[0, 2.13, 2.13], [2.13, 0, 2.13], [2.13, 2.13, 0]],
    species=["Mg", "O"],
    coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
)
eg_exp = 0.0

workflow_maker = ArtatopWorkflowMaker()
flow = workflow_maker.make(
    structure,
    prev_dir=None,
    additional_metadata={"eg_exp": eg_exp, "cif_name": "example.cif"},
)
run_locally(flow, create_folders=True)
```

Quick start (jobflow-remote)
----------------------------
```python
from jobflow_remote import submit_flow, set_run_config
from pymatgen.core import Structure
from atomate2_artatop.flows.artatop import ArtatopWorkflowMaker

structure = Structure(
    lattice=[[0, 2.13, 2.13], [2.13, 0, 2.13], [2.13, 2.13, 0]],
    species=["Mg", "O"],
    coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
)

flow = ArtatopWorkflowMaker().make(
    structure=structure,
    prev_dir=None,
    additional_metadata={"cif_name": "example.cif"},
)

flow = set_run_config(
    flow,
    name_filter="optics",
    resources={"nodes": 1, "ntasks": 1, "partition": ",mars"},
)

submit_flow(
    flow,
    project="atlas",
    worker="xy_worker",
    resources={"nodes": 1, "ntasks": 64, "partition": "mars"},
)
```

Workflow stages
---------------
1. Geometry: double relaxation to optimize the structure.
2. Static preconverge: wavefunction and charge density.
3. Optics NSCF: dense k-points and increased bands for frequency-dependent dielectric response.
4. HSE06 static: refined band structure and bandgap for scissor corrections.
5. ARTATOP computation: linear and nonlinear responses (dielectric, refractive indices, birefringence, d_eff) with atomic and orbital analysis.
6. Results packaging: aggregates ARTATOP outputs into structured files and a JSON summary.

Outputs and plotting
--------------------
Key ARTATOP outputs include optical tensors, chi2, d_eff, and atomic/orbital SHG contributions. Example plotting for text outputs:
```python
import numpy as np
import matplotlib.pyplot as plt
import os

def load_data(filename):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"{filename} not found.")
    data = np.loadtxt(filename)
    return data[:, 0], data[:, 1]

E_v, zeta_v = load_data("nshgv")
E_c, zeta_c = load_data("nshgc")
E_dv, dzeta_v = load_data("dshgv")
E_dc, dzeta_c = load_data("dshgc")

plt.figure(figsize=(8, 3.5))
plt.plot(E_v, zeta_v, color="darkgreen", label="zetaV(E_B) (VB)")
plt.plot(E_c, zeta_c, color="mediumpurple", label="zetaC(E_B) (CB)")
plt.axvline(0, color="gray", linestyle="--", linewidth=0.8)
plt.xlabel("E_B (eV)")
plt.ylabel("zeta(E_B)")
plt.title("zeta(E_B) vs E_B")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig("zeta_combined.png", dpi=300)
plt.close()

plt.figure(figsize=(8, 3.5))
plt.plot(E_dv, dzeta_v, color="darkgreen", label="delta zetaV(E_B) (VB)")
plt.plot(E_dc, dzeta_c, color="mediumpurple", label="delta zetaC(E_B) (CB)")
plt.axvline(0, color="gray", linestyle="--", linewidth=0.8)
plt.xlabel("E_B (eV)")
plt.ylabel("delta zeta(E_B)")
plt.title("delta zeta(E_B) vs E_B")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig("dzeta_combined.png", dpi=300)
plt.close()
```

Tips and notes
--------------
- Set `ARTATOP_CMD` and `VASP_CMD` to cluster-accessible executables.
- Align job store and Mongo configurations across submit and query environments.
- Tune `resources` for your partitions; optics steps may need dense meshes and more bands.
- HSE06 improves the bandgap used for scissor corrections; ensure sufficient wall time.

License and citation
--------------------
The project is an add-on / extension. It uses the same license as Atomate2. The license text is available in the repository [LICENSE](LICENSE).
