# atomate2-artatop
Atomate2 extension integrating the (Atomic Response Theory Toolkit for Optical Properties) [ARTATOP](http://www.artatop.net/) for optical properties analysis, including full atomic and orbital-resolved SHG contributions powered by atomate2, jobflow, VASP, and ARTATOP. Built and tested with Python 3.10+.

Badges
-----------------
[![python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/downloads/)

Overview
--------

Taking advantage of Atomate2’s systematic workflow framework, ARTATOP is integrated into [Atomate2](https://github.com/materialsproject/atomate2) to automate and streamline the calculation of nonlinear optical properties of materials, enabling atom-resolved analysis of second-harmonic generation (SHG) responses using advanced theoretical approaches. Workflow execution is managed remotely via [jobflow-remote](https://github.com/Matgenix/jobflow-remote), providing robust and scalable job handling. 
Developed within the Atomate2 framework, the ARTATOP workflow supports fully automated, high-throughput calculations of linear and nonlinear optical properties in crystalline materials. Its purpose is to ensure reproducible and efficient computation of optical tensors, including dielectric constants, refractive indices, birefringence, and SHG coefficients, based on VASP and atomic response theory.

Installation
------------
```bash
git clone https://github.com/xiyue-cheng/atomate2-artatop.git
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
    lattice=[[0, a / 2, a / 2],[a / 2, 0, a / 2],[a / 2, a / 2, 0]],
    species=["Si", "C"],
    coords=[[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
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
    lattice=[[0, a / 2, a / 2],[a / 2, 0, a / 2],[a / 2, a / 2, 0]],
    species=["Si", "C"],
    coords=[[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
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
    worker="mars_worker",
    resources={"nodes": 1, "ntasks": 64, "partition": "mars"},
)
```

Workflow stages
---------------
1. Relaxation: double relaxation to optimize the structure.
2. Static: wavefunction and charge density and bandgap.
3. Optic: uses dense k-points and increased bands for frequency-dependent optical response.
4. HSE06: refined bandgap for scissor correction.
5. ARTATOP: linear and nonlinear responses optical response( dielectric function, refractive indices, birefringence, SHG coefficient) with atomic and orbital contribution analysis.
6. Results: structured ARTATOP output files and JSON summary.

Outputs and plotting
--------------------
 Example plotting for text outputs:
```python
from jobflow import SETTINGS
from atomate2_artatop.plotting import plot_d_energy

store = SETTINGS.JOB_STORE
store.connect()

result = store.query_one({"name": "artatop"}, load=True)

plot_d_energy(result)

```

Tips and notes
--------------
- Set `ARTATOP_CMD` and `VASP_CMD` to cluster-accessible executables.
- Align job store and Mongo configurations across submit and query environments.
- Tune `resources` for your partitions.
- HSE06 improves the bandgap used for scissor corrections; ensure sufficient wall time.

License and citation
--------------------
The project is an add-on / extension. It uses the same license as Atomate2. The license text is available in the repository [LICENSE](LICENSE).
If you use atomate2-artatop, please cite the following repository:

```@software{artatop_workflow,
  author = {Sher, Asma and Cheng, Xiyue and Deng, Shuiquan},
  title  = {ARTATOP Workflow for Optical and Nonlinear Optical Property Calculations},
  year   = {2025},
  note   = {Manuscript in preparation},
  url    = {https://github.com/xiyue-cheng/atomate2-artatop}
}
```
