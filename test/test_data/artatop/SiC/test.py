from atomate2_artatop.artatop_parser import parse_artatop_outputs

case = "/home/asma/atomate2/tests/test_data/artatop/SiC/output"
job_paths = {
  "relax1": "/home/asma/atomate2/tests/test_data/artatop/SiC/jobs/relax1",
  "relax2": "/home/asma/atomate2/tests/test_data/artatop/SiC/jobs/relax2",
  "static": "/home/asma/atomate2/tests/test_data/artatop/SiC/jobs/static",
  "optics": "/home/asma/atomate2/tests/test_data/artatop/SiC/jobs/optics",
  "hse06":  "/home/asma/atomate2/tests/test_data/artatop/SiC/jobs/hse06",
}
out = parse_artatop_outputs(
  dir_name=case,
  input_file="unused",
  job_paths=job_paths,
  additional_metadata={"cif_name": "SiC.cif", "eg_exp": 1.17},
)
