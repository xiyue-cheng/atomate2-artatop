"""Input file generation for ARTATOP calculations."""
from pathlib import Path


class InputFileHandler:
    """
    A handler for generating ARTATOP input files: input_lin, input_nlin, and input_art.
    """

    def get_input_set(
        self,
        calc_type: str,
        calc_dir: Path,
        component: str = None,
        scissor: float = 0.00,
    ) -> Path:
        """
        Generate the appropriate input file for the given calculation type.

        Parameters
        ----------
        calc_type : str
            Type of calculation (e.g., "lin", "nlin", "art").
        calc_dir : Path
            Directory where the input file will be created.
        component : str, optional
            ART component code (e.g., "111", "122", "144") representing Cartesian directions.
            Required for ART calculations.
        scissor : float, optional
            Scissor correction in eV, by default 0.00.

        Returns
        -------
        Path
            Path to the generated input file.
        """
        calc_dir.mkdir(parents=True, exist_ok=True)

        if calc_type == "lin":
            content = """LO 77
$dft_src lvasp=T $end
$opc maxomega = 15 domega = 0.01167  scissor={scissor:.3f}  ecutmin = 0.03  smear = 0.03 $end
""".format(
                scissor=scissor
            )
            (calc_dir / "out_lin").mkdir(exist_ok=True)
        elif calc_type == "nlin":
            content = """NO 777
$dft_src lvasp=T $end
$opc maxomega = 1.167  domega = 0.01167  scissor={scissor:.3f}  ecutmin = 0.03  smear = 0.03 $end
""".format(
                scissor=scissor
            )
            (calc_dir / "out_nonlin").mkdir(exist_ok=True)
        elif calc_type == "art":
            if not component:
                raise ValueError("Component is required for ART calculations.")
            content = f"""AR {component}
$dft_src lvasp=T $end
$opc maxomega = 0  domega = 0.01167  scissor={scissor:.3f}  ecutmin = 0.03  smear = 0.03 $end
""".format(
                scissor=scissor
            )
        else:
            raise ValueError(f"Unsupported calculation type: {calc_type}")

        input_file = calc_dir / f"input_{calc_type}"
        with open(input_file, "w") as f:
            f.write(content)

        return input_file

    @staticmethod
    def determine_highest_component(result_re_file: Path) -> str:
        """
        Determine the highest component from any valid non-zero T-matrix in result.re.
        Tries all "d at omega = ..." blocks (0 eV, 0.65 eV, 1.167 eV, etc).
        """
        if not result_re_file.exists():
            raise FileNotFoundError(f"{result_re_file} not found.")

        matrix_to_component = {
            (1, 1): "111",
            (1, 2): "122",
            (1, 3): "144",
            (1, 4): "124",
            (1, 5): "114",
            (1, 6): "112",
            (2, 1): "211",
            (2, 2): "222",
            (2, 3): "244",
            (2, 4): "224",
            (2, 5): "214",
            (2, 6): "212",
            (3, 1): "411",
            (3, 2): "422",
            (3, 3): "444",
            (3, 4): "424",
            (3, 5): "414",
            (3, 6): "412",
        }

        with open(result_re_file) as f:
            lines = f.readlines()

        max_value = None
        max_position = None

        for idx, line in enumerate(lines):
            if line.strip().lower().startswith("d at omega"):
                try:
                    matrix = []
                    for i in range(1, 4):
                        row_line = lines[idx + i].strip()
                        row = [float(x) for x in row_line.split()]
                        if len(row) != 6:
                            raise ValueError(f"Expected 6 values, got {len(row)}")
                        matrix.append(row)

                    if all(all(val == 0.0 for val in row) for row in matrix):
                        continue

                    for row_idx, row in enumerate(matrix):
                        for col_idx, value in enumerate(row):
                            if max_value is None or abs(value) > abs(max_value):
                                max_value = value
                                max_position = (row_idx + 1, col_idx + 1)

                except Exception as e:
                    print(f"Skipping block at line {idx} due to error: {e}")
                    continue

        if max_position and max_position in matrix_to_component:
            return matrix_to_component[max_position]

        raise ValueError("No valid non-zero T-matrix block found in result.re.")
