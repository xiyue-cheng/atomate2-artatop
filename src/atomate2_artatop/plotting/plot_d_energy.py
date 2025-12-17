import pandas as pd
import matplotlib.pyplot as plt


def plot_d_energy(result, savepath=None):
    output_data = result["output"]["output_data"]
    d_energy_data = output_data["d_energy"]

    df = pd.DataFrame(d_energy_data)
    df_vb = df[df["E_repr"] < 0]
    df_cb = df[df["E_repr"] >= 0]

    plt.figure(figsize=(8, 4))

    plt.plot(df_vb["E_repr"], df_vb["value"], label="VB")
    plt.plot(df_cb["E_repr"], df_cb["value"], label="CB")

    plt.scatter(df_vb["E_repr"], df_vb["value"], s=10, alpha=0.6)
    plt.scatter(df_cb["E_repr"], df_cb["value"], s=10, alpha=0.6)

    plt.axvline(0, linestyle=":", linewidth=1)

    plt.xlabel("E_B (eV)")
    plt.ylabel(r"$\zeta(E_B)$ (pm/V)")
    plt.title("Energy-Dependent Nonlinear Response")

    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.3)
    plt.tight_layout()

    if savepath:
        plt.savefig(savepath, dpi=300)
    else:
        plt.show()
