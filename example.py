"""Run the synthetic Dang Van demonstration and save its result chart."""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent / "src"))
from dang_van import DangVanMaterial, evaluate_dang_van, synthetic_non_proportional_history


def main() -> None:
    phase, stress = synthetic_non_proportional_history()
    material = DangVanMaterial(reversed_bending_limit=220.0, reversed_torsion_limit=130.0)
    result = evaluate_dang_van(stress, material)

    print(f"Dang Van parameters: a={material.a:.3f}, b={material.b:.1f} MPa")
    print(f"Maximum indicator: {result.equivalent_stress:.2f} MPa")
    print(f"Safety factor: {result.safety_factor:.2f}")

    output = Path("results")
    output.mkdir(exist_ok=True)
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    axes[0].plot(np.degrees(phase), stress[:, 0], label=r"$\sigma_{xx}$")
    axes[0].plot(np.degrees(phase), stress[:, 3], label=r"$\tau_{xy}$")
    axes[0].set_ylabel("Stress [MPa]")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    hydrostatic = np.sum(stress[:, :3], axis=1) / 3.0
    axes[1].plot(np.degrees(phase), hydrostatic, color="tab:green")
    axes[1].set(xlabel="Cycle phase [deg]", ylabel="Hydrostatic stress [MPa]")
    axes[1].grid(alpha=0.3)
    fig.suptitle(f"Synthetic non-proportional loading — safety factor = {result.safety_factor:.2f}")
    fig.tight_layout()
    fig.savefig(output / "synthetic_loading.png", dpi=180)


if __name__ == "__main__":
    main()
