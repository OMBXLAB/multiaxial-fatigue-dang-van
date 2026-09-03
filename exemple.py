"""Exécute la démonstration de Dang Van et enregistre son graphique."""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent / "src"))
from dang_van import MateriauDangVan, creer_historique_non_proportionnel, evaluer_dang_van


def main() -> None:
    phase, contraintes = creer_historique_non_proportionnel()
    materiau = MateriauDangVan(220.0, 130.0)
    resultat = evaluer_dang_van(contraintes, materiau)

    print(f"Paramètres de Dang Van : a={materiau.a:.3f}, b={materiau.b:.1f} MPa")
    print(f"Indicateur maximal : {resultat.contrainte_equivalente:.2f} MPa")
    print(f"Coefficient de sécurité : {resultat.coefficient_securite:.2f}")

    dossier_resultats = Path("resultats")
    dossier_resultats.mkdir(exist_ok=True)
    figure, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    axes[0].plot(np.degrees(phase), contraintes[:, 0], label=r"$\sigma_{xx}$")
    axes[0].plot(np.degrees(phase), contraintes[:, 3], label=r"$\tau_{xy}$")
    axes[0].set_ylabel("Contrainte [MPa]")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    hydrostatique = np.sum(contraintes[:, :3], axis=1) / 3.0
    axes[1].plot(np.degrees(phase), hydrostatique, color="tab:green")
    axes[1].set(xlabel="Phase du cycle [degrés]", ylabel="Contrainte hydrostatique [MPa]")
    axes[1].grid(alpha=0.3)
    figure.suptitle(
        "Chargement non proportionnel synthétique — "
        f"coefficient de sécurité = {resultat.coefficient_securite:.2f}"
    )
    figure.tight_layout()
    figure.savefig(dossier_resultats / "chargement_synthetique.png", dpi=180)


if __name__ == "__main__":
    main()
