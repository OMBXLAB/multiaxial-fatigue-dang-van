import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))
from dang_van import MateriauDangVan, contrainte_hydrostatique, evaluer_dang_van, voigt_vers_tenseur


class TestsDangVan(unittest.TestCase):
    """Vérifications unitaires des principales fonctions du démonstrateur."""

    def test_calibrage_materiau(self) -> None:
        materiau = MateriauDangVan(220.0, 130.0)
        self.assertAlmostEqual(materiau.a, 3.0 * 130.0 / 220.0 - 1.5)
        self.assertEqual(materiau.b, 130.0)

    def test_conversion_voigt_et_contrainte_hydrostatique(self) -> None:
        historique = np.array([[90.0, 60.0, 30.0, 10.0, 20.0, 15.0]])
        tenseurs = voigt_vers_tenseur(historique)
        self.assertEqual(tenseurs.shape, (1, 3, 3))
        self.assertTrue(np.allclose(tenseurs[0], tenseurs[0].T))
        self.assertAlmostEqual(contrainte_hydrostatique(tenseurs)[0], 60.0)

    def test_historique_nul_coefficient_infini(self) -> None:
        resultat = evaluer_dang_van(
            np.zeros((20, 6)), MateriauDangVan(220.0, 130.0), 3, 4, 4
        )
        self.assertAlmostEqual(resultat.contrainte_equivalente, 0.0)
        self.assertTrue(np.isinf(resultat.coefficient_securite))

    def test_format_historique_invalide(self) -> None:
        with self.assertRaises(ValueError):
            voigt_vers_tenseur(np.zeros((10, 5)))


if __name__ == "__main__":
    unittest.main()
