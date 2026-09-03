import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))
from dang_van import DangVanMaterial, evaluate_dang_van, hydrostatic_stress, voigt_to_tensor


class DangVanTests(unittest.TestCase):
    def test_material_calibration(self) -> None:
        material = DangVanMaterial(220.0, 130.0)
        self.assertAlmostEqual(material.a, 3.0 * 130.0 / 220.0 - 1.5)
        self.assertEqual(material.b, 130.0)

    def test_voigt_conversion_and_hydrostatic_stress(self) -> None:
        history = np.array([[90.0, 60.0, 30.0, 10.0, 20.0, 15.0]])
        tensors = voigt_to_tensor(history)
        self.assertEqual(tensors.shape, (1, 3, 3))
        self.assertTrue(np.allclose(tensors[0], tensors[0].T))
        self.assertAlmostEqual(hydrostatic_stress(tensors)[0], 60.0)

    def test_zero_history_has_infinite_safety_factor(self) -> None:
        result = evaluate_dang_van(
            np.zeros((20, 6)), DangVanMaterial(220.0, 130.0), 3, 4, 4
        )
        self.assertAlmostEqual(result.equivalent_stress, 0.0)
        self.assertTrue(np.isinf(result.safety_factor))

    def test_invalid_history_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            voigt_to_tensor(np.zeros((10, 5)))


if __name__ == "__main__":
    unittest.main()
