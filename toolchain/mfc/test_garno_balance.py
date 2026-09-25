"""Independent balance and reactor checks for the uniform Garno regression."""

import ast
import math
import unittest
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parents[2] / "tests"
CASE_IDS = ("C384ADDC", "91AE660D")


class TestGarnoBalance(unittest.TestCase):
    def setUp(self):
        self.case_dir = TESTS_DIR / "C384ADDC"
        self.load_case()

    def load_case(self):
        tree = ast.parse((self.case_dir / "case.py").read_text())
        self.case = ast.literal_eval(next(node.value for node in tree.body if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "case" for t in node.targets)))
        self.golden = {parts[0]: float(parts[1]) for parts in (line.split() for line in (self.case_dir / "golden.txt").read_text().splitlines())}

    def value(self, index, step):
        return self.golden[f"D/cons.{index}.00.{step:06d}.dat"]

    def test_material_and_energy_balance(self):
        for case_id in CASE_IDS:
            with self.subTest(case_id=case_id):
                self.case_dir = TESTS_DIR / case_id
                self.load_case()
                final = self.case["t_step_stop"]
                initial = [self.value(i, 0) for i in range(1, 9)]
                result = [self.value(i, final) for i in range(1, 9)]
                self.assertTrue(all(math.isfinite(value) for value in result))
                self.assertAlmostEqual(result[0], initial[0], delta=1e-10)  # air mass
                self.assertAlmostEqual(result[1] + result[2], initial[1] + initial[2], delta=1e-8)
                self.assertAlmostEqual(result[5], initial[5], delta=1e-12)  # air volume
                self.assertAlmostEqual(result[6] + result[7], initial[6] + initial[7], delta=1e-12)
                self.assertAlmostEqual(result[4] - initial[4], self.case["rburn%q"] * (initial[1] - result[1]), delta=1e-3)

    def test_reactant_fraction_matches_independent_ode(self):
        for case_id in CASE_IDS:
            with self.subTest(case_id=case_id):
                self.case_dir = TESTS_DIR / case_id
                self.load_case()
                case = self.case
                rho = sum(self.value(i, 0) for i in range(1, 4))
                y = self.value(2, 0) / rho
                dt = case["dt"] * case["t_step_stop"] / 40000

                def derivative(value):
                    ratio = rho / case["rburn%rho0"]
                    ignition = case["rburn%ki"] * value ** case["rburn%m1"] * abs(1 - ratio) ** case["rburn%m2"]
                    growth = case["rburn%kg"] * value ** case["rburn%n1"] * (1 - value) ** case["rburn%n2"] * ratio ** case["rburn%n3"]
                    return -(ignition + growth)

                for _ in range(40000):
                    k1 = derivative(y)
                    k2 = derivative(y + 0.5 * dt * k1)
                    k3 = derivative(y + 0.5 * dt * k2)
                    k4 = derivative(y + dt * k3)
                    y += dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6

                self.assertAlmostEqual(self.value(2, case["t_step_stop"]) / rho, y, delta=1e-9)


if __name__ == "__main__":
    unittest.main()
