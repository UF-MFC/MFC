"""Independent balance and reactor checks for the uniform I&G regression."""

import math
import re
import unittest
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parents[2] / "tests"
CASES = {"5FE634A0": 1.0, "3E6116D6": 2.0}
RHO0 = 1800.0
Q = 4.0e6
KI = 1.0e9
KG = 1.0e8
M1 = 1.0
M2 = 4.0
N1 = 1.0
N3 = 2.0
DT = 1.0e-9


class TestIgnitionGrowthBalance(unittest.TestCase):
    def load_golden(self, case_id):
        path = TESTS_DIR / case_id / "golden.txt"
        return {parts[0]: float(parts[1]) for parts in (line.split() for line in path.read_text().splitlines())}

    def value(self, golden, index, step):
        return golden[f"D/cons.{index}.00.{step:06d}.dat"]

    def test_material_and_energy_balance(self):
        for case_id in CASES:
            with self.subTest(case_id=case_id):
                golden = self.load_golden(case_id)
                final = max(int(step) for step in re.findall(r"D/cons\.1\.00\.(\d+)\.dat", "\n".join(golden)))
                initial = [self.value(golden, i, 0) for i in range(1, 9)]
                result = [self.value(golden, i, final) for i in range(1, 9)]
                self.assertTrue(all(math.isfinite(value) for value in result))
                self.assertAlmostEqual(result[0], initial[0], delta=1e-10)  # air mass
                self.assertAlmostEqual(result[1] + result[2], initial[1] + initial[2], delta=1e-8)
                self.assertAlmostEqual(result[5], initial[5], delta=1e-12)  # air volume
                self.assertAlmostEqual(result[6] + result[7], initial[6] + initial[7], delta=1e-12)
                self.assertAlmostEqual(result[4] - initial[4], Q * (initial[1] - result[1]), delta=1e-3)

    def test_reactant_fraction_matches_independent_ode(self):
        for case_id, n2 in CASES.items():
            with self.subTest(case_id=case_id):
                golden = self.load_golden(case_id)
                final = max(int(step) for step in re.findall(r"D/cons\.1\.00\.(\d+)\.dat", "\n".join(golden)))
                rho = sum(self.value(golden, i, 0) for i in range(1, 4))
                y = self.value(golden, 2, 0) / rho
                dt = DT * final / 40000

                def derivative(value):
                    ratio = rho / RHO0
                    ignition = KI * value**M1 * abs(1 - ratio) ** M2
                    growth = KG * value**N1 * (1 - value) ** n2 * ratio**N3
                    return -(ignition + growth)

                for _ in range(40000):
                    k1 = derivative(y)
                    k2 = derivative(y + 0.5 * dt * k1)
                    k3 = derivative(y + 0.5 * dt * k2)
                    k4 = derivative(y + dt * k3)
                    y += dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6

                self.assertAlmostEqual(self.value(golden, 2, final) / rho, y, delta=1e-9)


if __name__ == "__main__":
    unittest.main()
