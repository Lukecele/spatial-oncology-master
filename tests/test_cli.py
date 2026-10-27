#!/usr/bin/env python3
"""
Unit tests for the spatial-sanctuary CLI entry point.
"""

import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class TestSpatialSanctuaryCLI(unittest.TestCase):
    def test_cli_version(self):
        from spatial_sanctuary import __version__
        cmd = [sys.executable, "-m", "spatial_sanctuary.cli", "--version"]
        res = subprocess.run(cmd, cwd=str(ROOT), env={"PYTHONPATH": str(ROOT / "src")}, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertIn(f"spatial-sanctuary {__version__}", res.stdout)

    def test_cli_simulate(self):
        cmd = [sys.executable, "-m", "spatial_sanctuary.cli", "simulate", "--json"]
        res = subprocess.run(cmd, cwd=str(ROOT), env={"PYTHONPATH": str(ROOT / "src")}, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertIn('"sanctuary_exists": true', res.stdout)
        self.assertIn('"critical_sanctuary_radius_r_crit_um": 165.9', res.stdout)

    def test_cli_lead_time(self):
        cmd = [sys.executable, "-m", "spatial_sanctuary.cli", "lead-time", "--json"]
        res = subprocess.run(cmd, cwd=str(ROOT), env={"PYTHONPATH": str(ROOT / "src")}, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertIn('"lead_time_months": 5.4', res.stdout)

    def test_cli_validate(self):
        cmd = [sys.executable, "-m", "spatial_sanctuary.cli", "validate", "--json"]
        res = subprocess.run(cmd, cwd=str(ROOT), env={"PYTHONPATH": str(ROOT / "src")}, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertIn('"model1_extreme_contrast"', res.stdout)


if __name__ == "__main__":
    unittest.main()
