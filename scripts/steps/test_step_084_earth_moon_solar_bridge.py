"""Controls for the physical-scale matched-profile bridge diagnostic."""
import math
import unittest

import numpy as np
from scipy.optimize import brentq

import step_084_earth_moon_solar_bridge as bridge


class EarthMoonSolarBridgeTests(unittest.TestCase):
    def test_finite_body_surface_continuity_and_central_limit(self):
        for sector in ("baseline", "two_branch"):
            at_surface = bridge.radial_gradient(bridge.GM_EARTH, bridge.EARTH_RADIUS_M,
                                                bridge.EARTH_RADIUS_M, sector)
            inside = bridge.radial_gradient(bridge.GM_EARTH, bridge.EARTH_RADIUS_M,
                                            bridge.EARTH_RADIUS_M * (1 - 1e-6), sector)
            outside = bridge.radial_gradient(bridge.GM_EARTH, bridge.EARTH_RADIUS_M,
                                             bridge.EARTH_RADIUS_M * (1 + 1e-6), sector)
            self.assertAlmostEqual(inside / at_surface, 1, delta=1e-6)
            self.assertAlmostEqual(outside / at_surface, 1, delta=1e-6)
            self.assertEqual(bridge.radial_gradient(bridge.GM_EARTH,
                                                    bridge.EARTH_RADIUS_M, 0, sector), 0)

    def test_transverse_solar_field_prevents_zero_gradient(self):
        for sector in ("baseline", "two_branch"):
            aligned = bridge.bridge(sector, 0)
            transverse = bridge.bridge(sector, 90)
            self.assertLess(aligned["bridge_min_gradient_over_gt"], 1e-9)
            self.assertAlmostEqual(transverse["bridge_min_gradient_over_gt"],
                                   transverse["solar_gradient_over_gt"], places=5)
            self.assertGreaterEqual(aligned["bridge_min_distance_from_earth_km"],
                                    bridge.EARTH_RADIUS_M / 1000)
            self.assertLessEqual(aligned["bridge_min_distance_from_earth_km"],
                                 (bridge.SEPARATION_M - bridge.MOON_RADIUS_M) / 1000)

    def test_narrow_low_stiffness_alignment_window(self):
        solar = bridge.bridge("two_branch", 0)["solar_gradient_over_gt"]
        q_limit = brentq(lambda q: bridge.kinetic_coefficient(q, "two_branch") - 1,
                         0, 1)
        half_width = math.degrees(math.asin(q_limit / solar))
        self.assertLess(half_width, 0.03)
        self.assertGreater(half_width, 0.01)

    def test_flux_law_at_lunar_orbit(self):
        for sector in ("baseline", "two_branch"):
            q = bridge.radial_gradient(bridge.GM_EARTH, bridge.EARTH_RADIUS_M,
                                       bridge.SEPARATION_M, sector)
            flux = bridge.kinetic_coefficient(q, sector) * q
            target = bridge.GM_EARTH / (bridge.G_T * bridge.SEPARATION_M ** 2)
            np.testing.assert_allclose(flux, target, rtol=1e-10)


if __name__ == "__main__":
    unittest.main()
