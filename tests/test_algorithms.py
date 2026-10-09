"""Numerical contracts checked against NumPy, including reported regressions."""

from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from QR import QR, QR_decomposition
from arnoldi import arnoldi_iteration
from bisection import gershgorin_bound, sturm_bisection, sturm_evaluate
from conjugate_gradient import conjugate_gradient
from divide_and_conquer import div_conq, partition
from folded_spectrum import folded_spectrum_iteration
from homotopy import f, homotopy_eigenpairs
from inverse_iter import inverse_iteration
from jacobi import jacobi_eigenvalue_algorithm
from laguerre_iter import construct_characteristic_polynomial, laguerre_method
from lanczos import lanczos
from power_iter import power_iteration
from rayleigh_quotient_iter import rayleigh


ROOT = Path(__file__).resolve().parents[1]


class NumericalTests(unittest.TestCase):
    def setUp(self):
        self.random_state = np.random.get_state()
        np.random.seed(0)
        self.rng = np.random.default_rng(0)
        self.error_state = np.seterr(divide="raise", invalid="raise", over="raise")

    def tearDown(self):
        np.random.set_state(self.random_state)
        np.seterr(**self.error_state)

    def assertEigenvector(self, A, x, tol=1e-8):
        self.assertTrue(np.all(np.isfinite(x)))
        np.testing.assert_allclose(np.linalg.norm(x), 1., atol=1e-12)
        eigval = np.vdot(x, A @ x)
        scale = max(np.linalg.norm(A, ord=np.inf), np.finfo(float).tiny)
        self.assertLessEqual(np.linalg.norm(A @ x - eigval * x), tol * scale)
        return eigval

    def assertEigenbasis(self, A, eigvals, Q, tol=1e-8):
        np.testing.assert_allclose(np.sort(eigvals), np.linalg.eigvalsh(A), atol=tol, rtol=tol)
        np.testing.assert_allclose(Q.conj().T @ Q, np.eye(len(A)), atol=tol)
        np.testing.assert_allclose(A @ Q, Q @ np.diag(eigvals), atol=tol, rtol=tol)

    def test_qr_opposite_equal_magnitude_eigenvalues(self):
        np.testing.assert_allclose(np.sort(QR([[0., 1.], [1., 0.]])), [-1., 1.], atol=1e-12)

    def test_qr_zero_and_repeated_eigenvalues(self):
        for values in ([0., 1.], [0., 0., 0.], [-2., -2., 1., 1.]):
            with self.subTest(values=values):
                Q, _ = np.linalg.qr(self.rng.normal(size=(len(values), len(values))))
                A = Q @ np.diag(values) @ Q.T
                np.testing.assert_allclose(np.sort(QR(A)), np.sort(values), atol=1e-10)

    def test_qr_random_indefinite_and_scaled(self):
        for n in (2, 5, 12):
            M = self.rng.normal(size=(n, n))
            A = M + M.T
            for scale in (1e-100, 1., 1e100):
                with self.subTest(n=n, scale=scale):
                    np.testing.assert_allclose(np.sort(QR(scale * A)) / scale,
                                               np.linalg.eigvalsh(A), atol=1e-9)

    def test_qr_decomposition_rank_deficient_and_complex(self):
        for A in (np.array([[0., 1.], [0., 2.], [0., 3.]]),
                  np.array([[1., 2j, 3.], [2., 1., -1j]])):
            with self.subTest(shape=A.shape):
                Q, R = QR_decomposition(A)
                np.testing.assert_allclose(Q @ R, A, atol=1e-12)
                np.testing.assert_allclose(Q.conj().T @ Q, np.eye(Q.shape[1]), atol=1e-12)

    def test_qr_exhaustion_and_invalid_domain(self):
        with self.assertRaises(RuntimeError):
            QR([[0., 1.], [1., 0.]], max_iter=0)
        with self.assertRaises(ValueError):
            QR([[0., 1.], [0., 0.]])

    def test_arnoldi_real_and_complex_relations(self):
        for complex_input in (False, True):
            A = self.rng.normal(size=(6, 6))
            if complex_input:
                A = A + 1j * self.rng.normal(size=A.shape)
            Q, H = arnoldi_iteration(A, np.ones(6), m=4)
            np.testing.assert_allclose(Q.conj().T @ Q, np.eye(Q.shape[1]), atol=1e-12)
            np.testing.assert_allclose(A @ Q[:, :H.shape[1]], Q @ H, atol=1e-12)

    def test_arnoldi_breakdown_and_full_basis(self):
        Q, H = arnoldi_iteration(np.eye(4), np.ones(4), m=3)
        self.assertEqual(Q.shape, (4, 1))
        self.assertEqual(H.shape, (1, 1))
        A = np.diag([-3., -1., 1., 5.])
        Q, H = arnoldi_iteration(A, np.ones(4), m=4)
        np.testing.assert_allclose(A @ Q, Q @ H, atol=1e-12)
        np.testing.assert_allclose(Q.T @ Q, np.eye(4), atol=1e-12)

    def test_arnoldi_invalid_start_and_steps(self):
        for x0, m in ((np.zeros(3), 2), (np.ones(3), 0), (np.ones(3), 4)):
            with self.subTest(m=m), self.assertRaises(ValueError):
                arnoldi_iteration(np.eye(3), x0, m)

    def test_lanczos_signed_eigenvalues(self):
        A = np.diag([-3., -1., 2.])
        np.testing.assert_allclose(np.linalg.eigvalsh(lanczos(A)), np.linalg.eigvalsh(A), atol=1e-12)

    def test_lanczos_breakdown_and_multiplicities(self):
        np.testing.assert_allclose(lanczos(np.eye(4)), [[1.]])
        T = lanczos(np.diag([1., 1., 2., 2.]))
        self.assertEqual(T.shape, (2, 2))
        np.testing.assert_allclose(np.linalg.eigvalsh(T), [1., 2.], atol=1e-12)

    def test_lanczos_random_real_and_complex(self):
        for complex_input in (False, True):
            M = self.rng.normal(size=(8, 8))
            if complex_input:
                M = M + 1j * self.rng.normal(size=M.shape)
            A = M + M.conj().T
            np.testing.assert_allclose(np.linalg.eigvalsh(lanczos(A)), np.linalg.eigvalsh(A), atol=1e-10)

    def test_lanczos_partial_ritz_values(self):
        A = np.diag([-3., -1., 1., 5.])
        values = np.linalg.eigvalsh(lanczos(A, m=2))
        self.assertEqual(len(values), 2)
        self.assertGreaterEqual(values[0], -3.)
        self.assertLessEqual(values[-1], 5.)

    def test_lanczos_invalid_domain(self):
        with self.assertRaises(ValueError):
            lanczos([[0., 1.], [0., 0.]])
        with self.assertRaises(ValueError):
            lanczos(np.eye(3), x0=np.zeros(3))

    def bisection_values(self, d, e, tol=1e-9):
        lower, upper = gershgorin_bound(d, e)
        return np.array([sturm_bisection(i, d, e, lower, upper, tol=tol)
                         for i in range(1, len(d) + 1)])

    def test_bisection_reducible_matrix(self):
        np.testing.assert_allclose(self.bisection_values([2., 1., 3.], [0., 0.]), [1., 2., 3.], atol=1e-9)

    def test_bisection_scalar_and_repeated_spectrum(self):
        for d in ([2.], [1., 1., 1.], [0., 0., 0.], [-2., 0., -2., 0.]):
            with self.subTest(d=d):
                np.testing.assert_allclose(self.bisection_values(d, np.zeros(len(d) - 1)), np.sort(d), atol=1e-9)

    def test_sturm_zero_pivots_and_equality(self):
        for z, expected in ((-2., 3), (0., 1), (2., 0)):
            self.assertEqual(sturm_evaluate(z, [0., 0., 0.], [1., 1.]), expected)
        self.assertEqual(sturm_evaluate(2., [2., 1., 3.], [0., 0.]), 1)
        self.assertEqual(sturm_evaluate(-np.inf, [1.], []), 1)
        self.assertEqual(sturm_evaluate(np.inf, [1.], []), 0)

    def test_bisection_random_signed_scaled(self):
        for scale in (1e-100, 1., 1e100):
            d, e = self.rng.normal(size=7), self.rng.normal(size=6)
            e[2] = 0
            T = np.diag(d) + np.diag(e, 1) + np.diag(e, -1)
            np.testing.assert_allclose(self.bisection_values(scale * d, scale * e, tol=scale * 1e-10) / scale,
                                       np.linalg.eigvalsh(T), atol=1e-9)

    def test_bisection_already_narrow_bracket(self):
        result = sturm_bisection(1, [1., 2.], [0.], 1., 1. + 1e-8)
        self.assertLess(abs(result - 1.), 1e-6)

    def test_bisection_invalid_bracket_index_and_exhaustion(self):
        for index, lower, upper in ((0, 1., 2.), (3, 1., 2.), (1, 1.5, 2.), (1, 2., 1.)):
            with self.subTest(index=index, lower=lower), self.assertRaises(ValueError):
                sturm_bisection(index, [1., 2.], [0.], lower, upper)
        with self.assertRaises(RuntimeError):
            sturm_bisection(1, [1., 2.], [0.], 1., 2., max_iter=0)

    def test_conjugate_gradient_already_solved(self):
        np.testing.assert_array_equal(conjugate_gradient(np.diag([2., 3.]), [2., 3.], [1., 1.]), [1., 1.])
        np.testing.assert_array_equal(conjugate_gradient(np.eye(3), np.zeros(3), np.zeros(3)), np.zeros(3))

    def test_conjugate_gradient_random_spd_and_input_preservation(self):
        M = self.rng.normal(size=(10, 10))
        A, b, x0 = M.T @ M + np.eye(10), self.rng.normal(size=10), np.zeros(10)
        before = [a.copy() for a in (A, b, x0)]
        x = conjugate_gradient(A, b, x0, tol=1e-10)
        np.testing.assert_allclose(x, np.linalg.solve(A, b), atol=1e-10)
        for actual, expected in zip((A, b, x0), before):
            np.testing.assert_array_equal(actual, expected)

    def test_conjugate_gradient_detected_breakdown_and_exhaustion(self):
        with self.assertRaises(ValueError):
            conjugate_gradient(-np.eye(2), np.ones(2), np.zeros(2))
        with self.assertRaises(RuntimeError):
            conjugate_gradient(np.diag([2., 3.]), np.ones(2), np.zeros(2), max_iter=1)

    def test_divide_conquer_partition_preserves_input(self):
        T = np.array([[2., 1.], [1., 3.]])
        before = T.copy()
        B, update = partition(T)
        np.testing.assert_allclose(B + update, before)
        eigvals, Q = div_conq(T)
        np.testing.assert_array_equal(T, before)
        self.assertEigenbasis(before, eigvals, Q)

    def test_divide_conquer_signed_zero_couplings_and_scalar(self):
        for n in (1, 2, 5, 10):
            d, e = self.rng.normal(size=n), self.rng.normal(size=n - 1)
            if n > 3:
                e[2] = 0
            T = np.diag(d) + np.diag(e, 1) + np.diag(e, -1)
            eigvals, Q = div_conq(T)
            self.assertEqual(eigvals.shape, (n,))
            self.assertEigenbasis(T, eigvals, Q)

    def test_divide_conquer_rejects_nontridiagonal(self):
        with self.assertRaises(ValueError):
            div_conq(np.ones((3, 3)))

    def test_folded_spectrum_original_wrong_target(self):
        A = np.diag([1., 2., 2.9, 4.])
        eigval = self.assertEigenvector(A, folded_spectrum_iteration(A, shift=3.))
        self.assertAlmostEqual(eigval, 2.9, places=9)

    def test_folded_spectrum_exact_target_and_equal_distances(self):
        for A, shift in ((np.diag([1., 3., 5.]), 3.), (np.diag([2., 4.]), 3.), (3. * np.eye(3), 3.)):
            with self.subTest(shift=shift, n=len(A)):
                eigval = self.assertEigenvector(A, folded_spectrum_iteration(A, shift))
                self.assertAlmostEqual(abs(eigval - shift), np.min(abs(np.linalg.eigvalsh(A) - shift)), places=8)

    def test_folded_spectrum_complex_hermitian(self):
        A = np.array([[2., 1j], [-1j, 4.]])
        eigval = self.assertEigenvector(A, folded_spectrum_iteration(A, shift=4.))
        self.assertAlmostEqual(eigval.real, np.linalg.eigvalsh(A)[-1], places=8)

    def test_folded_spectrum_failure_is_explicit(self):
        with self.assertRaises(RuntimeError):
            folded_spectrum_iteration(np.diag([1., 2.]), 1.1, max_iter=0)
        with self.assertRaises(ValueError):
            folded_spectrum_iteration([[0., 1.], [0., 0.]], 1.)

    def test_power_negative_dominant_and_zero_matrix(self):
        A = np.diag([-5., 2., 1.])
        self.assertAlmostEqual(self.assertEigenvector(A, power_iteration(A, 100)), -5., places=10)
        self.assertEigenvector(np.zeros((3, 3)), power_iteration(np.zeros((3, 3)), 10))

    def test_inverse_real_and_complex(self):
        for A, shift in ((np.diag([-3., 1., 4.]), 1.1), (np.array([[2., 1j], [-1j, 4.]]), 1.5)):
            x = inverse_iteration(A, 100, shift)
            eigval = self.assertEigenvector(A, x)
            expected = np.linalg.eigvalsh(A)
            self.assertAlmostEqual(eigval.real, expected[np.argmin(abs(expected - shift))], places=9)

    def test_inverse_singular_shift_and_iteration_validation(self):
        with self.assertRaisesRegex(np.linalg.LinAlgError, "nonsingular"):
            inverse_iteration(np.diag([1., 2.]), 10, shift=1.)
        with self.assertRaises(ValueError):
            power_iteration(np.eye(2), -1)
        with self.assertRaises(ValueError):
            inverse_iteration(np.eye(2), 1, shift=np.nan)

    def test_jacobi_random_and_repeated_spectrum(self):
        for values in ([-3., -1., 2., 4.], [0., 0., 1., 1.]):
            Q, _ = np.linalg.qr(self.rng.normal(size=(4, 4)))
            A = Q @ np.diag(values) @ Q.T
            eigvals, Q, n_iter = jacobi_eigenvalue_algorithm(A)
            self.assertEigenbasis(A, eigvals, Q)
            self.assertGreater(n_iter, 0)

    def test_jacobi_scalar_zero_and_exhaustion(self):
        for A in (np.array([[2.]]), np.zeros((3, 3))):
            eigvals, Q, n_iter = jacobi_eigenvalue_algorithm(A)
            self.assertEigenbasis(A, eigvals, Q)
            self.assertEqual(n_iter, 0)
        with self.assertRaises(RuntimeError):
            jacobi_eigenvalue_algorithm([[0., 1.], [1., 0.]], max_iter=0)

    def test_rayleigh_real_and_complex(self):
        for A, shift, x0 in ((np.diag([-3., 1., 4.]), 1.1, np.array([.1, 1., .2])),
                            (np.array([[2., 1j], [-1j, 4.]]), 1.5, np.array([1., .1j]))):
            self.assertEigenvector(A, rayleigh(A, tol=1e-10, shift=shift, x0=x0))

    def test_rayleigh_already_solved_and_exact_shift(self):
        A = np.diag([1., 2.])
        self.assertEigenvector(A, rayleigh(A, 1e-10, 1., [1., 0.]))
        self.assertEigenvector(A, rayleigh(A, 1e-10, 1., [1., 1.]))

    def test_rayleigh_zero_projection_reports_nonconvergence(self):
        with self.assertRaises(RuntimeError):
            rayleigh(np.diag([-1., 1.]), 1e-10, 0., np.ones(2), max_iter=5)

    def test_characteristic_polynomial_matches_determinant(self):
        for n in (1, 3, 6):
            d, e = self.rng.normal(size=n), self.rng.normal(size=n - 1)
            T = np.diag(d) + np.diag(e, 1) + np.diag(e, -1)
            p = construct_characteristic_polynomial(d, e)
            for z in (-2., .5, 2.):
                np.testing.assert_allclose(p(z), np.linalg.det(T - z * np.eye(n)), atol=1e-10)

    def test_laguerre_array_and_polynomial_inputs(self):
        for p in (np.array([1., 0., -2.]), np.poly1d([1., 0., -2.])):
            z = laguerre_method(p, z0=3., tol=1e-12)
            self.assertAlmostEqual(z, np.sqrt(2.), places=10)

    def test_laguerre_complex_and_stationary_start(self):
        for coefficients in ([1., 0., 1.], [1., -2j, -2.], [1., 0., 0., -1.]):
            z = laguerre_method(coefficients, z0=0., tol=1e-12)
            self.assertLess(abs(np.polyval(coefficients, z)), 1e-10)

    def test_laguerre_scaled_coefficients(self):
        for scale in (1e-120, 1., 1e120):
            z = laguerre_method(scale * np.array([1., 0., -1.]), 2., tol=1e-12)
            self.assertAlmostEqual(z, 1., places=10)

    def test_laguerre_exhaustion_and_invalid_polynomial(self):
        with self.assertRaises(RuntimeError):
            laguerre_method([1., 0., 0., -2.], 3., max_iter=1)
        with self.assertRaises(ValueError):
            laguerre_method([1.], 0.)
        self.assertEqual(laguerre_method([1., -1.], 1., max_iter=0), 1.)

    def test_laguerre_small_tridiagonal_deflation(self):
        d, e = np.array([-3., -1., 1., 2., 4.]), np.array([.2, -.3, .1, .4])
        p = construct_characteristic_polynomial(d, e)
        values = []
        for _ in d:
            z = laguerre_method(p, .5, tol=1e-12)
            values.append(z)
            p = np.polydiv(p, [-1., z])[0]
        T = np.diag(d) + np.diag(e, 1) + np.diag(e, -1)
        np.testing.assert_allclose(np.sort(values), np.linalg.eigvalsh(T), atol=1e-9)

    def test_homotopy_verified_endpoint(self):
        A = np.array([[2., .2, .1], [.2, 4., .3], [.1, .3, 6.]])
        before = A.copy()
        eigvals, Q = homotopy_eigenpairs(A)
        self.assertEigenbasis(A, eigvals, Q, tol=1e-7)
        np.testing.assert_array_equal(A, before)

    def test_homotopy_diagonal_repeated_endpoint(self):
        A = np.diag([3., 1., 1.])
        eigvals, Q = homotopy_eigenpairs(A)
        self.assertEigenbasis(A, eigvals, Q)

    def test_homotopy_rejects_complex_target_and_bad_start(self):
        with self.assertRaises(ValueError):
            homotopy_eigenpairs([[0., -1.], [1., 0.]])
        with self.assertRaises(ValueError):
            homotopy_eigenpairs([[2., .1], [.1, 3.]], D=np.eye(2))

    def test_homotopy_crossing_and_work_limit(self):
        D = np.diag([-1., 1.])
        with self.assertRaisesRegex(RuntimeError, "singular"):
            f(.5, np.array([1., 0., 0.]), -D, D)
        with self.assertRaisesRegex(RuntimeError, "max_iter"):
            homotopy_eigenpairs([[2., .1], [.1, 3.]], max_iter=0)

    def test_homotopy_failed_solver_is_reported(self):
        failed = SimpleNamespace(success=False, y=[], message="test integration failure")
        with patch("homotopy.solve_ivp", return_value=failed):
            with self.assertRaisesRegex(RuntimeError, "test integration failure"):
                homotopy_eigenpairs([[2., .1], [.1, 3.]])

    def test_invalid_dimensions_finite_values_and_tolerances(self):
        for A in (np.ones((2, 3)), np.array([[np.nan]]), np.empty((0, 0))):
            with self.subTest(shape=A.shape), self.assertRaises(ValueError):
                power_iteration(A, 1)
        with self.assertRaises(ValueError):
            gershgorin_bound([1., 2.], [])
        for tol in (0., -1., np.nan, np.inf):
            with self.subTest(tol=tol), self.assertRaises(ValueError):
                QR(np.eye(2), tol=tol)


class DemonstrationTests(unittest.TestCase):
    def test_imports_do_not_run_demonstrations(self):
        modules = [path.stem for path in ROOT.glob("*.py")]
        command = "import " + ", ".join(modules)
        result = subprocess.run([sys.executable, "-B", "-c", command], cwd=ROOT,
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "")

    def test_all_demonstrations_finish_with_finite_output(self):
        for path in ROOT.glob("*.py"):
            if path.stem in ("arnoldi", "_validation"):
                continue
            with self.subTest(module=path.stem):
                result = subprocess.run([sys.executable, "-B", str(path)], cwd=ROOT,
                                        capture_output=True, text=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stderr, "")
                self.assertNotIn("nan", result.stdout.lower())
                self.assertNotIn("inf", result.stdout.lower())


if __name__ == "__main__":
    unittest.main()
