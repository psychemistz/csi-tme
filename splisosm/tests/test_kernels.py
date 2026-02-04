"""Tests for kernel implementations."""

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_less

from splisosm.kernels import (
    ICARKernel,
    CompositionalKernel,
    GaussianKernel,
    LinearKernel,
    clr_transform,
)


class TestLinearKernel:
    def test_compute_shape(self):
        X = np.random.randn(50, 5)
        kernel = LinearKernel()
        K = kernel.compute(X)

        assert K.shape == (50, 50)

    def test_symmetry(self):
        X = np.random.randn(50, 5)
        kernel = LinearKernel()
        K = kernel.compute(X)

        assert_allclose(K, K.T)

    def test_positive_semidefinite(self):
        X = np.random.randn(50, 5)
        kernel = LinearKernel()
        K = kernel.compute(X)

        eigenvalues = np.linalg.eigvalsh(K)
        # All eigenvalues should be >= -1e-10 (numerical tolerance)
        assert np.all(eigenvalues >= -1e-10)


class TestICARKernel:
    @pytest.fixture
    def grid_coords(self):
        # 10x10 grid
        x = np.arange(10)
        y = np.arange(10)
        xx, yy = np.meshgrid(x, y)
        return np.column_stack([xx.ravel(), yy.ravel()])

    def test_compute_shape(self, grid_coords):
        kernel = ICARKernel(k_neighbors=8, rho=0.95)
        L = kernel.compute(grid_coords)

        n = len(grid_coords)
        assert L.shape == (n, n)

    def test_symmetry(self, grid_coords):
        kernel = ICARKernel(k_neighbors=8, rho=0.95)
        L = kernel.compute(grid_coords)

        assert_allclose(L, L.T, rtol=1e-10)

    def test_positive_definite(self, grid_coords):
        kernel = ICARKernel(k_neighbors=8, rho=0.95)
        L = kernel.compute(grid_coords)

        eigenvalues = np.linalg.eigvalsh(L)
        # All eigenvalues should be positive (with regularization)
        assert np.all(eigenvalues > 0)

    def test_rho_validation(self):
        with pytest.raises(ValueError):
            ICARKernel(rho=1.0)
        with pytest.raises(ValueError):
            ICARKernel(rho=0.0)
        with pytest.raises(ValueError):
            ICARKernel(rho=-0.5)

    def test_adjacency_symmetry(self, grid_coords):
        kernel = ICARKernel(k_neighbors=8)
        W = kernel.build_adjacency(grid_coords)

        # Adjacency should be symmetric
        assert_allclose(W.toarray(), W.T.toarray())

    def test_boundary_weighting(self, grid_coords):
        kernel = ICARKernel(k_neighbors=8, rho=0.95)
        L = kernel.compute(grid_coords)

        # Create boundary scores (high at edges)
        n = len(grid_coords)
        boundary_scores = np.ones(n) * 0.1
        boundary_scores[:10] = 1.0  # First row

        L_boundary = kernel.apply_boundary_weighting(boundary_scores)

        # Boundary-weighted kernel should be smaller in interior
        interior_mask = np.arange(n) >= 10
        assert np.mean(L_boundary[np.ix_(interior_mask, interior_mask)]) < \
               np.mean(L[np.ix_(interior_mask, interior_mask)])


class TestCompositionalKernel:
    @pytest.fixture
    def proportions(self):
        # Random proportions that sum to 1
        n, q = 50, 4
        raw = np.random.rand(n, q)
        return raw / raw.sum(axis=1, keepdims=True)

    def test_compute_shape(self, proportions):
        kernel = CompositionalKernel()
        K = kernel.compute(proportions)

        assert K.shape == (50, 50)

    def test_symmetry(self, proportions):
        kernel = CompositionalKernel()
        K = kernel.compute(proportions)

        assert_allclose(K, K.T)

    def test_clr_rows_sum_zero(self, proportions):
        clr_data = clr_transform(proportions)

        # CLR-transformed rows should sum to approximately zero
        row_sums = clr_data.sum(axis=1)
        assert_allclose(row_sums, 0, atol=1e-10)

    def test_clr_handles_zeros(self):
        # Data with some zeros
        proportions = np.array([
            [0.5, 0.3, 0.2, 0.0],
            [0.4, 0.4, 0.1, 0.1],
        ])
        # Should not raise
        clr_data = clr_transform(proportions)
        assert not np.any(np.isnan(clr_data))
        assert not np.any(np.isinf(clr_data))


class TestGaussianKernel:
    def test_compute_shape(self):
        X = np.random.randn(50, 2)
        kernel = GaussianKernel()
        K = kernel.compute(X)

        assert K.shape == (50, 50)

    def test_diagonal_ones(self):
        X = np.random.randn(50, 2)
        kernel = GaussianKernel()
        K = kernel.compute(X)

        # Diagonal should be 1 (K(x, x) = exp(0) = 1)
        assert_allclose(np.diag(K), 1.0)

    def test_values_in_range(self):
        X = np.random.randn(50, 2)
        kernel = GaussianKernel()
        K = kernel.compute(X)

        # All values should be in (0, 1]
        assert np.all(K > 0)
        assert np.all(K <= 1)

    def test_median_heuristic(self):
        X = np.random.randn(50, 2)
        kernel = GaussianKernel(sigma=None)  # Use median heuristic
        K = kernel.compute(X)

        assert kernel.bandwidth is not None
        assert kernel.bandwidth > 0


class TestKernelCentering:
    def test_centering_row_col_means(self):
        X = np.random.randn(50, 5)
        kernel = LinearKernel()
        K = kernel.compute(X)
        K_centered = kernel.center(K)

        # Centered kernel should have zero row/column means
        assert_allclose(K_centered.mean(axis=0), 0, atol=1e-10)
        assert_allclose(K_centered.mean(axis=1), 0, atol=1e-10)

    def test_centering_matrix_identity(self):
        n = 50
        X = np.random.randn(n, 5)
        kernel = LinearKernel()
        K = kernel.compute(X)

        # Manual centering with H matrix
        H = np.eye(n) - np.ones((n, n)) / n
        K_centered_manual = H @ K @ H

        K_centered = kernel.center(K)

        assert_allclose(K_centered, K_centered_manual, rtol=1e-10)
