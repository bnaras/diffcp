"""Derivative parity fixtures (mode='dense').

Run from the project root with:

    uv run python inst/python_parity/derivative.py

Emits R-ready ``c(...)`` literals for forward (D) and adjoint (DT)
derivatives applied to fixed perturbations, on the same fixtures used
by `inst/python_parity/forward_solve.py`.

Pinning these in `tests/testthat/test-derivative.R` ensures that the
R port of `solve_and_derivative` matches Python diffcp at the level of
a single (M, dQ, pi_z, dz) round-trip.
"""

import numpy as np
import scipy.sparse as sparse

import diffcp.cone_program as cone_prog


def emit(label, vec):
    s = ", ".join(f"{v:.12g}" for v in np.asarray(vec).flatten())
    print(f"## {label}: c({s})")


def fixture_lp():
    n = 3
    c = np.array([1.0, 2.0, 3.0])
    A_eq = np.ones((1, n))
    A_l = -np.eye(n)
    A = sparse.csc_matrix(np.vstack([A_eq, A_l]))
    b = np.array([1.0, 0.0, 0.0, 0.0])
    cone_dict = {"z": 1, "l": n}
    return A, b, c, cone_dict


def fixture_soc():
    n = 3
    nv = 1 + n
    c = np.zeros(nv); c[0] = 1.0
    A_eq = np.zeros((1, nv));  A_eq[0, 1:] = 1.0
    A_soc = -np.eye(nv)
    A = sparse.csc_matrix(np.vstack([A_eq, A_soc]))
    b = np.zeros(1 + nv); b[0] = 1.0
    cone_dict = {"z": 1, "q": [nv]}
    return A, b, c, cone_dict


def run(name, fixture, tol=1e-10):
    A, b, c, cone_dict = fixture()
    x, y, s, D, DT = cone_prog.solve_and_derivative(
        A, b, c, cone_dict,
        solve_method="Clarabel",
        mode="dense",
        tol_gap_abs=tol, tol_gap_rel=tol, tol_feas=tol)
    print(f"# === Problem {name} ===")
    emit(f"{name}.x", x); emit(f"{name}.y", y); emit(f"{name}.s", s)

    # Forward derivative perturbations: deterministic, easy to pin.
    np.random.seed(42)
    dA = sparse.csc_matrix(
        (np.random.randn(A.nnz) * 1e-3, A.nonzero()),
        shape=A.shape,
    )
    db = np.arange(1, b.size + 1) * 1e-3
    dc = np.arange(1, c.size + 1) * 1e-3

    dx, dy, ds = D(dA, db, dc)
    ## dA was built from CSC constructor; .data is already CSC-order.
    emit(f"{name}.dA_data", dA.tocsc().data)
    emit(f"{name}.db", db)
    emit(f"{name}.dc", dc)
    emit(f"{name}.dx", dx)
    emit(f"{name}.dy", dy)
    emit(f"{name}.ds", ds)

    # Adjoint perturbations.
    np.random.seed(43)
    dxa = np.random.randn(c.size)
    dya = np.random.randn(b.size)
    dsa = np.random.randn(b.size)
    dAa, dba, dca = DT(dxa, dya, dsa)
    emit(f"{name}.dxa", dxa); emit(f"{name}.dya", dya); emit(f"{name}.dsa", dsa)
    ## Emit dAa values in CSC (column-major) storage order so R's
    ## `Matrix::summary()` round-trip matches without reordering.
    dAa_csc = dAa.tocsc()
    emit(f"{name}.dAa_data", dAa_csc.data)
    emit(f"{name}.dba", dba)
    emit(f"{name}.dca", dca)
    ## Likewise for the forward-derivative input dA: emit in CSC order
    ## so the R test can reconstruct without ambiguity.
    print(f"## {name}.dA_csc_indptr_indices: indptr={list(dA.tocsc().indptr)} indices={list(dA.tocsc().indices)}")


if __name__ == "__main__":
    run("lp", fixture_lp)
    run("soc", fixture_soc)
