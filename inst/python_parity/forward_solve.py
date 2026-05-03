"""Forward-solve parity fixtures for the diffcp R package.

Run from the project root with:

    uv run python inst/python_parity/forward_solve.py

Emits R-ready ``c(...)`` literals for the expected `(x, y, s)` of a few
small cone programs, which `tests/testthat/test-forward-solve.R` pins
against.

Each fixture also dumps the problem data (`A`, `b`, `c`, `cone_dict`)
in a form the R test can reconstruct verbatim, so the parity bar is
bit-exact (modulo solver tolerance).
"""

import numpy as np
import scipy.sparse as sparse

import diffcp.cone_program as cone_prog


def emit(label, vec):
    s = ", ".join(f"{v:.12g}" for v in np.asarray(vec).flatten())
    print(f"## {label}: c({s})")


def emit_problem(name, A, b, c, cone_dict):
    """Emit R code that reconstructs (A, b, c, cone_dict)."""
    print(f"# === Problem {name} ===")
    print(f"## A.shape = {A.shape}")
    rows, cols = A.nonzero()
    vals = np.asarray(A.todense()).flatten("F")
    A_dense = np.asarray(A.todense())
    print(f"# A as a dense matrix (column-major flatten):")
    flat = ", ".join(f"{v:.12g}" for v in A_dense.flatten("F"))
    print(f"# A_data <- c({flat})")
    print(f"# A <- matrix(A_data, nrow = {A.shape[0]}, ncol = {A.shape[1]})")
    bs = ", ".join(f"{v:.12g}" for v in b)
    cs = ", ".join(f"{v:.12g}" for v in c)
    print(f"# b <- c({bs})")
    print(f"# c <- c({cs})")
    print(f"# cone_dict <- {cone_dict!r}")


def fixture_lp_eq_constrained():
    """LP with one equality constraint and a nonneg variable.

    minimise c^T x
    s.t.     1^T x = 1,  x >= 0
    in standard form (A x + s = b, s in K = {0} x R_+):
        A = [ 1 1 1 ;
              -I_3 ]
        b = [ 1, 0, 0, 0 ]
        cone_dict = {z: 1, l: 3}
    """
    n = 3
    c = np.array([1.0, 2.0, 3.0])
    A_eq = np.ones((1, n))
    A_l = -np.eye(n)
    A = sparse.csc_matrix(np.vstack([A_eq, A_l]))
    b = np.array([1.0, 0.0, 0.0, 0.0])
    cone_dict = {"z": 1, "l": n}
    return A, b, c, cone_dict


def fixture_soc_program():
    """SOC: minimise t s.t. ||x||_2 <= t, sum(x) = 1, n=3.

    Variables ordered (t, x1, x2, x3).  Only c_t = 1 in the objective.

    Standard form (A z + s = b, z = (t, x), s in K = {0} x SOC^4):
        eq: 0*t + sum(x) = 1     (z cone)
        soc: (t; x) in SOC       => -I @ z + s = 0, s in SOC^4
    """
    n = 3
    nv = 1 + n  # (t, x)
    c = np.zeros(nv); c[0] = 1.0

    A_eq = np.zeros((1, nv));  A_eq[0, 1:] = 1.0
    A_soc = -np.eye(nv)
    A = sparse.csc_matrix(np.vstack([A_eq, A_soc]))
    b = np.zeros(1 + nv); b[0] = 1.0
    cone_dict = {"z": 1, "q": [nv]}
    return A, b, c, cone_dict


def fixture_least_squares():
    """A small least-1-norm problem with simplex constraints.

    Mirrors `diffcp.utils.least_squares_eq_scs_data(m=4, n=2, seed=0)`
    in shape but with hand-pinned (A, b, c) so the R test does not
    depend on CVXPY canonicalisation.
    """
    np.random.seed(0)
    m, n = 4, 2
    Ad = np.random.randn(m, n)
    bd = np.random.randn(m)

    # Variables: x (n), t_i (m) for |Ax-b|_1 = sum t_i, with
    # -t <= Ax-b <= t, x >= 0, sum(x) = 1.
    nv = n + m
    c = np.concatenate([np.zeros(n), np.ones(m)])

    rows = []
    # eq: sum(x) = 1
    eq = np.concatenate([np.ones(n), np.zeros(m)])
    rows.append((eq, 1.0, "z"))
    # x >= 0  ->  -x + s = 0, s in R_+
    for i in range(n):
        r = np.zeros(nv); r[i] = -1.0
        rows.append((r, 0.0, "l"))
    # t_i + (Ax-b)_i >= 0  ->  -(t_i + A_i x - b_i) + s = 0
    for i in range(m):
        r = np.zeros(nv); r[:n] = -Ad[i]; r[n + i] = -1.0
        rows.append((r, -bd[i], "l"))
    # t_i - (Ax-b)_i >= 0  ->  -(t_i - A_i x + b_i) + s = 0
    for i in range(m):
        r = np.zeros(nv); r[:n] = Ad[i]; r[n + i] = -1.0
        rows.append((r, bd[i], "l"))

    A = sparse.csc_matrix(np.array([r[0] for r in rows]))
    b = np.array([r[1] for r in rows])
    cone_dict = {"z": 1, "l": n + 2 * m}
    return A, b, c, cone_dict


def run(name, fixture, tol=1e-10):
    A, b, c, cone_dict = fixture()
    x, y, s = cone_prog.solve_only(
        A, b, c, cone_dict,
        solve_method="Clarabel",
        tol_gap_abs=tol, tol_gap_rel=tol, tol_feas=tol)
    emit_problem(name, A, b, c, cone_dict)
    emit(f"{name}.x", x)
    emit(f"{name}.y", y)
    emit(f"{name}.s", s)


if __name__ == "__main__":
    run("lp_eq", fixture_lp_eq_constrained)
    run("soc", fixture_soc_program)
    run("ls", fixture_least_squares)
