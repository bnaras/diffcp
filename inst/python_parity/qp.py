"""QP-forward parity fixture for the diffcp R package.

Mirrors how cvxpylayers calls `solve_only(P=P)` for forward QP solves
(quadratic objective `0.5 x'Px + c'x`).  Python diffcp's QP support is
forward-only: `solve_and_derivative(P=P, mode='dense'/'lsqr'/'lsmr')`
errors with "Dense, lsqr, and lsmr modes currently do not support
quadratic objectives. Consider switching to 'lpgd' mode."  Only `lpgd`
modes (which we do not port) handle QP derivatives.

Run with:

    uv run python inst/python_parity/qp.py
"""

import numpy as np
import scipy.sparse as sparse

import diffcp.cone_program as cone_prog


def emit(label, vec):
    s = ", ".join(f"{v:.16g}" for v in np.asarray(vec).flatten())
    print(f"## {label}: c({s})")


def emit_csc(label, M):
    M_csc = sparse.csc_matrix(M)
    print(f"## {label}.shape: c({M_csc.shape[0]}, {M_csc.shape[1]})")
    print(f"## {label}.A_i: c({', '.join(str(r + 1) for r in M_csc.indices)})")
    print(f"## {label}.A_p: c({', '.join(str(p + 1) for p in M_csc.indptr)})")
    emit(f"{label}.A_x", M_csc.data)


def main():
    ## Small QP: min 0.5 x'Px + c'x  s.t. x >= 0, sum(x) = 1, n = 3
    n = 3
    P = sparse.csc_matrix(np.array([
        [2.0, 0.5, 0.0],
        [0.5, 2.0, 0.5],
        [0.0, 0.5, 2.0],
    ]))
    c = np.array([-1.0, 0.5, 1.0])
    A_eq = np.ones((1, n))
    A_ineq = -np.eye(n)
    A = sparse.csc_matrix(np.vstack([A_eq, A_ineq]))
    b = np.array([1.0, 0.0, 0.0, 0.0])
    cone_dict = {"z": 1, "l": n}

    print("# === QP forward fixture ===")
    emit_csc("qp.A", A)
    emit_csc("qp.P", P)
    emit("qp.b", b)
    emit("qp.c", c)
    print("## qp.cone_dict: list(z = 1L, l = 3L)")

    for solver in ("Clarabel", "SCS"):
        kwargs = {}
        if solver == "Clarabel":
            kwargs.update(tol_gap_abs=1e-12, tol_gap_rel=1e-12, tol_feas=1e-12)
        else:
            kwargs.update(eps_abs=1e-12, eps_rel=1e-12, max_iters=100000)
        x, y, s = cone_prog.solve_only(A, b, c, cone_dict,
                                       solve_method=solver, P=P, **kwargs)
        tag = solver.lower()
        emit(f"qp.{tag}.x", x)
        emit(f"qp.{tag}.y", y)
        emit(f"qp.{tag}.s", s)


if __name__ == "__main__":
    main()
