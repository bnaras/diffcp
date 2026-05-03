"""EXP cone parity fixture, mirroring
  diffcp/tests/test_clarabel.py::test_expcone

Builds the entropy-maximisation problem `min -sum(entr(y))` s.t.
`sum(y) = 1` for n=10, canonicalises to SCS form via
`scs_data_from_cvxpy_problem`, and emits R-ready CSC triplets and a
deterministic finite-difference perturbation (with the same random
seed Python uses) so the R test can pin against Python's
solve_and_derivative + finite-difference assertion.

Run with:

    uv run python inst/python_parity/expcone.py
"""

import numpy as np
import scipy.sparse as sparse
import cvxpy as cp


def scs_data_from_cvxpy_problem(problem):
    data = problem.get_problem_data(cp.SCS)[0]
    cone_dims = cp.reductions.solvers.conic_solvers.scs_conif.dims_to_solver_dict(
        data["dims"]
    )
    return data["A"], data["b"], data["c"], cone_dims


def emit(label, vec):
    s = ", ".join(f"{v:.16g}" for v in np.asarray(vec).flatten())
    print(f"## {label}: c({s})")


def emit_csc(label, A):
    A = sparse.csc_matrix(A)
    print(f"## {label}.shape: c({A.shape[0]}, {A.shape[1]})")
    print(f"## {label}.A_i: c({', '.join(str(r + 1) for r in A.indices)})")
    print(f"## {label}.A_p: c({', '.join(str(p + 1) for p in A.indptr)})")
    emit(f"{label}.A_x", A.data)


def emit_cone_dict(label, cd):
    parts = []
    for k, v in cd.items():
        if isinstance(v, (list, tuple)):
            if len(v) == 0:
                parts.append(f"{k} = c()")
            else:
                parts.append(f"{k} = c({', '.join(str(int(x)) + 'L' for x in v)})")
        else:
            parts.append(f"{k} = {int(v)}L")
    print(f"## {label}: list({', '.join(parts)})")


def main():
    np.random.seed(0)
    n = 10
    y = cp.Variable(n)
    obj = cp.Minimize(- cp.sum(cp.entr(y)))
    const = [cp.sum(y) == 1]
    prob = cp.Problem(obj, const)
    A, b, c, cone_dict = scs_data_from_cvxpy_problem(prob)

    print(f"# === expcone fixture (n={n}) ===")
    emit_csc("expcone", A)
    emit("expcone.b", b)
    emit("expcone.c", c)
    emit_cone_dict("expcone.cone_dict", cone_dict)

    import diffcp.cone_program as cone_prog
    x, yv, s, D, DT = cone_prog.solve_and_derivative(
        sparse.csc_matrix(A), b, c, cone_dict,
        solve_method="Clarabel",
        mode="dense",
        tol_gap_abs=1e-13, tol_gap_rel=1e-13,
        tol_feas=1e-13, tol_ktratio=1e-13)
    emit("expcone.x", x)
    emit("expcone.y", yv)
    emit("expcone.s", s)

    ## Same finite-difference perturbations Python uses (np.random.seed(0)
    ## was set above and consumed by cvxpy; we re-seed for the perturbation
    ## block separately to make the R test reproducible).
    np.random.seed(0)
    dA = sparse.csc_matrix(
        (np.random.normal(0, 1e-6, size=A.nnz), A.nonzero()),
        shape=A.shape,
    )
    db = np.random.normal(0, 1e-6, size=b.size)
    dc = np.random.normal(0, 1e-6, size=c.size)
    emit("expcone.dA_data", dA.tocsc().data)
    emit("expcone.db", db)
    emit("expcone.dc", dc)

    dx, dy, ds = D(dA, db, dc)
    emit("expcone.dx", dx)
    emit("expcone.dy", dy)
    emit("expcone.ds", ds)


if __name__ == "__main__":
    main()
