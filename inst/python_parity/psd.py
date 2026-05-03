"""PSD-cone parity fixtures for the diffcp R package.

Mirrors the test cases in
  /Users/naras/research/cvxr/new_design/diffing/diffcp/tests/test_clarabel_psd.py

For each Python test, we build the same cvxpy Problem, canonicalize it
to SCS form via `scs_data_from_cvxpy_problem`, then emit R code that
reconstructs (A, b, c, cone_dict) plus the cvxpy reference objective.
The R port then runs `diffcp::solve_only` via both Clarabel and SCS
and asserts the same agreement bounds the Python tests use (atol=1e-4
on objective and x; atol=1e-5 on constraint residual).

Run with:

    uv run python inst/python_parity/psd.py
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


def emit_problem(name, A, b, c, cone_dict, cvxpy_obj):
    """Emit R code that reconstructs (A, b, c, cone_dict) + obj_ref."""
    A = sparse.csc_matrix(A)
    print(f"# === {name} ===")
    print(f"## {name}.shape: c({A.shape[0]}, {A.shape[1]})")
    rows = ", ".join(str(r + 1) for r in A.indices)             # 1-based
    ptrs = ", ".join(str(p + 1) for p in A.indptr)              # 1-based; first = 1
    data = ", ".join(f"{v:.16g}" for v in A.data)
    print(f"## {name}.A_i: c({rows})")
    print(f"## {name}.A_p: c({ptrs})")
    print(f"## {name}.A_x: c({data})")
    print(f"## {name}.b: c({', '.join(f'{v:.16g}' for v in b)})")
    print(f"## {name}.c: c({', '.join(f'{v:.16g}' for v in c)})")
    cd_parts = []
    for k, v in cone_dict.items():
        if isinstance(v, (list, tuple)):
            cd_parts.append(f"{k} = c({', '.join(str(int(x)) + 'L' for x in v)})")
        else:
            cd_parts.append(f"{k} = {int(v)}L")
    print(f"## {name}.cone_dict: list({', '.join(cd_parts)})")
    print(f"## {name}.obj_ref: {cvxpy_obj:.16g}")


# ---- 1. test_multiple_psd_cones_objective_match -----------------
def fixture_multiple_psd_cones():
    A = np.array([[1, 2, 3],
                  [2, 4, 5],
                  [3, 5, 6]], dtype=float)
    B = np.array([[7, 8, 9],
                  [8, 10, 11],
                  [9, 11, 12]], dtype=float)
    X = cp.Variable((3, 3), symmetric=True)
    y = cp.Variable(2)
    constraints = [
        y[0] * A + y[1] * B >> 0,
        X >> 0,
        cp.trace(A @ X) == 1,
        y >= 0,
    ]
    obj = cp.Minimize(cp.trace(X) + np.ones(2) @ y)
    prob = cp.Problem(obj, constraints)
    cvxpy_obj = prob.solve(solver=cp.CLARABEL)
    A_, b_, c_, cones_ = scs_data_from_cvxpy_problem(prob)
    return A_, b_, c_, cones_, cvxpy_obj


# ---- 2. test_single_psd_cone ------------------------------------
def fixture_single_psd_cone():
    n = 3
    C = np.eye(n)
    X = cp.Variable((n, n), symmetric=True)
    constraints = [X >> 0, cp.trace(X) == 1]
    obj = cp.Minimize(cp.trace(C @ X))
    prob = cp.Problem(obj, constraints)
    cvxpy_obj = prob.solve(solver=cp.CLARABEL)
    A_, b_, c_, cones_ = scs_data_from_cvxpy_problem(prob)
    return A_, b_, c_, cones_, cvxpy_obj


# ---- 3. test_mixed_cones ----------------------------------------
def fixture_mixed_cones():
    n = 2
    X = cp.Variable((n, n), symmetric=True)
    t = cp.Variable()
    A = np.array([[1, 0.5], [0.5, 2]])
    constraints = [
        X >> 0,
        t >= 0,
        cp.trace(A @ X) == 1,
        t <= 5,
    ]
    obj = cp.Minimize(cp.trace(X) + t)
    prob = cp.Problem(obj, constraints)
    cvxpy_obj = prob.solve(solver=cp.CLARABEL)
    A_, b_, c_, cones_ = scs_data_from_cvxpy_problem(prob)
    return A_, b_, c_, cones_, cvxpy_obj


# ---- 4. test_constraint_satisfaction ----------------------------
def fixture_constraint_satisfaction():
    A_mat = np.array([[1, 2, 3],
                      [2, 4, 5],
                      [3, 5, 6]], dtype=float)
    B_mat = np.array([[7, 8, 9],
                      [8, 10, 11],
                      [9, 11, 12]], dtype=float)
    X = cp.Variable((3, 3), symmetric=True)
    y_var = cp.Variable(2)
    constraints = [
        y_var[0] * A_mat + y_var[1] * B_mat >> 0,
        X >> 0,
        cp.trace(A_mat @ X) == 1,
        y_var >= 0,
    ]
    obj = cp.Minimize(cp.trace(X) + np.ones(2) @ y_var)
    prob = cp.Problem(obj, constraints)
    cvxpy_obj = prob.solve(solver=cp.CLARABEL)
    A_, b_, c_, cones_ = scs_data_from_cvxpy_problem(prob)
    return A_, b_, c_, cones_, cvxpy_obj


# ---- 5. test_derivative_lsqr_mode -------------------------------
def fixture_derivative_lsqr():
    n = 2
    C = np.array([[1.0, 0.3], [0.3, 2.0]])
    X = cp.Variable((n, n), symmetric=True)
    constraints = [X >> 0, cp.trace(X) == 1]
    obj = cp.Minimize(cp.trace(C @ X))
    prob = cp.Problem(obj, constraints)
    cvxpy_obj = prob.solve(solver=cp.CLARABEL)
    A_, b_, c_, cones_ = scs_data_from_cvxpy_problem(prob)
    return A_, b_, c_, cones_, cvxpy_obj


# ---- 6. test_psd_permutation_logic ------------------------------
def fixture_psd_permutation_logic():
    np.random.seed(0)  # match Python test's lack-of-seed deterministically
    n = 3
    C = np.random.randn(n, n)
    C = C @ C.T
    X = cp.Variable((n, n), symmetric=True)
    constraints = [X >> 0, cp.trace(X) == 1]
    obj = cp.Minimize(cp.trace(C @ X))
    prob = cp.Problem(obj, constraints)
    cvxpy_obj = prob.solve(solver=cp.SCS)
    A_, b_, c_, cones_ = scs_data_from_cvxpy_problem(prob)
    return A_, b_, c_, cones_, cvxpy_obj


def main():
    for name, fixture in [
        ("multi_psd",    fixture_multiple_psd_cones),
        ("single_psd",   fixture_single_psd_cone),
        ("mixed",        fixture_mixed_cones),
        ("constr_sat",   fixture_constraint_satisfaction),
        ("deriv_lsqr",   fixture_derivative_lsqr),
        ("perm_logic",   fixture_psd_permutation_logic),
    ]:
        A, b, c, cones, obj_ref = fixture()
        emit_problem(name, A, b, c, cones, obj_ref)


if __name__ == "__main__":
    main()
