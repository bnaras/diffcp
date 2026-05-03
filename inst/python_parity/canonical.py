"""Canonical Clarabel + SCS test fixtures for the diffcp R package.

Mirrors the test cases in:
  diffcp/tests/test_clarabel.py::test_solve_and_derivative
  diffcp/tests/test_clarabel.py::test_psdcone
  diffcp/tests/test_scs.py::test_solve_and_derivative
  diffcp/tests/test_scs.py::test_warm_start
  diffcp/tests/test_scs.py::test_expcone

For each Python test we replicate the *exact* RNG sequence (seed +
order of consumption) so the perturbations Python uses are bit-exact
in the dumped data, allowing the R-side test to reconstruct them
without needing Python at run time.

Run with:

    uv run python inst/python_parity/canonical.py
"""

import numpy as np
import scipy.sparse as sparse
import cvxpy as cp

import diffcp.cone_program as cone_prog
import diffcp.utils as utils


def _r_name(label):
    """Convert a Python-style label like 'cl_lseq.A.x' into an R-safe
    variable name 'cl_lseq_A_x'."""
    return label.replace(".", "_")


def emit(label, vec):
    """Emit `label <- c(...)` as a valid R assignment."""
    s = ", ".join(f"{v:.16g}" for v in np.asarray(vec).flatten())
    print(f"{_r_name(label)} <- c({s})")


def emit_csc(label, M):
    M_csc = sparse.csc_matrix(M)
    nm = _r_name(label)
    print(f"{nm}_shape <- c({M_csc.shape[0]}L, {M_csc.shape[1]}L)")
    print(f"{nm}_A_i <- c({', '.join(str(r + 1) + 'L' for r in M_csc.indices)})")
    print(f"{nm}_A_p <- c({', '.join(str(p + 1) + 'L' for p in M_csc.indptr)})")
    emit(f"{label}_A_x", M_csc.data)


def emit_cone_dict(label, cd):
    parts = []
    for k, v in cd.items():
        if isinstance(v, (list, tuple)):
            if len(v) == 0:
                parts.append(f"{k} = integer(0)")
            else:
                parts.append(f"{k} = c({', '.join(str(int(x)) + 'L' for x in v)})")
        else:
            parts.append(f"{k} = {int(v)}L")
    print(f"{_r_name(label)} <- list({', '.join(parts)})")


# --- 1. test_clarabel.py::test_solve_and_derivative ---------------
def emit_lseq_clarabel():
    np.random.seed(0)
    A, b, c, cone_dims = utils.least_squares_eq_scs_data(20, 10)
    print("# === clarabel_lseq fixture ===")
    emit_csc("clarabel_lseq.A", A)
    emit("clarabel_lseq.b", b)
    emit("clarabel_lseq.c", c)
    emit_cone_dict("clarabel_lseq.cone_dict", cone_dims)

    for mode in ("lsqr", "dense"):
        x, y, s, D, DT = cone_prog.solve_and_derivative(
            sparse.csc_matrix(A), b, c, cone_dims, mode=mode,
            solve_method="Clarabel")
        dA = utils.get_random_like(
            A, lambda n: np.random.normal(0, 1e-6, size=n))
        db = np.random.normal(0, 1e-6, size=b.size)
        dc = np.random.normal(0, 1e-6, size=c.size)
        dx, dy, ds = D(dA, db, dc)
        emit(f"clarabel_lseq.{mode}.dA_data", dA.tocsc().data)
        emit(f"clarabel_lseq.{mode}.db", db)
        emit(f"clarabel_lseq.{mode}.dc", dc)
        emit(f"clarabel_lseq.{mode}.dx", dx)
        emit(f"clarabel_lseq.{mode}.dy", dy)
        emit(f"clarabel_lseq.{mode}.ds", ds)
        # Python re-seeds inside utils.least_squares_eq_scs_data, so
        # the second iteration would use stale RNG. We skip and just
        # dump per-mode data; the R test runs each mode independently
        # against its own dA/db/dc.
        np.random.seed(0)
        utils.least_squares_eq_scs_data(20, 10)  # re-consume RNG


# --- 2. test_clarabel.py::test_psdcone ----------------------------
def emit_psdcone_clarabel():
    DIM = 5
    X = cp.Variable(shape=(DIM, DIM), PSD=True)
    C = np.zeros((DIM, DIM))
    C[0, 0] = 1
    C[4, 4] = -1
    objective = cp.Minimize(cp.trace(C @ X))
    constraint = cp.trace(X) == 1
    problem = cp.Problem(objective, [constraint])
    A, b, c, cone_dims = utils.scs_data_from_cvxpy_problem(problem)
    print("# === clarabel_psdcone fixture ===")
    emit_csc("clarabel_psdcone.A", A)
    emit("clarabel_psdcone.b", b)
    emit("clarabel_psdcone.c", c)
    emit_cone_dict("clarabel_psdcone.cone_dict", cone_dims)
    sol_vec, _, _, _, _ = cone_prog.solve_and_derivative(
        sparse.csc_matrix(A), b, c, cone_dims, solve_method="Clarabel")
    emit("clarabel_psdcone.sol_vec", sol_vec)


# --- 3. test_scs.py::test_solve_and_derivative --------------------
def emit_lseq_scs():
    np.random.seed(0)
    A, b, c, cone_dims = utils.least_squares_eq_scs_data(20, 10)
    print("# === scs_lseq fixture ===")
    # A, b, c, cone_dict are the same as clarabel_lseq (same seed) but
    # we re-emit them under the scs name for clarity.
    emit_csc("scs_lseq.A", A)
    emit("scs_lseq.b", b)
    emit("scs_lseq.c", c)
    emit_cone_dict("scs_lseq.cone_dict", cone_dims)

    for mode in ("lsqr", "dense"):
        x, y, s, D, DT = cone_prog.solve_and_derivative(
            sparse.csc_matrix(A), b, c, cone_dims,
            eps=1e-10, mode=mode, solve_method="SCS")
        dA = utils.get_random_like(
            A, lambda n: np.random.normal(0, 1e-6, size=n))
        db = np.random.normal(0, 1e-6, size=b.size)
        dc = np.random.normal(0, 1e-6, size=c.size)
        dx, dy, ds = D(dA, db, dc)
        emit(f"scs_lseq.{mode}.dA_data", dA.tocsc().data)
        emit(f"scs_lseq.{mode}.db", db)
        emit(f"scs_lseq.{mode}.dc", dc)
        emit(f"scs_lseq.{mode}.dx", dx)
        emit(f"scs_lseq.{mode}.dy", dy)
        emit(f"scs_lseq.{mode}.ds", ds)
        np.random.seed(0)
        utils.least_squares_eq_scs_data(20, 10)


# --- 4. test_scs.py::test_warm_start ------------------------------
def emit_warm_start_scs():
    np.random.seed(0)
    A, b, c, cone_dims = utils.least_squares_eq_scs_data(20, 10)
    print("# === scs_warm_start fixture ===")
    emit_csc("scs_warm_start.A", A)
    emit("scs_warm_start.b", b)
    emit("scs_warm_start.c", c)
    emit_cone_dict("scs_warm_start.cone_dict", cone_dims)
    x, y, s, _, _ = cone_prog.solve_and_derivative(
        sparse.csc_matrix(A), b, c, cone_dims, eps=1e-9, solve_method="SCS")
    emit("scs_warm_start.x", x)
    emit("scs_warm_start.y", y)
    emit("scs_warm_start.s", s)


# --- 5. test_scs.py::test_expcone ---------------------------------
def emit_expcone_scs():
    np.random.seed(0)
    n = 10
    y = cp.Variable(n)
    obj = cp.Minimize(- cp.sum(cp.entr(y)))
    const = [cp.sum(y) == 1]
    prob = cp.Problem(obj, const)
    A, b, c, cone_dims = utils.scs_data_from_cvxpy_problem(prob)
    print("# === scs_expcone fixture ===")
    emit_csc("scs_expcone.A", A)
    emit("scs_expcone.b", b)
    emit("scs_expcone.c", c)
    emit_cone_dict("scs_expcone.cone_dict", cone_dims)

    for mode in ("lsqr", "dense"):
        x, y, s, D, DT = cone_prog.solve_and_derivative(
            sparse.csc_matrix(A), b, c, cone_dims,
            solve_method="SCS", mode=mode, eps=1e-10)
        dA = utils.get_random_like(
            A, lambda n: np.random.normal(0, 1e-6, size=n))
        db = np.random.normal(0, 1e-6, size=b.size)
        dc = np.random.normal(0, 1e-6, size=c.size)
        dx, dy, ds = D(dA, db, dc)
        emit(f"scs_expcone.{mode}.dA_data", dA.tocsc().data)
        emit(f"scs_expcone.{mode}.db", db)
        emit(f"scs_expcone.{mode}.dc", dc)
        emit(f"scs_expcone.{mode}.dx", dx)
        emit(f"scs_expcone.{mode}.dy", dy)
        emit(f"scs_expcone.{mode}.ds", ds)
        # Re-seed for next mode, mirroring how Python's outer
        # np.random.seed(0) plus the cvxpy canonicalization sets
        # state independent of mode iteration.
        np.random.seed(0)
        utils.scs_data_from_cvxpy_problem(prob)


def main():
    emit_lseq_clarabel()
    emit_psdcone_clarabel()
    emit_lseq_scs()
    emit_warm_start_scs()
    emit_expcone_scs()


if __name__ == "__main__":
    main()
