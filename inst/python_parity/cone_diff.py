"""Parity fixtures for the test_cone_prog_diff.py port.

For tests where Python uses cvxpy as the reference (test_proj_soc,
test_proj_psd, test_proj_exp, test_dprojection_exp), we generate the
same problem in Python and emit the cvxpy / diffcp reference values
as R-sourceable assignments.

Run with:

    uv run python inst/python_parity/cone_diff.py
"""

import numpy as np
import scipy.sparse as sparse
import cvxpy as cp

import diffcp.cone_program as cone_prog
import diffcp.cones as cone_lib
import diffcp.utils as utils


def _r_name(label):
    return label.replace(".", "_")


def emit(label, vec):
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


# ---- test_proj_soc reference: 15 random projections, n=100 -------
def emit_proj_soc():
    np.random.seed(0)
    n = 100
    print("# === proj_soc fixture ===")
    xs, projs = [], []
    for k in range(15):
        x = np.random.randn(n)
        p = cone_lib._proj(x, cone_lib.SOC, dual=False)
        xs.append(x); projs.append(p)
    emit("proj_soc_x", np.concatenate(xs))
    emit("proj_soc_p", np.concatenate(projs))


# ---- test_proj_psd reference: 15 random projections, n=10 --------
def emit_proj_psd():
    np.random.seed(0)
    n = 10
    print("# === proj_psd fixture ===")
    xs, projs = [], []
    for k in range(15):
        x = np.random.randn(n, n)
        x = x + x.T
        x_vec = cone_lib.vec_symm(x)
        p_vec = cone_lib._proj(x_vec, cone_lib.PSD, dual=False)
        xs.append(x_vec); projs.append(p_vec)
    emit("proj_psd_x", np.concatenate(xs))
    emit("proj_psd_p", np.concatenate(projs))


# ---- test_proj_exp reference: 15 random projections, n=9 ---------
def emit_proj_exp():
    np.random.seed(0)
    print("# === proj_exp fixture ===")
    xs, projs, dprojs = [], [], []
    for k in range(15):
        x = np.random.randn(9)
        p = cone_lib._proj(x, cone_lib.EXP, dual=False)
        p_dual = cone_lib._proj(x, cone_lib.EXP_DUAL, dual=False)
        xs.append(x); projs.append(p); dprojs.append(p_dual)
    emit("proj_exp_x", np.concatenate(xs))
    emit("proj_exp_p", np.concatenate(projs))
    emit("proj_exp_p_dual", np.concatenate(dprojs))


# ---- test_dprojection_exp: cvxpy parameter problem ----------------
def emit_dprojection_exp():
    print("# === dprojection_exp fixture ===")
    x_ = cp.Variable()
    _lam = cp.Parameter(1, nonneg=True)
    _lam.value = np.ones(1)
    objective = cp.Maximize(x_ + _lam * (cp.log(1 + x_) + cp.log(1 - x_)))
    problem = cp.Problem(objective)
    A, b, c, cone_dims = utils.scs_data_from_cvxpy_problem(problem)
    emit_csc("dprojexp_A", A)
    emit("dprojexp_b", b)
    emit("dprojexp_c", c)
    emit_cone_dict("dprojexp_cone_dict", cone_dims)
    ## analytical[0] for lam=1:  -1 + 1/sqrt(2) = -0.292893...
    print(f"dprojexp_analytical <- {-1 + 1.0/np.sqrt(2):.16g}")


def main():
    emit_proj_soc()
    emit_proj_psd()
    emit_proj_exp()
    emit_dprojection_exp()


if __name__ == "__main__":
    main()
