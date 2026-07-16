#!/usr/bin/env python3
"""
read_substructure_stiffness.py
Parse an Abaqus *SUBSTRUCTURE MATRIX OUTPUT .mtx (lower-triangular),
rebuild the full symmetric stiffness, and report directional stiffness.
Standard library only.
"""

import os

# ------------------------------ INPUTS ------------------------------
# Resolved next to this script so the cwd you launch from doesn't matter.
MTX_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "perf_K.mtx")
N_NODES  = 2                # retained reference nodes
NDOF     = 6                # DOF per node (3 trans + 3 rot)

# For a stiffness readout, tell it what you pulled:
PULL_NODE_INDEX = 2         # 1 = first retained node, 2 = second, ...
PULL_DOF        = 1         # 1=X 2=Y 3=Z 4=RX 5=RY 6=RZ
APPLIED_U       = 0.002     # prescribed displacement you applied
MEASURED_RF     = 4.71e4    # reaction force you read (optional cross-check; None to skip)
# --------------------------------------------------------------------

N = N_NODES * NDOF
DOF_NAMES = ["X", "Y", "Z", "RX", "RY", "RZ"]


def parse_mtx(path):
    """Return full NxN symmetric matrix from a lower-triangular .mtx."""
    vals = []
    in_matrix = False
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("**"):
                continue
            if s.startswith("*"):
                # Data before *MATRIX (e.g. the *USER ELEMENT DOF list) is not
                # matrix data, and would otherwise be read as stiffness values.
                in_matrix = s.upper().startswith("*MATRIX")
                continue
            if not in_matrix:
                continue
            # numbers may be comma- and/or space-separated; Fortran 'E' exponents ok
            for tok in s.replace(",", " ").split():
                try:
                    vals.append(float(tok))
                except ValueError:
                    pass  # skip stray labels

    expected = N * (N + 1) // 2
    if len(vals) != expected:
        raise ValueError(
            "Expected %d lower-tri entries, found %d. Check N_NODES/NDOF."
            % (expected, len(vals)))

    K = [[0.0] * N for _ in range(N)]
    k = 0
    # lower-triangular fill: row i has (i+1) entries
    for i in range(N):
        for j in range(i + 1):
            if k >= len(vals):
                raise ValueError(
                    "Ran out of numbers: expected %d lower-tri entries, got %d. "
                    "Check N_NODES/NDOF." % (N * (N + 1) // 2, len(vals)))
            K[i][j] = vals[k]
            K[j][i] = vals[k]      # mirror to upper triangle
            k += 1
    return K


def fmt(x):
    """Compact scientific / dotted-zero formatting for display."""
    if abs(x) < 1e-3 * 1.0:      # relative noise threshold applied later
        return x
    return x


def print_matrix(K):
    # column headers
    hdr = []
    for a in range(N_NODES):
        for d in range(NDOF):
            hdr.append("n%d%s" % (a + 1, DOF_NAMES[d]))
    # noise threshold = 1e-6 of the largest magnitude
    mx = max(abs(K[i][j]) for i in range(N) for j in range(N))
    tol = 1e-6 * mx

    colw = 12
    print(" " * 6 + "".join("%*s" % (colw, h) for h in hdr))
    for i in range(N):
        row = "%4s  " % hdr[i]
        for j in range(N):
            v = K[i][j]
            if abs(v) < tol:
                cell = "."
            else:
                cell = "%.4g" % v
            row += "%*s" % (colw, cell)
        print(row)


def check_health(K):
    # symmetry
    mx = max(abs(K[i][j]) for i in range(N) for j in range(N))
    asym = max(abs(K[i][j] - K[j][i]) for i in range(N) for j in range(N)) / mx
    print("\nSymmetry error (max|K-K^T|/max|K|) = %.2e  (want ~0)" % asym)

    # node-block diagonal twins (only meaningful for a symmetric 2-node part)
    if N_NODES == 2:
        print("Diagonal twins (node1 vs node2):")
        for d in range(NDOF):
            a, b = K[d][d], K[NDOF + d][NDOF + d]
            rel = abs(a - b) / max(abs(a), 1e-30)
            flag = "ok" if rel < 1e-3 else "DIFFER"
            print("  %-3s : %-12.5g  %-12.5g   %s"
                  % (DOF_NAMES[d], a, b, flag))


def directional_stiffness(K):
    idx = (PULL_NODE_INDEX - 1) * NDOF + (PULL_DOF - 1)
    kdiag = K[idx][idx]
    name = "node%d %s" % (PULL_NODE_INDEX, DOF_NAMES[PULL_DOF - 1])
    print("\n--- Directional stiffness ---")
    print("Pulled DOF: %s  (global index %d)" % (name, idx + 1))
    print("K(%d,%d) from matrix       = %.6g lb/in" % (idx + 1, idx + 1, kdiag))
    if MEASURED_RF is not None and APPLIED_U:
        kmeas = MEASURED_RF / APPLIED_U
        rel = abs(kmeas - kdiag) / max(abs(kdiag), 1e-30)
        print("K = RF/U = %.4g / %.4g = %.6g lb/in" % (MEASURED_RF, APPLIED_U, kmeas))
        print("Agreement                 = %.3f%%  (want < 1%%)" % (100 * rel))
    return kdiag


if __name__ == "__main__":
    K = parse_mtx(MTX_FILE)
    print("Full %dx%d reduced stiffness matrix (. = |value| < 1e-6*max):\n" % (N, N))
    print_matrix(K)
    check_health(K)
    directional_stiffness(K)