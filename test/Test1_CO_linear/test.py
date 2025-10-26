"""
Test3_CO_linear_angles — trotter_linear (no SOC)

Li 2011 data on orientation dependent ionization yield
https://doi.org/10.1103/PhysRevA.84.043429
Spiewanowski and Madsen 2015 ASFA calculations
https://doi.org/10.1103/PhysRevA.91.043406

This test checks:
  (1) Angle trend: Y(0°) > Y(45°) > Y(90°) for CO with the bond along +z,
  (2) Dipole aligns with field for theta=0° (z-dominant),
  (3) Reference-vs-current HDF5 match for a representative group.
"""

from __future__ import annotations
import os
from testutils import (
    Check, read_group,
    dipole_rms_xyz, compare_to_reference
)

def run(root=".", tol_ref: float=1e-6, tol_phys: float=0.05) -> int:
  T = "Test1_CO_linear"
  check = Check(T)

  cur = os.path.join(root, T, "rundir/t/data.h5")
  ref = os.path.join(root, T, "t.ref.h5")

  for grpstr in ["direction_1/field_1", "direction_2/field_1", "direction_3/field_1"]:
    ok, details = compare_to_reference(
      cur, ref, grpstr,
      keys=["time","efieldx","efieldy","efieldz","norm2","rate","mux","muy","muz"],
      rtol=tol_ref
    )
  if ok: check.pass_("Reference match (direction_1/field_1).")
  else:  [check.fail(d) for d in details]

  # Y(0°) > Y(45°) > Y(90°)
  g0  = read_group(cur, "direction_1/field_1")  # (0°,0°)
  g45 = read_group(cur, "direction_2/field_1")  # (45°,0°)
  g90 = read_group(cur, "direction_3/field_1")  # (90°,0°)
  y0, y45, y90 = 1-float(g0["norm2"][0,-1]), 1-float(g45["norm2"][0,-1]), 1-float(g90["norm2"][0,-1])
  if (y0 > (1+tol_phys)*y45) and (y45 > (1+tol_phys)*y90):
    check.pass_(f"Angle trend OK: {y0:.3e} > {y45:.3e} > {y90:.3e}.")
  else:
    check.fail(f"Angle trend violated: 0°={y0:.3e}, 45°={y45:.3e}, 90°={y90:.3e}.")

  # Dipole response along field
  rx, ry, rz = dipole_rms_xyz(g0["mux"][0], g0["muy"][0], g0["muz"][0])
  ratio = rz / (max(rx, ry) + 1e-12)
  if ratio >= 2.0: check.pass_(f"Dipole aligned (0°): RMSz/RMSperp={ratio:.2f}.")
  else:            check.warn (f"Dipole weakly aligned (0°): RMSz/RMSperp={ratio:.2f}.")


  return check.summarize()

