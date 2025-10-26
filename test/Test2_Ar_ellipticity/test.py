"""
Test2_Ar_ellipticity — ellipticity suppression (Ar)

Ar ionization yield decreases with increasing ellipticity,
  with circular (ε=1) < elliptical (ε≈0.5) < linear (ε≈0).

Wang 2014
https://doi.org/10.1103/PhysRevA.90.013422


This test checks:
  reference HDF5
  Y(ε=1) < Y(ε=0.5) for Ar
  more elliptical field induces more elliptical dipole in x-z plane
  dipole path matches chirality of field
"""


from __future__ import annotations
import os
from testutils import Check, read_group, dipole_circularity, compare_to_reference, shoelace_signed_area


def run(root=".", tol_ref: float=1e-6, tol_phys: float=0.05) -> int:
    T = "Test2_Ar_ellipticity"
    check = Check(T)

    def get_muxz(runname: str):
      test = os.path.join(root, T, f"rundir/{runname}/",  "data.h5")
      ref = os.path.join(root, T, f"{runname}.ref.h5")
      ok, details = compare_to_reference(
          test, ref, "direction_1/field_1",
          keys=["time","efieldx","efieldz","norm2","rate","mux","muz"], rtol=tol_ref
      )
      if ok: check.pass_(f"{runname} vs ref OK.")
      else:  [check.fail(d) for d in details]
      mux = read_group(test, "direction_1/field_1")["mux"][0]
      muz = read_group(test, "direction_1/field_1")["muz"][0]
      norm2 = read_group(test, "direction_1/field_1")["norm2"][0,-1]
      return mux, muz, norm2


    mux_hi, muz_hi, norm2_hi = get_muxz("t_eps1")
    mux_lo, muz_lo, norm2_lo = get_muxz("t_eps05")

    # Y(ε=1) < Y(ε=0.5)
    yH, yL = 1.0-float(norm2_hi), 1.0-float(norm2_lo)
    if yH < (1 - tol_phys)*yL: check.pass_(f"Y(ε=1) < Y(ε=0.5): {yH:.3e} < {yL:.3e}.")
    else:                      check.fail(f"Expected Y(ε=1) < Y(ε=0.5) by {tol_phys*100:.0f}%: {yH:.3e} vs {yL:.3e}.")

    # Dipole circularity increases with ε
    cH = dipole_circularity(mux_hi, muz_hi)
    cL = dipole_circularity(mux_lo, muz_lo)
    if cH > cL + 0.10: check.pass_(f"Dipole circularity increased with field circularity: {cH:.2f} (ε=1) > {cL:.2f} (ε=0.5).")
    else:              check.warn (f"Dipole circularity weird: {cH:.2f} vs {cL:.2f}.")

    # Dipole path matches chirality of field
    # when field normal is along y+, cirr should give positive, cirl should give negative.
    A_ = shoelace_signed_area(mux_hi, muz_hi)
    if (A_ > 0.0): check.pass_(f"Dipole path matches chirality of field (shoelace signed area: {A_:.2f})")
    else:          check.warn (f"Dipole path does not match chirality of field (shoelace signed area: {A_:.2f})")


    return check.summarize()



