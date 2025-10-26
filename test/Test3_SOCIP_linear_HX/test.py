"""
Test3_SOCIP_linear_HX , SOC sensitivity in Ztrotter_linear

HCl+ exhibits a change in yield when SOC effects are included.
Data published by Lee, Hoerner, Li, Schlegel in 2020
https://doi.org/10.1063/5.0034807

The lighter HF+ should have negligible effect from SOC inclusion.

This test checks:
  (1) HF: |ΔY|/Y_nonSOC < 1% (small SOC effect)
  (2) HCl: |ΔY|/Y_nonSOC > 5% (non-negligible SOC effect)
  (3) regression vs reference HDF5 for SOC on/off
"""


from __future__ import annotations
import os
from testutils import Check, read_group, ion_yield_from_norm2, compare_to_reference

def run(root=".", tol_ref: float=1e-6, tol_phys: float=0.05) -> int:
  T = "Test3_SOCIP_linear_HX"
  check = Check(T)

  def yield_for(runname: str) -> float:
    test = os.path.join(root, T, dirpath, f"rundir/{runname}/data.h5")
    ref  = os.path.join(root, T, f"{runname}.ref.h5")
    ok, details = compare_to_reference(
      test, ref, "direction_1/field_1",
      keys=["time","efieldz","norm2","rate","mux","muz"], rtol=tol_ref
    )
    if ok: check.pass_(f"{runname} vs ref OK.")
    else:  [check.fail(d) for d in details]
    testnorm2 = read_group(test, "direction_1/field_1")["norm2"][0,-1]
    testyield = 1-float(testnorm2)
    return testyield

  y_HF_1  = yield_for("t_HF_soc1")
  y_HF_0  = yield_for("t_HF_soc0")
  y_HCl_1 = yield_for("t_HCl_soc1")
  y_HCl_0 = yield_for("t_HCl_soc0")

  fracHF  = abs(y_HF_1  - y_HF_0)  / max(y_HF_0,  1e-30)
  fracHCl = abs(y_HCl_1 - y_HCl_0) / max(y_HCl_0, 1e-30)

  if fracHF <= 0.01:
    check.pass_(f"HF SOC yield change small as expected: {y_HF_1:.3f} SOC vs {y_HF_0:.3f} No SOC; (<1% yield change)")
  else:
    check.warn(f"HF SOC yield change larger than expected: {y_HF_1:.3f} SOC vs {y_HF_0:.3f} No SOC; {fracHF:.2%} >1% yield change")

  if fracHCl >= max(0.05, tol_phys):
    check.pass_(f"HCl SOC yield change >5% as expected: {y_HCl_1:.3f} SOC vs {y_HCl_0:.3f} No SOC ({fracHCl:.2%} change)")
  else:
    check.fail(f"HCl SOC yield change too small: {y_HCl_1:.3f} SOC vs {y_HCl_0:.3f} No SOC ({fracHCl:.2%} <5% change)")

  return check.summarize()

