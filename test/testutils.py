#!/usr/bin/env python3
from __future__ import annotations
import math
from dataclasses import dataclass, field
from typing import Dict, List, Tuple
import numpy as np
import h5py

FS_PER_AU = 2.4188843265857e-2  # fs per atomic unit of time

# PASS/WARN/FAIL collector
@dataclass
class Check:
  name: str
  messages: List[Tuple[str, str]] = field(default_factory=list)
  failures: int = 0
  warnings: int = 0
  def pass_(self, msg: str): self.messages.append(("PASS", msg))
  def warn(self, msg: str):  self.messages.append(("WARN", msg)); self.warnings += 1
  def fail(self, msg: str):  self.messages.append(("FAIL", msg)); self.failures += 1
  def status(self) -> str:
    if self.failures: return "FAIL"
    if self.warnings: return "WARN"
    return "PASS"
  def summarize(self) -> int:
    print(f"\n=== {self.name} [{self.status()}] ===")
    for lvl, msg in self.messages:
      print(f"- {lvl}: {msg}")
    return 2 if self.failures else (1 if self.warnings else 0)


def read_group(h5_path: str, group: str) -> Dict[str, np.ndarray]:
  with h5py.File(h5_path, "r") as h5:
    g = h5[group]  # raises if missing
    return {k: v[...] for k, v in g.items() if isinstance(v, h5py.Dataset)}



def dipole_rms_xyz(mux: np.ndarray, muy: np.ndarray, muz: np.ndarray):
  x = mux - float(np.mean(mux))
  y = muy - float(np.mean(muy))
  z = muz - float(np.mean(muz))
  rms = lambda v: float(np.sqrt(np.mean(v*v)))
  return rms(x), rms(y), rms(z)

# During an elliptical pulse, the dipole vector may lag the field
# therefore the ellips(-ish shape) traced by the dipole vector may not be aligned with the field axes
# Borrowing the idea of covariance ellipses from machine vision, we can determine the ellipticity
#   of the dipole vector path by comparing the eigenvalues of the covariance matrix of (mux, muz)
# https://carstenschelp.github.io/2018/09/14/Plot_Confidence_Ellipse_001.html
# https://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix
#
# returns the "circularity" (λ_min/λ_max), 1 is circular, 0 is linear (very elliptical).
def dipole_circularity(mux: np.ndarray, muz: np.ndarray) -> float:
  # center data
  x = mux - float(np.mean(mux))
  z = muz - float(np.mean(muz))

  C = np.cov(np.vstack([x, z]))
  evals = np.linalg.eigvalsh(C); evals.sort()

  circ = float(0.0 if evals[-1] <= 0 else max(0.0, min(1.0, evals[0]/evals[1])))
  return circ

# Use shoelace formula to calculate the signed area sweeped by dipole vector
# positive is counter-clockwise, negative is clockwise (along y axis)
# when field normal is along y+, cirr should give positive, cirl should give negative.
# https://en.wikipedia.org/wiki/Shoelace_formula
def shoelace_signed_area(mux: np.ndarray, muz: np.ndarray) -> float:
  # center data
  x = mux - float(np.mean(mux))
  z = muz - float(np.mean(muz))
  x2 = np.r_[x, x[0]] # append first point, wraparound term in shoelace formula
  z2 = np.r_[z, z[0]]
  area = 0.5*np.sum(x2[:-1]*z2[1:] - x2[1:]*z2[:-1])
  return float(area)


# Compare test h5 to reference h5
# Returns (ok: bool, details: List[str])
#   where "details" is a list of descriptions of what is mismatched
def compare_to_reference(cur_h5: str, ref_h5: str, group: str, keys: List[str],
                         rtol: float = 1e-6, atol: float = 1e-12):

  details: List[str] = []
  try:
    cur = read_group(cur_h5, group)
    ref = read_group(ref_h5, group)
  except Exception as e:
    return False, [f"Error opening {group}: {e}"]
  ok = True
  for k in keys:
    if k not in cur or k not in ref:
      ok = False; details.append(f"{group}/{k}: dataset missing.")
      continue
    a, b = np.asarray(cur[k]), np.asarray(ref[k])
    if a.shape != b.shape:
      ok = False; details.append(f"{group}/{k}: shape {a.shape} vs {b.shape}."); continue
    if np.issubdtype(a.dtype, np.floating) or np.issubdtype(b.dtype, np.floating):
      if not np.allclose(a, b, rtol=rtol, atol=atol):
        ok = False
        diff = float(np.max(np.abs(a - b) / (np.abs(b) + 1e-16)))
        details.append(f"{group}/{k}: max rel err {diff:.2e} > tol {rtol:g}.")
    else:
      if not np.array_equal(a, b):
        ok = False; details.append(f"{group}/{k}: values differ.")
  return ok, details

