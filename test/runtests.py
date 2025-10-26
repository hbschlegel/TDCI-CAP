#!/usr/bin/env python3
from __future__ import annotations
import argparse, importlib.util, os, sys


# Returns CompletedProcess which has members stdout stderr and returncode
def run(cmd, *, env=None, text=True):
  #print(f"Running : [ {cmd} ]\n    from directory: [ {os.getcwd()} ]")
  out = subprocess.run(cmd, shell=True, capture_output=True, text=text,
                        errors="replace", env=env if env is not None else os.environ )
  #if out.stdout != "":
  #  print(f"stdout : \n{out.stdout}")
  #if out.stderr != "":
  #  print(f"stderr : \n{out.stderr}")
  return out


def load_and_run(test_dir: str) -> int:
  # Load test.py from the test directory and execute its run() function
  test_py = os.path.join(test_dir, "test.py")
  if not os.path.isfile(test_py):
    print(f"Skipping {test_dir} (no test.py).")
    return 0
  # Allow "import testutils" from inside test.py
  root = os.path.abspath(os.path.dirname(__file__))
  if root not in sys.path:
    sys.path.insert(0, root)
  spec = importlib.util.spec_from_file_location("tdci_test_module", test_py)
  mod = importlib.util.module_from_spec(spec)
  assert spec.loader is not None
  spec.loader.exec_module(mod)
  run = getattr(mod, "run", None)
  if run is None:
    print(f"{test_py} has no run() function."); return 2
  return int(run(root=os.path.dirname(test_dir)))

def main():
  ap = argparse.ArgumentParser()
  ap.add_argument("--gaussian", default="", help="Path to folder with gdv executible (REQUIRED)")
  ap.add_argument("--root", default=".", help="Folder containing Test*/ subdirs")
  ap.add_argument("--only", default="", help="Comma-separated subset of test directory names")
  ap.add_argument("--no-execute", action='store_true', default=False, help="Do not execute gaussian and tdci, just test existing jobfiles.")
  args = ap.parse_args()

  tdcitestdir = os.path.dirname(os.path.realpath(__file__)) # directory this script resides in
  tdcibin = f"{tdcitestdir}/../bin/tdci"
  gaussdir = args.gaussian # directory containing gdv

  root = os.path.abspath(args.root)
  subdirs = sorted([d for d in os.listdir(root) if d.startswith("Test") and os.path.isdir(os.path.join(root, d))])
  if args.only:
    want = set([s.strip() for s in args.only.split(",") if s.strip()])
    subdirs = [d for d in subdirs if d in want]
  if not subdirs:
    print("No Test*/ subdirectories found."); sys.exit(1)

  overall = 0
  for d in subdirs:
    if not args.no_execute:
      os.chdir(d) # step into test dir
      run(f"bash run.sh {args.gaussian} {tdcibin}")
      os.chdir("..") # return to root dir
    
    code = load_and_run(os.path.join(root, d))
    overall = max(overall, code)

  print("\n==== Overall:", "PASS" if overall==0 else ("WARN" if overall==1 else "FAIL"))
  sys.exit(overall)

if __name__ == "__main__":
  main()

