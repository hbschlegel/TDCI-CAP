


import os
import sys

# Get testing root
up = os.path.dirname # shorthand
testroot = up(up(up(os.path.abspath(__file__))))
sys.path.insert(0, testroot) # allow testutils import from testroot
import testutils
  

def run_test(refdir, testdir):
  # 2 Fields
  all_pass = True
  for i in range(1,3):
    # 2 E-strengths
    for e in range(1,3):
      ref = f"{refdir}/RESULTS-e{e}-d{i}.dat"
      test = f"{testdir}/RESULTS-e{e}-d{i}.dat"
      # 1=time (fs), 4=field1, 9=norm2, 10=rate(fs-1) 
      cols = [1, 4, 9, 10]
      pass_test = testutils.CompareTable(ref, test, cols)
      if not pass_test: all_pass = False
  return all_pass
  

