#!/usr/bin/env python

# Run a simulation for each propagation type, and make sure the data outputted is identical


import sys, os, shutil, subprocess, time
import importlib.util

class Proctor:
  def __init__(self):
    self.testroot = os.path.dirname(os.path.realpath(__file__)) # Path of this script
    # Get path of subdirectories that contain .ref
    self.refdirs = []
    for d in os.listdir(self.testroot+"/tests/"):
      dp = os.path.join(self.testroot+"/tests/", d)
      if os.path.isdir(dp) and d.endswith(".ref"):
        self.refdirs.append(dp)

    self.binpath = self.testroot+"/../bin/tdci"
    print(self.refdirs)

  def RunTests(self):
    print("Starting Tests!")
    all_passed = True
    for d in self.refdirs:
      print("Starting : "+d)
      testdir = d[:-4]
      self.create_test(testdir)
      self.ExecuteTDCI(testdir)
      passed = self.ExecuteTest(d)
      if passed:
        print(f"Test Passed : {testdir}")
      else:
        all_passed = False
        print(f"Test FAILED! : {testdir}")
    pass
    if all_passed:
      print("All tests passed! Hooray!!")
    else:
      print("At least one test failed T__T ")

  def ExecuteTDCI(self, testdir):
    os.chdir(testdir)
    p = subprocess.Popen( self.binpath , shell=True)
    p.wait()
    os.chdir(self.testroot)

  def ExecuteTest(self, refdir):
    testdir = refdir[:-4]
    testjobpy = os.path.join(refdir, "testjob.py")
    # Dynamically load module
    spec = importlib.util.spec_from_file_location("testjob", testjobpy)
    testjob = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(testjob)

    # Run the test
    passed = False
    try:
      passed = testjob.run_test(refdir, testdir)
    finally: # Executed regardless of exceptions
      # Unload our dynamically loaded module
      if "testjob" in sys.modules:
        del sys.modules["testjob"]
    return passed
    
    
  # Create the files to run a test.
  #   Argument testdir must NOT have a trailing /, and must be the test directory, not the reference (.ref)
  def create_test(self, testdir):
    d = testdir[:]
    # Sanitize input
    if d == "": raise ValueError("Empty Input.")
    if d[-1] == "/" : raise ValueError(f"Trailing slash: {d}")
    if ".ref" in d: raise ValueError(f"Gave create_test the reference directory: {d}")
    if not os.path.exists(d+".ref"): raise ValueError(f"Can't find reference directory: {d}.ref")
    
    # Delete old test
    if os.path.exists(d): shutil.rmtree(d)
    os.mkdir(d)
    # Copy the symlink, do NOT copy the full TDCI.dat file
    shutil.copy( d+".ref/TDCI.dat" , d+"/TDCI.dat" , follow_symlinks=False)
    shutil.copy( d+".ref/input" , d+"/input" )
    return 0



a = Proctor()
a.RunTests()

