


import os


# fref: path of reference table (plaintext)
# ftest: path of test table
# columns: list of column indices to compare, ex. [0,3]
def CompareTable(fref, ftest, columns, delta=1E-6, skiprows=3):
  maxcol = max(columns)
  fr = open(fref, 'r')
  ft = open(ftest, 'r')
  for _ in range(0,skiprows):
    fr.readline()
    ft.readline()
  refline = fr.readline()
  testline = ft.readline()
  while ( len( refline.split()) >= maxcol and
          len(testline.split()) >= maxcol ):
    for c in columns:
      ref =  float( refline.split()[c])
      test = float(testline.split()[c])
      if (abs(ref-test) > delta):
        print(f'Mismatch between {fref} and {ftest}:')
        print(refline)
        print(testline)
        print(f"Column {c} mismatched ( {str(ref)} vs {str(test)} )")
        return False
    refline = fr.readline()
    testline = ft.readline()

  return True


