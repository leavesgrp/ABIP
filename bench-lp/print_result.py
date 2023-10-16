import sys
import numpy as np
import pandas as pd


# usage: python *.py mittelmann.220125.xlsx
def format_object(s):
  if type(s) == str:
    return s
  return s.__format__('.3e')


df = pd.read_excel(sys.argv[1]).set_index(['name',
                                           'method']).applymap(format_object)
print(df.to_latex(longtable=True, multirow=True, multicolumn=True))
