import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
argvs = sys.argv #Store arguments, in this case the path of the data files of thefoci ID to be analyzed
if len(argvs) <2:
	print('This script needs a cell data file and foci data file')
	sys.exit
plt.style.use('ggplot')
cellfilepath = argvs[1]
focifilepath = argvs[2]
celldf = pd.read_csv(cellfilepath,sep='\t')
print(celldf)