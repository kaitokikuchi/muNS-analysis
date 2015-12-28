import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys,re
argvs = sys.argv #Store arguments, in this case the path of the data files of thefoci ID to be analyzed
if len(argvs) <3:
	print('This script needs a cell data file, foci data file, and a foci brightness data file')
	sys.exit
plt.style.use('ggplot')
cellfilepath = argvs[1]
focifilepath = argvs[2]
brightnessfilepath = argvs[3]
#Set regex to extract foci ID
regex = re.compile(r"^[^_]*_[^_]*_(.*)[.]xls$")
#Store foci ID
fociID = re.sub(regex, r"\1",cellfilepath)

resultdir = '/Volumes/Kikuchi-SSD/150701F3-muNS-Per2/F3-muNS-30deg-M9glycas-0.2mMIPTG1h-bin2-TxRed300-GFP500-1minint_1/Python_Processed/'

celldf = pd.read_csv(cellfilepath,sep='\t')
celldf  = celldf.drop('Mean',1)
celldf = celldf.drop('X',1)
celldf = celldf.drop('Y',1)
celldf = celldf.drop('Perim.',1)
celldf = celldf.drop('Major',1)
celldf = celldf.drop('Minor',1)


#Original data is in reverse-chronological because of backwards tracking in ImageJ
#Reverse order here so dataframe matches chronological order
celldf =celldf.reindex(index=celldf.index[::-1])
celldf = celldf.reset_index(drop=True)
celldf.Area.plot(label = 'Cell Area').legend(loc='upper right')

#Read in foci data
focidf = pd.read_csv(focifilepath,sep='\t')
#Remove irrevelant data
focidf = focidf.drop('Label', 1)
focidf = focidf.drop('Min', 1)
focidf = focidf.drop('Max', 1)
focidf = focidf.drop('Slice', 1)
focidf = focidf.drop(399) #the cell ROI data lacks the first slice, so delete it here to match
#Reindex
focidf =focidf.reindex(index=focidf.index[::-1])
focidf = focidf.reset_index(drop=True)

#Read in foci brightness data
brightness = pd.read_csv(brightnessfilepath,sep='\t')
brightness = brightness.drop('Min',1)
#brightness data is in chronological order
brightness = brightness.drop(399)

#CELL CENTER DISPLACEMENT
cellx = celldf['XM']
celly = celldf['YM']
#Calcluate distance between cell center and muNS focus
cellfocidisp = ((celldf.XM-focidf.XM)**2+(celldf.YM-focidf.YM)**2)**(1/2)
#Calculate per-frame displacement of cell center
celldisp = pd.Series(((cellx[i+1]-cellx[i])**2+(celly[i+1]-celly[i])**2)**(1/2) for i in range (len(cellx)-1))



#Calculate per-frame displacement of muNS foci in polar coordinates
#convert cell ROI long axis angle to radians
rad = np.radians(celldf['Angle'])
#convert foci polar coordinate theta taking into account the ROI's angle
radfoci = rad-np.arctan2(focidf['YM'],focidf['XM'])

polardisp= pd.Series((cellfocidisp[i+1]**2+cellfocidisp[i]**2-2*cellfocidisp[i]*cellfocidisp[i+1]*np.cos(radfoci[i+1]-radfoci[i]))**(1/2) for i in range (len(cellfocidisp)-1))
#polardisp.plot(label = 'Polar Displacement').legend(loc='center left', bbox_to_anchor=(1, 0.5))

#Replace the outliers in Polar displacements (that comes from cell division) with mean
for x in polardisp:
	if x>8:
		polardisp[polardisp == x]=polardisp.mean()


#calculate cartesian distance
focix = focidf['XM']
fociy = focidf['YM']
cartdisp = pd.Series(((focix[i+1]-focix[i])**2+(fociy[i+1]-fociy[i])**2)**(1/2) for i in range (len(focix)-1))

#extract division timepoints
division = celldisp[celldisp > 10].index

#Set generation number for each timepoint
generation = pd.Series(range(len(polardisp)+1))

for i in range(len(division)+1):
    if i == 0: #Generation1
        generation[:division[i]+1]= 1
    elif i == len(division): #Last Generation
        generation[division[i-1]:]=i+1
    else:
        generation[division[i-1]+1:division[i]+1]= i+1


celldf['Foci X'] = focix
celldf['Foci Y'] = fociy
celldf['Cartesian_Displacement'] = cartdisp
celldf['Polar_Displacement'] = polardisp
celldf['Fluorescence_Intensity'] = brightness['Max']
celldf['Generation'] = generation
celldf.to_csv(resultdir+'data_%s.csv'%fociID, sep=',')


#Set first 10 brightness measurements as chaaracteristic foci brightness
foci_brightness = brightness.Max.head(10).mean()

#Set Per-Generation DataFrame
gendf = pd.DataFrame()

gendf['Generation_Time']=celldf.groupby('Generation').size()
gendf['Elongation_Rate'] =(np.log(celldf.groupby('Generation')["Area"].last() /celldf.groupby('Generation')["Area"].first()))/celldf.groupby('Generation').size()
gendf['Fluorescence_Intensity'] = foci_brightness

gendf['Mean_Cartesian_Displacement'] = celldf.groupby('Generation')['Cartesian_Displacement'].mean()
gendf['CV_Cartesian_Displacement'] = celldf.groupby('Generation')['Cartesian_Displacement'].std()/celldf.groupby('Generation')['Cartesian_Displacement'].mean()

gendf['Mean_Polar_Displacement'] = celldf.groupby('Generation')['Polar_Displacement'].mean()
gendf['CV_Polar_Displacement'] = celldf.groupby('Generation')['Polar_Displacement'].std()/celldf.groupby('Generation')['Polar_Displacement'].mean()

gendf.index.name = 'Generation'

#omit first and last generation
gendf = gendf.drop(1)
gendf = gendf.drop(gendf.tail(1).index)

gendf.to_csv(resultdir+'generationdata_%s.csv'%fociID, sep=',')
