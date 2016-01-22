import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import cm
import sys,re
pd.options.mode.chained_assignment = None  # default='warn'
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
#1pixel=0.0435µm
pix = 0.0435

resultdir = '/Volumes/Kikuchi-SSD/150701F3-muNS-Per2/F3-muNS-30deg-M9glycas-0.2mMIPTG1h-bin2-TxRed300-GFP500-1minint_1/Python_Processed/'

celldf = pd.read_csv(cellfilepath,sep='\t')
celldf = celldf.drop('Mean',1)
celldf = celldf.drop('X',1)
celldf = celldf.drop('Y',1)
celldf = celldf.drop('Perim.',1)
celldf = celldf.drop('Major',1)
celldf = celldf.drop('Minor',1)


#Original data is in reverse-chronological because of backwards tracking in ImageJ
#Reverse order here so dataframe matches chronological order
celldf = celldf.reindex(index=celldf.index[::-1])
celldf = celldf.reset_index(drop=True)
celldf.Area = celldf.Area*(pix**2)
celldf.XM = celldf.XM*pix
celldf.YM = celldf.YM*pix
#celldf.Area.plot(label = 'Cell Area').legend(loc='upper right')

#Read in foci data
focidf = pd.read_csv(focifilepath,sep='\t')
#Remove irrevelant data
focidf = focidf.drop('Label', 1)
focidf = focidf.drop('Min', 1)
focidf = focidf.drop('Max', 1)
focidf = focidf.drop('Slice', 1)
focidf = focidf.drop(399) #the cell ROI data lacks the first slice, so delete it here to match
focidf.XM = focidf.XM*pix
focidf.YM = focidf.YM*pix
#Reindex
focidf =focidf.reindex(index=focidf.index[::-1])
focidf = focidf.reset_index(drop=True)

#Read in foci brightness data
brightness = pd.read_csv(brightnessfilepath,sep='\t')
brightness = brightness.drop('Min',1)
brightness = brightness.drop(' ',1)
#brightness data is in chronological order
brightness = brightness.drop(399)

#CELL CENTER DISPLACEMENT
cellx = celldf['XM']
celly = celldf['YM']
focix = focidf['XM']
fociy = focidf['YM']
#Calcluate distance between cell center and muNS focus
cellfocidisp = ((celldf.XM-focidf.XM)**2+(celldf.YM-focidf.YM)**2)**(1/2)
#Calculate displacement of cell center
celldisp = pd.Series(((cellx[i+1]-cellx[i])**2+(celly[i+1]-celly[i])**2)**(1/2) for i in range (len(cellx)-1))
#Calculate per-frame displacement of foci in Cartesian, corrected by cell center displacement
cartdisp = pd.Series(((focix[i+1]-focix[i]-cellx[i+1]+cellx[i])**2+(fociy[i+1]-fociy[i]-celly[i+1]+celly[i])**2)**(1/2) for i in range (len(cellx)-1))

#Calculate per-frame displacement of muNS foci in polar coordinates
#convert cell ROI long axis angle to radians
rad = np.radians(celldf['Angle'])
#convert foci polar coordinate theta taking into account the ROI's angle
radfoci = rad-np.arctan2(focidf['YM'],focidf['XM'])

polardisp= pd.Series((cellfocidisp[i+1]**2+cellfocidisp[i]**2-2*cellfocidisp[i]*cellfocidisp[i+1]*np.cos(radfoci[i+1]-radfoci[i]))**(1/2) for i in range (len(cellfocidisp)-1))
#polardisp.plot(label = 'Polar Displacement').legend(loc='center left', bbox_to_anchor=(1, 0.5))

#extract division timepoints
division = celldisp[celldisp > 0.4].index

#Replace the outliers in Polar displacements (that comes from cell division) with mean
for x in division:
        polardisp[x] = polardisp.mean()
        cartdisp[x] = cartdisp.mean()

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
celldf['Fluorescence_Intensity'] = brightness.Max.head(10).mean()
celldf['Generation'] = generation
celldf['rcostheta'] = cellfocidisp*np.cos(radfoci)
celldf.to_csv(resultdir+'GeneralData/data_%s.csv'%fociID, sep=',')

## create independent dataframe for each generation
gendict = {x: celldf[celldf.Generation == x] for x in range(1,len(division)+2)}
for x in range(1,len(division)+2):
    gendict[x]['t'] = np.arange(len(gendict[x]))
radigidict = {}
genradigidict = {}
linradigydict = {}
for i in gendict:
	radigidict[i] = pd.Series((x-gendict[i].rcostheta.mean())**2 for x in gendict[i].rcostheta)
	genradigidict[i] = np.sqrt(radigidict[i].mean())
	linradigydict[i] = radigidict[i].sum()

genradigi = pd.DataFrame.from_dict(genradigidict, orient="index")
linradigydf = pd.DataFrame.from_dict(linradigydict, orient='index')
linradigydf = linradigydf.drop(1)
linradigydf = linradigydf.drop(linradigydf.tail(1).index)


#Set first 10 measurements from max intensity as characteristic foci brightness
#Representative of muNS particle size
foci_brightness = brightness.Max.head(10).mean()

#Set Per-Generation DataFrame
gendf = pd.DataFrame()

gendf['Generation_Time']=celldf.groupby('Generation').size()
gendf['Elongation_Rate'] =(np.log(celldf.groupby('Generation')["Area"].last() /celldf.groupby('Generation')["Area"].first()))/celldf.groupby('Generation').size()
gendf['Fluorescence_Intensity'] = foci_brightness
gendf['Max_Cartesian_Displacement']= celldf.groupby('Generation')['Cartesian_Displacement'].max()
gendf['Mean_Cartesian_Displacement'] = celldf.groupby('Generation')['Cartesian_Displacement'].mean()
gendf['CV_Cartesian_Displacement'] = celldf.groupby('Generation')['Cartesian_Displacement'].std()/celldf.groupby('Generation')['Cartesian_Displacement'].mean()
gendf['Total_Cartesian_Displacement'] = celldf.groupby('Generation')['Cartesian_Displacement'].sum()
gendf['Mean_Polar_Displacement'] = celldf.groupby('Generation')['Polar_Displacement'].mean()
gendf['Max_Polar_Displacement']= celldf.groupby('Generation')['Polar_Displacement'].max()
gendf['CV_Polar_Displacement'] = celldf.groupby('Generation')['Polar_Displacement'].std()/celldf.groupby('Generation')['Polar_Displacement'].mean()
largedf = celldf[celldf['Polar_Displacement'] > 0.4]
gendf['Polar_Displacements_Over_04'] = largedf.groupby('Generation')['Polar_Displacement'].count()
gendf['Foci_ID'] = fociID

gendf['Radius_Gyration'] = genradigi
gendf.index.name = 'Generation'

#omit first and last generation
gendf2 = gendf.drop(1)
gendf2 = gendf2.drop(gendf.tail(1).index)

gendf.to_csv(resultdir+'Generation/generationdata_%s.csv'%fociID, sep=',')
gendf2.to_csv(resultdir+'MiddleGeneration/middlegenerationdata_%s.csv'%fociID, sep=',')


lineage = pd.DataFrame()
#create dataframe omitting 1st and last generations
lindf = celldf[(celldf.Generation>1)&(celldf.Generation<(len(division)+1))]
largelindf = lindf[lindf['Polar_Displacement'] > 0.4]

lineage['ID'] = pd.Series(fociID)
lineage['Divisions'] = len(division)
lineage['Lineage_Length'] =gendf2.Generation_Time.sum()
lineage['Fluorescence_Intensity'] =foci_brightness
lineage['Lineage_Elongation_Rate'] = (np.log(lindf.groupby('Generation')["Area"].last() /lindf.groupby('Generation')["Area"].first())).sum()/lineage.Lineage_Length
lineage['Lineage_Fitness'] = (1/400)*len(division)*np.log(2)
lineage['Mean_Cartesian_Displacement'] = lindf['Cartesian_Displacement'].mean()
lineage['CV_Cartesian_Displacement'] = lindf['Cartesian_Displacement'].std(ddof=1)/lindf['Cartesian_Displacement'].mean()
lineage['Max_Cartesian_Displacement'] = lindf['Cartesian_Displacement'].max()
lineage['Mean_Polar_Displacement'] = lindf['Polar_Displacement'].mean()
lineage['CV_Polar_Displacement'] = lindf['Polar_Displacement'].std(ddof=1)/lindf['Polar_Displacement'].mean()
lineage['SEM_Polar_Displacement']= lindf['Polar_Displacement'].std(ddof=1)/np.sqrt(len(lindf['Polar_Displacement']))
lineage['Max_Polar_Displacement']= lindf['Polar_Displacement'].max()
lineage['Polar_Displacements_Over_04'] = largelindf['Polar_Displacement'].count()
lineage['Mean_Radius_Gyration'] = genradigi.mean()
lineage['Radius_Gyration'] = np.sqrt(linradigydf.sum()/gendf2.Generation_Time.sum())

lineage.to_csv(resultdir+'LineageData/lineagedata_%s.csv'%fociID, sep=',')




#Create time-coloured trajectory
from matplotlib.collections import LineCollection
from matplotlib import cm

xy = focidf[['XM','YM']].values
z = np.linspace(0, 1, len(focidf))

min, max = (0, len(focidf))
step = 10 #step in color bar

#Using contourf to provide colorbar info, then clearing the figure
Q = [[0,0],[0,0]]
levels = range(min,max+step,step)
CS3 = plt.contourf(Q, levels, cmap=cm.jet)
plt.clf()
lc = LineCollection(zip(xy[:-1], xy[1:]), array=z, cmap=cm.jet)
fig,  ax = plt.subplots(1, 1)
ax.add_collection(lc)
ax.margins(0.1)
ax2= fig.colorbar(CS3)
ax2.set_label('time [min]')
#plt.show() #required to generate trajectory plot in Jupyter

#Calculate MSD
t_step = 1

msddf = focidf.drop(' ',1)
msddf.index=pd.Index(np.arange(len(msddf)), name='t_stamp')
msddf['t'] = np.arange(len(msddf)) * t_step

def compute_msd(trajectory, t_step, coords=['XM', 'YM']):
    numberofDeltaT = int(np.floor(len(msddf)/4))
    delays = trajectory.t[:numberofDeltaT]
    shifts = np.floor(delays[:numberofDeltaT]/t_step).astype(np.int)
    msds = np.zeros(numberofDeltaT)
    for i, shift in enumerate(shifts):
        diffs = trajectory[coords] - trajectory[coords].shift(-shift)
        sqdist = np.square(diffs).sum(axis=1)
        msds[i] = sqdist.mean()
    return delays, msds

delays, msds = compute_msd(msddf, t_step=t_step)

#Generate Summary Fig

#plt.suptitle('Summary of Foci %s' %fociID,  fontsize=20, y = 1.02)

ax1 = plt.subplot2grid((3,2), (0,0), colspan=2)
ax2 = plt.subplot2grid((3,2), (1,0), colspan=2)
ax3 = plt.subplot2grid((3,2), (2, 0))
ax4 = plt.subplot2grid((3,2), (2, 1))

celldf.Area.plot(ax=ax1, figsize = (12,15))
ax1.set_title('Cell Area')
ax1.set_xlabel('Time [min]')
ax1.set_ylabel('Pixels')

polardisp.plot(label = 'Polar Displacement', alpha = 0.7, ax=ax2)
cartdisp.plot( label = 'Cartesian Displacement', alpha = 0.7, ax=ax2)
ax2.set_title('Per-Frame Displacement')
ax2.legend(loc='upper left')
ax2.set_xlabel('Time [min]')
ax2.set_ylabel('Pixels')


ax3.add_collection(lc,autolim=True)
ax3.autoscale_view()
ax3.set_title('Particle Trajectory')

ax4.plot(delays, msds, 'r')
ax4.set_title('Mean Squared Displacement')
ax4.set_xlabel('Steps')
ax4.set_ylabel('MSD')

plt.tight_layout()

filename = "Summary_%s.pdf"%fociID
savedir = resultdir+'Graphs/'
#plt.savefig(savedir+filename)
