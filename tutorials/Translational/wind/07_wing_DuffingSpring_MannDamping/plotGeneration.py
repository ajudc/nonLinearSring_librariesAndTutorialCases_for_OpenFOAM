#!/usr/bin/env python
# Make a legend for specific lines.
from __future__ import division


import numpy as np

import os, sys
from os import walk

import copy

import matplotlib
matplotlib.use('Agg') #para uso en cluster
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import rc
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker as mtick
from matplotlib.legend_handler import HandlerLine2D
import matplotlib.lines as mlines
import matplotlib as mpl
from matplotlib.patches import Polygon
from matplotlib.patches import Rectangle
from matplotlib.patches import ConnectionPatch
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import ticker

plt.style.use("classic")

##################################################################################################################

##################################################################################################################

##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##############                                                                                      ##############
##############                                                                                      ##############
##############                 FUNCTIONS                                                            ##############
##############                                                                                      ##############
##############                                                                                      ##############
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################

##################################################################################################################

##################################################################################################################

def checkFolder(orFolder,freeDecayCase):

	folders = []
	for (dirpath, dirnames, filenames) in walk(orFolder):
		folders.extend(dirnames)
		break

	if "postProcessing" in folders:

		orFolder=orFolder+"/postProcessing"

		folders = []
		for (dirpath, dirnames, filenames) in walk(orFolder):
			folders.extend(dirnames)
			break

		if not freeDecayCase:
			if "forceCoeffsWing" in folders:
				folderFC=orFolder+"/forceCoeffsWing"

			else:
				sys.exit("\n\nThere is no forceCoeffsWing folder. Scrpipt exiting\n\n")
		else:
			folderFC=""

		if "sixDoFRigidBodyState" in folders:
			folderDisp=orFolder+"/sixDoFRigidBodyState"
		else:
			sys.exit("\n\nThere is no sixDoFRigidBodyState folder. Scrpipt exiting\n\n")

		return folderFC,folderDisp

	else:

		sys.exit("\n\nThere is no postProcessing folder. Scrpipt exiting\n\n")

##################################################################################################################

##################################################################################################################

def read_forceCoeffs(filename,time_ini):
	#reads the openfoam generated file for the force coefficients
	#INPUTS
	#filename the full path to the file containing the time-history data of the force coefficients
	#time_ini the time from which the data will be read
	#OUTPUTS
	#time the time instants of the time-history
	#Cd,Cl,Cm the drag, lift and moment coefficients


	time=[];Cd=[];Cl=[];Cm=[]

	comp=True
	cont=0

	f = open(filename, 'r')

	while comp:

		line = f.readline()
		line = line.strip()
		columns=" ".join(line.split())
		columns=columns.split(" ")
		#print("columns= ",columns)
	
		if len(columns)>1:
			if columns[1]=='Time' or columns[0]=='Time':
				cont=cont+1
				comp=False
			else:
				cont=cont+1
		else:
			cont=cont+1

	f.close()

	columns=columns[1:]
	#print("columns= ",columns)

	numCols=len(columns)
	posFC=[]
	for i in range(3):#I am looking for Cd, Cl and Cm
		for j in range(numCols):
			if i==0:
				if columns[j]=='Cd':
					posFC.append(j)
			elif i==1:
				if columns[j]=='Cl':
					posFC.append(j)
			elif i==2:
				if columns[j]=='Cm' or columns[j]=='CmPitch':
					posFC.append(j)

	#print("posFC= ",posFC)

	#read the data
	###########################################

	#read headers
	f = open(filename, 'r')
	for i in range(cont):
		header1 = f.readline()

	#read data
	for line in f:
		line = line.strip()
		columns=" ".join(line.split())
		columns=columns.split(" ")
		#columns = line.split("\t")
		time.append(float(columns[0]))
		Cd.append(float(columns[posFC[0]]))
		Cl.append(float(columns[posFC[1]]))
		Cm.append(float(columns[posFC[2]]))

	f.close()

	time=np.array(time)
	Cd=np.array(Cd)
	Cl=np.array(Cl)
	Cm=np.array(Cm)

	#find postions with time greater than time_ini
	pos=np.argwhere(time>=time_ini)
	if len(pos)>0:
		if pos.ndim>1:
			pos=pos[:,0]

		time=time[pos]
		Cd=Cd[pos]
		Cl=Cl[pos]
		Cm=Cm[pos]
	else:
		time=[]
		Cd=[]
		Cl=[]
		Cm=[]

	return (time,Cd,Cl,Cm)

##################################################################################################################

##################################################################################################################

def read_displacements(filename,time_ini):
	#reads the openfoam generated file for the displacements
	#INPUTS
	#filename the full path to the file containing the time-history data of the displacements
	#time_ini the time from which the data will be read
	#OUTPUTS
	#time the time instants of the time-history
	#y,alpha the heave and pitch displacements


	time=[];y=[];alpha=[]

	comp=True
	cont=0

	f = open(filename, 'r')

	while comp:

		line = f.readline()
		line = line.strip()
		columns=" ".join(line.split())
		columns=columns.split(" ")
		#print("columns= ",columns)
	
		if len(columns)>1:
			if columns[1]=='Time' or columns[0]=='Time':
				cont=cont+1
				comp=False
			else:
				cont=cont+1
		else:
			cont=cont+1

	f.close()

	#read the data
	###########################################

	#read headers
	f = open(filename, 'r')
	for i in range(cont):
		header1 = f.readline()

	#read data
	for line in f:
		line = line.strip()
		line = line.replace("("," ")
		line = line.replace(")"," ")
		columns=" ".join(line.split())
		columns=columns.split(" ")
		#columns = line.split("\t")
		time.append(float(columns[0]))
		y.append(float(columns[2]))
		alpha.append(float(columns[9]))

	f.close()

	time=np.array(time)
	y=np.array(y)
	alpha=np.array(alpha)

	#find postions with time greater than time_ini
	pos=np.argwhere(time>=time_ini)
	if len(pos)>0:
		if pos.ndim>1:
			pos=pos[:,0]

		time=time[pos]
		y=y[pos]
		alpha=alpha[pos]
	else:
		time=[]
		y=[]
		alpha=[]

	return (time,y,alpha)

##################################################################################################################

##################################################################################################################

##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##############                                                                                      ##############
##############                                                                                      ##############
##############                 PROGRAMME EXECUTION                                                  ##############
##############                                                                                      ##############
##############                                                                                      ##############
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################

##################################################################################################################

##################################################################################################################

actualFolder=os.getcwd()

freeDecayCase=False

translationalSpring=True

folderFC,folderDisp=checkFolder(actualFolder,freeDecayCase)

#read the displacements file
dispFile=folderDisp+'/0/sixDoFRigidBodyState.dat'

if not os.path.isfile(dispFile):
	string="\n\nThere is no displacement file for 0 time folder inside the sixDoFRigidBodyState."
	string=string+"\nIf in your case study the time folder is different, please modify line 369 of the current script to adjust it to your case"
	sys.exit(string)

else:
	tDisp,y,alpha=read_displacements(dispFile,0.0)

if not freeDecayCase:
	#read force coefficients and displacements
	fcFile=folderFC+'/0/coefficient.dat'

	if not os.path.isfile(fcFile):
		string="\n\nThere is no force coefficient files for 0 time folder inside the forceCoeffsWing."
		string=string+'\nIf in your case study the time folder is different or the file name containing them is different from "coefficient.dat",'
		string=string+' please modify line 381 of the current script to adjust it to your case'
		sys.exit(string)

	tFC,Cd,Cl,Cm=read_forceCoeffs(fcFile,0.0)

if translationalSpring:
	disp=copy.deepcopy(y)

	if not freeDecayCase:
		force=copy.deepcopy(Cl)

else:
	disp=copy.deepcopy(alpha)

	if not freeDecayCase:
		force=copy.deepcopy(Cm)

C=1.0#[m] the chord length of the airfoil
U=100#[m/s] wind speed for the cases not in free decay

#############################################################################
###############General characteristics of the images

#que la fuente la renderice latex
thickness=0.75
plt.rcParams.update({
    "text.usetex": False,
    "font.size":9,
    'axes.linewidth':thickness,
    'xtick.major.width':0.9*thickness,
    'xtick.minor.width':0.75*thickness,
    'ytick.major.width':0.9*thickness,
    'ytick.minor.width':0.75*thickness,
    'xtick.major.size':2,
    'xtick.minor.size':1,
    'ytick.major.size':2,
    'ytick.minor.size':1,})

#controls the number of markers shown in the legend
plt.rcParams['legend.numpoints'] = 1
#controls the length of the line shown in the legend
plt.rcParams['legend.handlelength']=3
#controls the number of open windows
plt.rcParams['figure.max_open_warning']=30

#controls the size of the axis numbers
label_size = 7
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size 

fontSize=9
linew=0.5
markSize=4

#for contour plots
plt.rcParams["lines.linewidth"]=linew/2
plt.rcParams['contour.negative_linestyle'] = 'solid'


fig = plt.figure(constrained_layout=True)

horSizeCm=17#the horizontal dimension of the figure in cm

if freeDecayCase:
	verSizeCm=6.5#the vertical dimension of the figure in cm
	
	fig.set_size_inches(horSizeCm/2.54,verSizeCm/2.54)#first number is the width and second the height
	
	gs = gridspec.GridSpec(nrows=1,ncols=1,height_ratios=[1],width_ratios=[1],figure=fig)

else:
	verSizeCm=13#the vertical dimension of the figure in cm
	
	fig.set_size_inches(horSizeCm/2.54,verSizeCm/2.54)#first number is the width and second the height
	
	gs = gridspec.GridSpec(nrows=2,ncols=1,height_ratios=[0.5,0.5],width_ratios=[1],figure=fig)

########################################################################################

########################################################################################

ax0=fig.add_subplot(gs[0,0])
	
#plot the data
if freeDecayCase:
	if translationalSpring:
		l1 = ax0.plot(tDisp,disp,c='k',linestyle='-',linewidth=linew)
		labelY='x [m]'
	else:
		l1 = ax0.plot(tDisp,disp-disp[0],c='k',linestyle='-',linewidth=linew)
		labelY=u'\u03B1'+' ['+u'\N{DEGREE SIGN}'+']'
	labelX='t [s]'

else:
	if translationalSpring:
		l1 = ax0.plot(tDisp*U/C,disp/C,c='k',linestyle='-',linewidth=linew)
		labelY='x/C'
	else:
		l1 = ax0.plot(tDisp*U/C,disp-disp[0],c='k',linestyle='-',linewidth=linew)
		labelY=u'\u03B1'+' ['+u'\N{DEGREE SIGN}'+']'
	labelX='tU/C'

ax0.set_xlabel(labelX,fontsize=fontSize)
ax0.set_ylabel(labelY,fontsize=fontSize)

if freeDecayCase:
	ax0.set_xlim(0,10)
else:
	ax0.set_xlim(0,300)

ax0.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

########################################################################################

########################################################################################

if not freeDecayCase:
	ax1=fig.add_subplot(gs[1,0])
	
	l2 = ax1.plot(tFC*U/C,force,c='k',linestyle='-',linewidth=linew)

	if translationalSpring:
		labelY='Cl'
	else:
		labelY='Cm'

	labelX='tU/C'

	ax1.set_xlabel(labelX,fontsize=fontSize)
	ax1.set_ylabel(labelY,fontsize=fontSize)
	ax1.set_xlim(0,300)

########################################################################################

########################################################################################


#plt.savefig(pdf_name, dpi=100)
plt.savefig('image.pdf')
plt.savefig('image.png',dpi=1500)

plt.close()

######################################################################################################

######################################################################################################

######################################################################################################
