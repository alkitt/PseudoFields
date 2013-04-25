# Script to plot the PMF for a triangular graphene sealed microchamber
from scipy import *
from Constants import *
import Pseudo

def mask(xgrid,ygrid,Bs,h,lt,alpha,theta):
	# Function which masks out the pixels which are not over the microchamber
	# Only have to use for concave shapes like the bowtie
	# h the total extent of the bowie in y
	# lt the width of the gap
	# theta the angle to the bowtie was rotated by

	for i in range(xgrid.shape[1]):
		for j in range(xgrid.shape[0]):
			# Rotate the x y values back to their original orientation
			xr=cos(-theta)*xgrid[0,i]-sin(-theta)*ygrid[j,0]
			yr=sin(-theta)*xgrid[0,i]+cos(-theta)*ygrid[j,0]
			# Check if they are over the hole
			if \
			yr>h/2 or \
			yr<-h/2 or \
			(yr<= tan(alpha)*(xr-1./2.) and yr>=-tan(alpha)*(xr-1./2.)) or \
			(yr<=-tan(alpha)*(xr+1./2.) and yr>= tan(alpha)*(xr+1./2.)):#Above, Below, Right Wedge, Left Wedge
				Bs[j,i]=nan

	return Bs

for i in arange(0,181,2.5):
	# Geometry of the modeled device
	theta=i*pi/180. # angle of graphene armchair relative to x axis
	dfiles=[\
		"./FEM/exx_bowt_20nm_2000psi_extrafine.dat",\
		"./FEM/exy_bowt_20nm_2000psi_extrafine.dat",\
		"./FEM/eyy_bowt_20nm_2000psi_extrafine.dat"]
	lt=20*nano # Fundamental length of the device (width of gap)
	h=100*nano/lt # Total extent of unrotaed bowtie in y normalized
	alpha=45*pi/180 # The half angle of the bowtie
	P=2000*psi # Applied pressure
	q=P*lt/(yg*tg) # Dimensionless loading parameter for thin plates

	# Grid to be interpolated onto--formatted as image
	dL=.02
	ext=4
	ygrid,xgrid=mgrid[-ext:ext:dL, -ext:ext:dL]

	# Calculates the PMF
	Bs=Pseudo.PMF(dfiles,lt,q,theta,xgrid,ygrid)
	
	#Masks out the data points not over the microchamber
	Bs=mask(xgrid,ygrid,Bs,h,lt,alpha,theta)

	Pseudo.plotPMF(xgrid,ygrid,Bs,save="./Bowtie_PMF/Rotate_"+str(int(i))+".png",show=False,limits=[-35,35])