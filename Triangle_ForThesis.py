# Script to plot the PMF for a triangular graphene sealed microchamber
from scipy import *
from Constants import *
import Pseudo

# Geometry of the modeled device
theta=pi # angle of graphene armchair relative to x axis
dfiles=[\
	"./FEM/exx_50nm_2000psi_extrafine.dat",\
	"./FEM/exy_50nm_2000psi_extrafine.dat",\
	"./FEM/eyy_50nm_2000psi_extrafine.dat"]
lt=50*nano # Fundamental length of the device
P=2000*psi # Applied pressure
q=P*lt/(yg*tg) # Dimensionless loading parameter for thin plates

# Grid to be interpolated onto--formatted as image
dL=.005;
ygrid,xgrid=mgrid[-1/sqrt(3):1/sqrt(3):dL, -1/sqrt(3):1/sqrt(3):dL]

# Calculates the PMF
Bs=Pseudo.PMF(dfiles,lt,q,theta,xgrid,ygrid)
Pseudo.plotPMFsingle(xgrid,ygrid,Bs,save="./Triangle_PMF/ForThesis.png",show=True)