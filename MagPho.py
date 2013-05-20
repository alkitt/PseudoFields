# MagPho.py
# Module to hold formulas which deal with the magneto phonon coupling in graphene
# Follows M.O. Goerbig Rev Mod Phys 83, 1193 (2011)
from scipy import *		# scipy scripts
import Constants as C 	# fundamental constants and graphene constants

def lB(B):
	# The magnetic length in meters for a given magnetic field B (T)
	return sqrt(C.hbar/(C.ec*B))

def cycw(B):
	# The cyclotron frquency for a given magnetic field (T)
	return sqrt(2)*C.vg/lB(B)

def En(lam,n,B):
	# Returns the energy, in eV, with band index lambda (+/-1) of the nth Landau level at magnetic field B (T)
	if n<0:
		print "negative Landau level index passed to En--absolute value taken"
	if lam!=1 and lam!=-1:
		print "lambda must be set to -1 or +1 in En"

	energy=lam*C.hbar*cycw(B)*sqrt(n)/C.ec

	return energy

def fill(nel,B):
	# The filling factor for an electron density nel (1/m^2) and magnetic field B (T)
	nB=1./(2.*pi*lB(B)**2) # The magnetic flux density
	return nel/nB

def pfill(nel,lam,n,B):
	# Partial filling factor of the (lam, n) Landau level at magnetic field B (T) and an electron density nel (1/m^2).
	# Ranges from 0 to 1- 0 is empty, 1 is filled
	
	# First get the absolute filling factor
	nufil=fill(nel,B)
	# Convert filling factor to energy level (eq 104)
	nfild=1./4.*(nufil-2.)
	# If rounded toward zero, the sign gives the value of lambda and the absolute value gives the value of n
	# The absolute value of the decimal gives the partial filling of the nth level
	if nfild >=0:
		nfil=floor(nfild)
	else:
		nfil=ceil(nfild)

	if nel*lam<nfil: 	# If the level is less than the fill, it is filled
		pf=1
	elif nel*lam>nfil:	# If the level is greater than the fill, it is empty
		pf=0
	else:				# If they are the same level, use the decimal point
		pf=abs(nfil%1)

	return pf

def kron(a,b):
	if a==b:
		return 1
	else:
		return 0

def gRH(nel,n,B):
	# Effective coupling constants (eqn 189) for doping nel (1/m^2), the level n, and magnetic field B (T)
	gamma=3.*sqrt(3)*C.ag**2/(2.*pi*lB(B))
	return C.gep*sqrt((1+kron(n,0))*gamma)*sqrt(pfil(nel,-1,n+1,B)-pfil(nel,+1,n,B))

def gLH(nel,n,B):
	# Effective coupling constants (eqn 189) for doping nel (1/m^2), the level n, and magnetic field B (T)
	gamma=3.*sqrt(3)*C.ag**2/(2.*pi*lB(B))
	return C.gep*sqrt((1+kron(n,0))*gamma)*sqrt(pfil(nel,-1,n,B)-pfil(nel,+1,n+1,B))

def MEinter(n,B):
	# The energy of the magneto exciton (eV) corresponding to the transitions between -,n and +,n+1 at magnetic field B (T)
	return C.hbar*cycw(B)*(sqrt(n+1)+sqrt(n))/C.ec

def MEintra(n,B):
	# The energy of the magneto exciton (eV) corresponding to the transitions between +,n and +,n+1 at magnetic field B (T)
	return C.hbar*cycw(B)*(sqrt(n+1)-sqrt(n))/C.ec
