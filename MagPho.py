# MagPho.py
# Module to hold formulas which deal with the magneto phonon coupling in graphene
from scipy import *
import Constants as C

def En(n,B):
	# Returns the energy, in eV, of the nth Landau level at magnetic field B (T)
	energy=.0324*copysign(1,n)*sqrt(abs(n)*B)

	return energy

def fermi(s,n,B,ndoping,T):
	# 
	fermilevel=1./(exp(s*En(n,B)-fermilevel(ndoping,B)))

	return fermilevel