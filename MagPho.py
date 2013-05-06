# MagPho.py
# Module to hold formulas which deal with the magneto phonon coupling in graphene
from scipy import *
import Constants as C

def En(n,B):
	# Returns the energy, in eV, of the nth Landau level at magnetic field B (T)
	# negative n gives negative energy
	energy=sqrt(2*C.hbar*C.ec*C.vg**2)*copysign(1,n)*sqrt(abs(n)*B)/C.ec

	return energy

def fermi(s,n,B,ndoping,T):
	# 
	fermilevel=1./(exp(s*En(n,B)-fermilevel(ndoping,B)))

	return fermilevel