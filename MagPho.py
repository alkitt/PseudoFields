# MagPho.py
# Module to hold formulas which deal with the magneto phonon coupling in graphene
# Follows M.O. Goerbig Rev Mod Phys 83, 1193 (2011)
from scipy import *					# scipy scripts
import Constants as C 				# fundamental constants and graphene constants
from scipy.optimize import root		# Numerically solve for roots of a function

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
	# Determine the heighest filled (lam,n) and get partial filling factor
	# Be careful, \nu -2 to 2 is n =0, then goes by \nu+4 increments
	if abs(nufil)<=2:								# n=0
		if lam<0 and n>=1:							# Level full
			pf=1
		elif lam>0 and n>=1:						# Level empty
			pf=0
		elif n==0:									# Partial fill
			pf=(nufil+2.)/4.
		else:
			print "Error 0 in pfill function"
			return nan
	elif abs(nufil)>2:								# n>=1
		# Determine (lamt,nt) for given nel
		lamt=copysign(1,nufil)
		if lamt>0:
			nt=1+floor((abs(nufil)-2.)/4.)
		elif lamt<0:
			nt=  ceil((abs(nufil)-2.)/4.)
		else:
			print "Error 1 in pfill function"
			return nan
		# Compare (lam,n) to (lamt,nt)
		if   lam*n< lamt*nt:						# Level full
			pf=1
		elif lam*n> lamt*nt:						# Level empty
			pf=0
		elif lam*n==lamt*nt:						# Partial fill
			if lamt>0:
				pf=((abs(nufil)-2.)%4)/4.
			elif lamt<0:
				pf=1-((abs(nufil)-2.)%4)/4.
			else:
				print "Error 4 in pfill function"
				return nan
		else:
			print "Error 2 in pfill function"
			return nan
	else:
		print "Error 3 in pfill function"
		return nan

	return pf

def kron(a,b):
	if a==b:
		return 1
	else:
		return 0

def g(A,nel,n,B):
	# Effective coupling constants (eqn 189) for doping nel (1/m^2), the level n, magnetic field B (T), and handedness A
	# A=+1 for RH, A=-1 for LH
	gamma=3.*sqrt(3)*C.ag**2/(2.*pi*lB(B)**2)
	if A==+1:
		Aout=C.gep*sqrt((1+kron(n,0))*gamma)*sqrt(pfill(nel,-1,n+1,B)-pfill(nel,+1, n ,B))
	elif A==-1:
		Aout=C.gep*sqrt((1+kron(n,0))*gamma)*sqrt(pfill(nel,-1, n ,B)-pfill(nel,+1,n+1,B))
	else:
		print "Error in g in MagPho.py"
		return nan
	return Aout

def MEinter(n,B):
	# The energy of the magneto exciton (eV) corresponding to the transitions between -,n and +,n+1 at magnetic field B (T)
	return C.hbar*cycw(B)*(sqrt(n+1)+sqrt(n))/C.ec+C.dLL

def MEintra(n,B):
	# The energy of the magneto exciton (eV) corresponding to the transitions between +,n and +,n+1 at magnetic field B (T)
	return C.hbar*cycw(B)*(sqrt(n+1)-sqrt(n))/C.ec+C.dLL

def eq195(eA,A,nel,B,Nc):
	# The roots of EQ 195 give the complex energy of the magneto-phonon mixed states, eA=[Re[eA],Im[eA]] (eV).
	# B magnetic field (T), nel electron doping in (1/m**2), Nc heighest LL to sum to
	# A handedness A=+1 for RH, A=-1 for LH
	ceA=complex(eA[0],eA[1])
	LHS=ceA**2-C.w0**2 													#LHS of eqn 195
	# Building the RHS of eqn 195
	Nf=0 # Restricting the sum would save some computation power (if necessary)
	inter=sum(MEinter(n ,B)*g(A,nel,n ,B)**2/(ceA**2-MEinter(n, B)**2) for n in xrange(Nf+1,Nc+1,1))
	intra=    MEintra(Nf,B)*g(A,nel,Nf,B)**2/(ceA**2-MEintra(Nf,B)**2)
	RHS=4*C.w0*(inter+intra) 											#RHS of eqn 195

	# Find root by finding where real and imaginary part go to zero
	metric=LHS-RHS
	return (metric.real,metric.imag)

def Seq195(stpoint,A,nel,B,Nc):
	# Solves for a root of eq195() near stpoint=[Real(G band energy),Imaginary(G band energy)] in eV
	# Return real and complex parts of G band energy in eV
	sol1=root(lambda x:eq195(x,A,nel,B,Nc),stpoint,jac=False)
	if sol1.success==True:
		out=sol1.x
	elif sol1.success==False:
		print sol1.message
		out=nan
	else:
		print "Error Seq195 solve"
		return nan

	return out

def MSeq195(A,nel,B,Nc):
	# Uses Seq195 to solve for each of the roots of eq195
	# Look on each side of each resonance and then eliminate duplicates
	# Expect a total 1 root for G band, 1 root for each magnetoexciton (inter and intra)
	roots=zeros((2*(1+Nc+1),2))
	# Look on both sides of each resonance by a factor of epsilon
	eps=.001
	roots[0]        =Seq195([(1+eps)*C.w0.real        ,C.w0.imag]        ,A,nel,B,Nc)		# G band
	roots[1]	    =Seq195([(1-eps)*C.w0.real        ,C.w0.imag]        ,A,nel,B,Nc)
	for n in range(Nc):																		# inter band transitions
		roots[2+n*2]=Seq195([(1+eps)*MEinter(n,B).real,MEinter(n,B).imag],A,nel,B,Nc)
		roots[3+n*2]=Seq195([(1-eps)*MEinter(n,B).real,MEinter(n,B).imag],A,nel,B,Nc)
	roots[2*Nc+2]   =Seq195([(1+eps)*MEintra(n,B).real,MEintra(n,B).imag],A,nel,B,Nc)		# intra band transitions
	roots[2*Nc+3]   =Seq195([(1-eps)*MEintra(n,B).real,MEintra(n,B).imag],A,nel,B,Nc)

	return roots