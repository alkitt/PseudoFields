# Raman.py
# Module to hold often used Raman fuctions
from scipy import *
import Constants as C

def Lor(k,wc,ga):
	# The value of a normalized Lorentzian with center frequency wc and width ga evaluated at k
	spec=1./pi*(ga/2.)/((k-wc)**2.+(1./2.*ga)**2.)

	return spec