# Trouble shooting MagPho.py by plotting some functions
import MagPho as MP 				# Scripts to be tested
import Constants as C 				# Fundamental and graphene constants
from scipy import * 				# scipy functionality
import matplotlib.pyplot as plt 	# Plotting functionality

# Reproducing the dotted lines on Fig 22 (a) in Goerbig Rev Mod Phys 83, 1223 (2011)----------
Bmin=3; Bmax=45;
B=arange(Bmin,Bmax,.2) # In Tesla

# fig=plt.figure(1)
# ax=fig.add_subplot(111,aspect=6)
# ax.plot(B,zeros(B.shape[0])+1,'k--')
# for n in range(0,3,1):
# 	ax.plot(B,MP.MEinter(n,B)/C.w0,'k--')
# ax.set_xlabel('Magnetic Field, B (T)')
# ax.set_ylabel('Energy in units of $\omega_{ph}$')
# ax.axis([Bmin,Bmax,0,2])
# # plt.show()

# Testing the partial filling factor code-----------------------------------------------------
nel=1e16 # The electron density in 1/m^2

# First consider the filling factor as magnetic field is tuned
fig=plt.figure(2)
ax=fig.add_subplot(321)
ax.plot(B,MP.fill(nel,B))
ax.set_ylabel(r'$\nu$')
ax.set_xlim(Bmin,Bmax)

# Intermediate manipulations toward Partial filling factor
nfild=1./4.*(MP.fill(nel,B)-2.)
ax=fig.add_subplot(323)
ax.plot(B,nfild)
ax.set_ylabel(r'$\frac{1}{4} (\nu-2)$')
ax.set_xlim(Bmin,Bmax)

nfil=zeros(nfild.shape[0])
for ni in range(nfild.shape[0]):
	# If rounded toward zero, the sign gives the value of lambda and the absolute value gives the value of n
	# The absolute value of the decimal gives the partial filling of the nth level
	if nfild[ni] >=0:
		nfil[ni]=floor(nfild[ni])
	else:
		nfil[ni]=ceil(nfild[ni])
ax=fig.add_subplot(325)
ax.plot(B,nfil)
ax.set_xlabel('Magnetic Field, B (T)')
ax.set_ylabel(r'$[\frac{1}{4} (\nu-2)]$')
ax.set_xlim(Bmin,Bmax)
plt.show()

# Finally checking the partial filling factors
fig=plt.figure(3)
ax=fig.add_subplot(111,aspect=10)
pstyles=['k-','r--','b:','g-.']
pstyle=0
for n in range(-3,3,1):
	pfplot=zeros(B.shape[0])
	for i in range(0,B.shape[0],1):
		pfplot[i]=MP.pfill(nel,copysign(1,n),abs(n),B[i])
	ax.plot(B,pfplot,pstyles[pstyle%len(pstyles)])
	pstyle+=1
ax.set_xlabel('Magnetic Field, B (T)')
ax.set_ylabel('Partial filling factor')
ax.axis([Bmin,Bmax,-.1,1.1])
# plt.show()