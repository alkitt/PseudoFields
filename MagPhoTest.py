# Trouble shooting MagPho.py by plotting some functions
import MagPho as MP 				# Scripts to be tested
import Constants as C 				# Fundamental and graphene constants
from scipy import * 				# scipy functionality
import matplotlib.pyplot as plt 	# Plotting functionality

# Reproducing the dotted lines on Fig 22 (a) in Goerbig Rev Mod Phys 83, 1223 (2011)----------
Bmin=3; Bmax=45;
B=arange(Bmin,Bmax,.2) # In Tesla

fig=plt.figure(1)
ax=fig.add_subplot(111,aspect=6)
ax.plot(B,zeros(B.shape[0])+1,'k:')
for n in range(0,3,1):
	ax.plot(B,MP.MEinter(n,B)/C.w0,'k:')
ax.set_xlabel('Magnetic Field, B (T)')
ax.set_ylabel('Energy in units of $\omega_{ph}$')
ax.axis([Bmin,Bmax,0,2])
# plt.show()

# Testing the partial filling factor code-----------------------------------------------------
nel=1e16 # The electron density in 1/m^2

# First consider the filling factor as magnetic field is tuned
fig=plt.figure(2)
fig.clear()
ax=fig.add_subplot(321)
ax.plot(B,MP.fill(nel,B))
ax.set_ylabel(r'$\nu$')
ax.set_xlim(Bmin,Bmax)

# As compared to the partial filling factor
# Lambda positive different n and different B
ax=fig.add_subplot(322)
pstyles=['k-','r--','b:','g-.','b-.']
pstyle=0
for n in [0,1,2,3,-1]:
	pfplot=zeros(B.shape[0])
	for i in range(0,B.shape[0],1):
		pfplot[i]=MP.pfill(nel,copysign(1,n),abs(n),B[i])
	ax.plot(B,pfplot,pstyles[pstyle%len(pstyles)])
	pstyle+=1
ax.set_ylabel('Partial filling factor')
ax.set_xlim(Bmin,Bmax)
ax.set_ylim(-.1,1.1)

nel=.5e16 # The electron density in 1/m^2

# First consider the filling factor as magnetic field is tuned
fig=plt.figure(2)
ax=fig.add_subplot(323)
ax.plot(B,MP.fill(nel,B))
ax.set_ylabel(r'$\nu$')
ax.set_xlim(Bmin,Bmax)

# As compared to the partial filling factor
ax=fig.add_subplot(324)
for n in [0,1,2,3,-1]:
	pfplot=zeros(B.shape[0])
	for i in range(0,B.shape[0],1):
		pfplot[i]=MP.pfill(nel,copysign(1,n),abs(n),B[i])
	ax.plot(B,pfplot,pstyles[pstyle%len(pstyles)])
	pstyle+=1
ax.set_ylabel('Partial filling factor')
ax.set_xlim(Bmin,Bmax)
ax.set_ylim(-.1,1.1)

nel=-1e16 # The electron density in 1/m^2

# First consider the filling factor as magnetic field is tuned
fig=plt.figure(2)
ax=fig.add_subplot(325)
# ax.plot(B,MP.fill(nel,B))
ax.plot(B,MP.fill(nel,B))
ax.set_ylabel(r'$\nu$')
ax.set_xlim(Bmin,Bmax)

# As compared to the partial filling factor
ax=fig.add_subplot(326)
for n in [0,-1,-2,-3,1]:
	pfplot=zeros(B.shape[0])
	for i in range(0,B.shape[0],1):
		pfplot[i]=MP.pfill(nel,copysign(1,n),abs(n),B[i])
	ax.plot(B,pfplot,pstyles[pstyle%len(pstyles)])
	pstyle+=1
ax.set_ylabel('Partial filling factor')
ax.set_xlim(Bmin,Bmax)
ax.set_ylim(-.1,1.1)

# Now test the effective coupling constants (eqn 189)-----------------------------------------
# fig=plt.figure(4)
# fig.clear()
# ax=fig.add_subplot(121)
# pstyle=0
# for n in [0,1,2,3]:
# 	ax.plot(B,sqrt((1+MP.kron(n,0))*3.*sqrt(3)*C.ag**2/(2.*pi*MP.lB(B)**2)),pstyles[pstyle%len(pstyles)])
# 	pstyle+=1
# ax=fig.add_subplot(122)
# pstyle=0
# for n in [0,-1,-2,-3]:
# 	ax.plot(B,sqrt(1/(MP.lB(B)**2)),pstyles[pstyle%len(pstyles)])
# 	pstyle+=1

fig=plt.figure(3)
fig.clear()
ax=fig.add_subplot(321)
pstyle=0
for n in [0,1,2,3]:
	pfplot=zeros(B.shape[0])
	for i in range(0,B.shape[0],1):
		pfplot[i]=MP.g(1,1e16,n,B[i])
	ax.plot(B,pfplot,pstyles[pstyle%len(pstyles)])
	pstyle+=1
ax=fig.add_subplot(322)
pstyle=0
for n in [0,1,2,3]:
	pfplot=zeros(B.shape[0])
	for i in range(0,B.shape[0],1):
		pfplot[i]=MP.g(-1,1e16,n,B[i])
	ax.plot(B,pfplot,pstyles[pstyle%len(pstyles)])
	pstyle+=1
ax=fig.add_subplot(323)
pstyle=0
for n in [0,1,2,3]:
	pfplot=zeros(B.shape[0])
	for i in range(0,B.shape[0],1):
		pfplot[i]=MP.g(1,0,n,B[i])
	ax.plot(B,pfplot,pstyles[pstyle%len(pstyles)])
	pstyle+=1
ax=fig.add_subplot(324)
pstyle=0
for n in [0,1,2,3]:
	pfplot=zeros(B.shape[0])
	for i in range(0,B.shape[0],1):
		pfplot[i]=MP.g(-1,0,n,B[i])
	ax.plot(B,pfplot,pstyles[pstyle%len(pstyles)])
	pstyle+=1
ax=fig.add_subplot(325)
pstyle=0
for n in [0,1,2,3]:
	pfplot=zeros(B.shape[0])
	for i in range(0,B.shape[0],1):
		pfplot[i]=MP.g(1,-1e16,n,B[i])
	ax.plot(B,pfplot,pstyles[pstyle%len(pstyles)])
	pstyle+=1
ax=fig.add_subplot(326)
pstyle=0
for n in [0,1,2,3]:
	pfplot=zeros(B.shape[0])
	for i in range(0,B.shape[0],1):
		pfplot[i]=MP.g(-1,-1e16,n,B[i])
	ax.plot(B,pfplot,pstyles[pstyle%len(pstyles)])
	pstyle+=1


plt.show()