# Second program to test MagPho.py
import MagPho as MP 				# Scripts to be tested
import Constants as C 				# Fundamental and graphene constants
from scipy import * 				# scipy functionality
import matplotlib.pyplot as plt 	# Plotting functionality
from scipy.optimize import root		# Numerically solve for roots of a function

Bset=6.4

eAtest=arange(.01,.36,.0001)
metric=zeros(eAtest.shape[0])
for i,eA in enumerate(eAtest):
	(metric[i],junk)=MP.eq195([eA,.01],1,0,Bset,4)
fig=plt.figure(1)
ax=fig.add_subplot(111)
ax.plot(eAtest,metric)
ax.grid(True)
plt.show()

# Bs=arange(1,30,.2)
# ws=zeros(Bs.shape[0])
# Gs=zeros(Bs.shape[0])
# for i,B in enumerate(Bs):
# 	(ws[i],Gs[i])=MP.Seq195([C.w0.real,C.w0.imag],1,0,B,4)
# fig=plt.figure(2)
# ax=fig.add_subplot(121)
# ax.plot(Bs,ws)
# ax=fig.add_subplot(122)
# ax.plot(Bs,Gs)
# plt.show()

print "n=0", MP.MEinter(0,Bset)
print "n=1", MP.MEinter(1,Bset)
print "n=2", MP.MEinter(2,Bset)
print "n=3", MP.MEinter(3,Bset)
print "intra", MP.MEintra(1,Bset)



print MP.Seq195([C.w0.real,C.w0.imag],1,0,Bset,4)

print MP.Seq195([.3105,.01],1,0,Bset,4)

bob=MP.MSeq195(1,0,Bset,4)
print bob

newlist=array([v for v in bob if len(set(v))==len(v)])

print newlist