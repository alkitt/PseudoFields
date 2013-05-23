# Second program to test MagPho.py
import MagPho as MP 				# Scripts to be tested
import Constants as C 				# Fundamental and graphene constants
from scipy import * 				# scipy functionality
import matplotlib.pyplot as plt 	# Plotting functionality
from scipy.optimize import root		# Numerically solve for roots of a function

Bset=6.33892

eAtest=arange(.18,.22,.0001)
metric=zeros(eAtest.shape[0])
for i,eA in enumerate(eAtest):
	(metric[i],junk)=MP.eq195([eA,.00856078],1,0,Bset,4)
fig=plt.figure(1)
ax=fig.add_subplot(111)
ax.plot(eAtest,metric)
ax.grid(True)
plt.show()
bob=MP.MSeq195(1,0,Bset,4)
print 'All zeros \n',bob
print 'No duplicates \n',MP.nodupl(bob)

# Bs=arange(5.8,6.8,.01)
# ws=zeros((5*Bs.shape[0],2))
# Gs=zeros((5*Bs.shape[0],2))
# count=0
# for B in Bs:
	# bob=MP.nodupl(MP.MSeq195(1,0,B,2))
	# for bobt in bob:
		# ws[count,0]=B; Gs[count,0]=B;
		# ws[count,1]=bobt[0]
		# Gs[count,1]=bobt[1]
		# count+=1
# fig=plt.figure(2)
# ax=fig.add_subplot(121)
# ax.plot(ws[:,0],ws[:,1],'k.')
# ax=fig.add_subplot(122)
# ax.plot(Gs[:,0],Gs[:,1],'k.')
# plt.show()