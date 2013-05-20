from math import pi

# Metric prefixes
mili=1e-3
micro=1e-6
nano=1e-9
kilo=1e3
Mega=1e6
Giga=1e9
Tera=1e12

# Fundamental Constants--stick with SI
kb=1.3806488e-23 		# Boltzmans constant in J/K
ec=1.602176565e-19 		# Charge on the electron in Coulomb
hbar=1.05457173e-34 	# hbar in Js

# Derived constants
fq=hbar*2*pi/ec 		# Flux quanta in T*m^2

# Unit conversions--all to SI units
psi=6894.75729

# Constants for graphene
yg=Tera 				# Young's modulus of graphene (Pa)
tg=.34*nano 			# Thickness of graphene (m)
ag=.142*nano 			# Nearest neighbor distance (m)
betag=3.37 				# Electron phonon coupling
thop=2.8*ec 			# Nearest neighbor hopping energy (J)
vg=3./2.*ag*thop/hbar 	# Fermi velocity (m/s)
gep=.26					# Graphene electron phonon coupling (eV)
w0=.2 					# G band phonon energy (eV)