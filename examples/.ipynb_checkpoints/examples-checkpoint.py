from MDAnalysisData import datasets
adk = datasets.fetch_adk_equilibrium()
from ramachandran import *
u = Universe(adk['topology'],adk['trajectory'])
rp = RamachandranPlots(u)
rp.run()
adk
datasets.fetch_adk_equilibrium()
datasets.fetch_yiip_equilibrium_long()
YiPP = datasets.fetch_yiip_equilibrium_long()
u = Universe(YiPP['topology'],YiPP['trajectory'])
u
u.residues
rp = RamachandranPlots(u)
u = Universe(YiPP.topology,YiPP.trajectory)
u.bonds
u.atoms
YiPP.DESCR
print(YiPP.DESCR)
datasets.fetch_nhaa_equilibrium()
nhaa = datasets.fetch_nhaa_equilibrium()
u = Universe(nhaa.topology,nhaa.trajectory)
u
u.trajectory
rp = RamachandranPlots(u)
mpd  = datasets.fetch_membrane_peptide()
u = Universe(mpd.topology,mpd.trajectory)
rp = RamachandranPlots(u)
rp.run()
pwd
