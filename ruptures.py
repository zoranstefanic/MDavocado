import os
import pickle
import pandas as pd
import numpy as np
import ruptures as rpt

class RupturePlots:
    """
    A class used for change point detection in Ramachandran
    angles time series. Uses ruptures library.
    """

    def __init__(self,file):
        "Reads pickled Pandas DataFrame file with Ramachandran angles"
        self.rama = pd.read_pickle(file)
        self.nframes = self.rama.shape[0]
        self.skip = 100
        self.penalty = 1e4
        self.width = 100
        self.model = "l2"
        self.n_bkps = 4
        self.change_points = {}

    def phi_psi_for_n(self,n):
        "Returns column with phi and psi for a given amino acid number"
        angles = self.rama.iloc[:self.nframes:self.skip,2*n-4:2*n-2]
        return angles
    
    def detect(self,n):
        "Detects the change points"
        print(n)
        angles = self.phi_psi_for_n(n)
        algo = rpt.Window(width=self.width, model=self.model).fit(angles)
        #result = algo.predict(n_bkps=n_bkps) # Fixed number of breakpoints
        result = algo.predict(pen=self.penalty)
        return angles, result

    def remove_spikes(self):
        d = self.rama.diff()
        d[d>180]-= 360
        d[d<-180]+= 360
        d.iloc[0,:] = self.rama.iloc[0,:]
        d = d.cumsum()
        # This code has the consequence of blowing angles for some residues
        # which leads to nonsense correlations latter and also to false breakpoints
        # As a remedy the final Ramachandran is cut mod 360, so that it can not exceed
        # +360 or go below -360. For many angles it will not change a thing but for some 
        # it will mean a lot of difference, for sure for terminal residues between chains
        # FIXME
        # It seems that this line below only worsens the stuff for some residues like this:
        # https://alokomp.irb.hr/md/correlations_circular/1714/1325
        # https://alokomp.irb.hr/md/correlations_circular/1697/1184
        d = d%(np.sign(d)*360)
        pickle.dump(d, open('rama_all_no_spikes.pkl','wb'))

    def save_png(self,n):
        angles, result = self.detect(n)
        self.change_points[n] = result
        fig, a = rpt.display(angles,result)
        phi_ax, psi_ax = a
        phi, psi = angles.iloc[:,0], angles.iloc[:,1]
        max_phi, max_psi = max([phi.max(),180]), max([psi.max(),180])
        min_phi, min_psi = min([phi.min(),-180]), min([psi.min(),-180])
        print(min_phi,max_phi,min_psi,max_psi)
        phi_ax.set_ylim([min_phi,max_phi])
        psi_ax.set_ylim([min_psi,max_psi])
        #phi_psi = ('\u03c6','\u03c8')[i%2]
        phi_ax.set_title(num_to_res_n(n) + ' \u03c6', pad=-10,loc="left")
        psi_ax.set_title(num_to_res_n(n) + ' \u03c8', pad=-10,loc="left")
        fig.savefig('%d_.png' %n)
    
    def make_figures(self):
        for n in range(2,1397):
            self.save_png(n)
