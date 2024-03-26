import os
import pickle
import pandas as pd
import numpy as np
import datashader as ds
import ruptures as rpt
from datashader.colors import viridis
from MDAnalysis import Universe
from MDAnalysis.analysis.dihedrals import Ramachandran
from string import ascii_uppercase

class RamachandranPlots:
    def __init__(self, universe):
        self.universe = universe
        self.nframes = universe.trajectory.n_frames
        self.protein = universe.select_atoms('protein')
        self.nresidues = len(self.protein.residues.resids)
        self.nchains = len(self.protein.fragments)
        self.res_per_chain = self.nresidues//self.nchains
        self.resnames = self.get_resnames()
        self.dir = os.path.dirname(self.universe.filename)
        self.n  = 10 # Number of pieces to cut the trajectory into
        self.delta = self.nframes // self.n
        self.cvs = ds.Canvas(plot_width=500, plot_height=500,x_range=(-180,180),y_range=(-180,180))
        self.partialdfs = []

    def get_resnames(self):
        return [r.resname for r in self.protein.residues]

    def n_to_residue(self,n):
        """Converts the serial number n (from 1 to nresidues)
        to a string "Res Ch num" e.g. PHE B 221
        """
        chains = list(ascii_uppercase)[:self.nchains]
        chain = chains[(n-1)//self.res_per_chain]
        resname = self.resnames[n-1]
        return u"%s %s %s" %(resname, chain, n%self.res_per_chain or self.res_per_chain)

    def slices(self):
        "A list of frame numbers to cut the trajectory"
        m, n = self.nframes, self.n
        l = list(range(0,m,m//n))
        if m%n: l[-1] = m
        else: l += [m]
        return l

    def make_partial_dataframe(self,start,end):
        os.chdir(self.dir)
        R = Ramachandran(self.protein,verbose=True) 
        R.run(start,end)
        angles = R.results['angles'].reshape(end-start,(self.nresidues-2)*2)
        df = pd.DataFrame(angles)
        print('Writing out rama_%d_%d.pkl' %(start,end))
        df.to_pickle('rama_%d_%d.pkl' %(start,end))
        self.partialdfs.append('rama_%d_%d.pkl' %(start,end))
        del df

    def make_all_dataframes(self):
        frames = self.slices()
        for i in range(len(frames)-1):
           self.make_partial_dataframe(frames[i],frames[i+1]) 
    
    def concat_dataframes(self):
        df = pd.concat([pd.read_pickle(f) for f in self.partialdfs])
        df.reset_index(inplace=True,drop=True)
        df.astype('float32')
        df.to_pickle('rama_all.pkl')
        self.rama = pd.read_pickle('rama_all.pkl')

    def read_rama(self,file):
        self.rama = pd.read_pickle(file)
            
    def make_avokado_images(self,aggregate=False):
        frames = self.slices()
        for i in range(len(frames)-1):
            directory = 'df' + str(i)
            print(directory)
            if aggregate:
                start, stop = 0, frames[i+1]
            else:
                start, stop = frames[i], frames[i+1]
            dfw = self.rama.iloc[start:stop,:]
            self.create_images(dfw,directory)
            print('Finished directory %s' %directory)
            os.chdir(self.dir)

    def create_images(self,df,d):
        if not os.path.exists(d):
           os.mkdir(d)
        os.chdir(d)
        for n in range(2,self.nresidues):
           self.create_phi_psi(df,n)

    def create_phi_psi(self,df,n):
        "Creates the phi_psi plot for residue number n during the trajectory"
        phi_psi = df[[2*n-4,2*n-3]]
        phi_psi.columns = ['Phi','Psi']
        agg = self.cvs.points(phi_psi,'Phi','Psi')
        img = ds.tf.shade(agg, cmap=viridis)
        img = img.to_pil()
        img.save('%s.png' %(n))

    def run(self):
        self.make_all_dataframes()
        self.concat_dataframes()
        self.make_avokado_images()
