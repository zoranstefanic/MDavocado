import os
import pickle
import pandas as pd
import numpy as np
import datashader as ds
import ruptures as rpt
from datashader.colors import viridis
from MDAnalysis import Universe
from MDAnalysis.analysis.dihedrals import Ramachandran, Janin
from string import ascii_uppercase
from shutil import copy
import json
from astropy.stats import circmean, circcorrcoef

class RamachandranPlots:
    def __init__(self, universe):
        "Takes MDAnalysis Universe as a input"
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

    def create_blanks(self):
        """Create two blank images for the first and last amino acid
        This relies on convert command from ImageMagick to be installed.
        """
        #FIXME: Could be implemented using PIL in pure Python.
        os.chdir(self.dir)
        n = self.nresidues
        os.system('convert -size 500x500 xc:white 1.png')
        os.system('convert -size 500x500 xc:white %d.png' % n)
        for i in range(self.n):
            print('Copying blanks to df%d' %i)           
            copy('1.png','df%d/1.png'%i) 
            copy('%d.png' %n,'df%d/%d.png'%(i,n)) 

    def annotate_images(self):
        "Annotate all images with their file numbers"
        #TODO: Could be implemented using PIL in pure Python.
        for i in range(self.n):
            print('Adding labels in dir: df%d' %i) 
            os.chdir('df%d' %i)
            os.system("mogrify -gravity South -annotate 0 '%t' -pointsize 50 *.png") 
            os.chdir(self.dir) 

    def create_gifs(self):
        """Create one gif for each amino acid stacking
        png images from directores df0 ... df9 i.e. times
        """
        if not os.path.exists('gifs'):
            print('Making directory gifs')
            os.mkdir('gifs')
        for i in range(1,self.nresidues+1):
            print('Making gif: %d' %i)  
            #TODO check if the command below respects the order df0 ... df9 (shell expends this)
            #                                  originally this was ./df{0..9}/
            command = "convert -dispose previous -delay 10 -loop 0 ./df*/%d.png gifs/%d.gif" %(i,i)
            os.system(command)
        os.chdir(self.dir) 
    
    def montage(self):
        """
        Combine all pngs into a tile for each chain.
        Finally combine all chain montages to a single gif.
        """
        chns = list(ascii_uppercase[:self.nchains]) # ['A','B','C'...,'F']
        N = self.res_per_chain 
        chains = dict(zip(chns,[(N*i+1,N*(i+1)) for i in range(self.nchains)]))
        for i in range(self.n):
            print('Making montage in dir: df%d' %i) 
            os.chdir('df%d' %i)
            for c,v in chains.items():
                pnglist = ' '.join([str(i)+'.png' for i in range(v[0],v[1])])
                os.system("montage -geometry 100x100 %s %s.png" %(pnglist,c))
                os.system("mogrify -gravity SouthEast -annotate 0 '%s' -pointsize 70 %s.png" % (c,c))
            os.chdir(self.dir) 
        for c in chns:
            flist = ' '.join(['./df'+str(i)+'/%s.png' %c for i in range(self.n)])
            os.system("convert -dispose previous -delay 10 -loop 0 %s %s.gif" % (flist,c))
        print('All done!')
    
    def run(self):
        self.make_all_dataframes()
        self.concat_dataframes()
        self.make_avokado_images()
        self.create_blanks()
        self.annotate_images()
        self.create_gifs()
        self.montage()

class JaninPlots:
    def __init__(self, universe):
        "Takes MDAnalysis Universe as a input"
        self.universe = universe
        self.nframes = universe.trajectory.n_frames
        self.protein = universe.select_atoms('protein and not resname ALA CYS* GLY PRO SER THR VAL')
        self.all_residues= universe.select_atoms('protein')
        self.resids = list(self.protein.residues.resids)
        self.nresidues = len(self.protein.residues.resids)
        self.nchains = len(self.protein.fragments)
        self.res_per_chain = self.nresidues//self.nchains
        self.resnames = self.get_resnames()
        self.all_resnames = self.get_all_resnames()
        self.num_res_per_chain_all = len(self.all_resnames)//self.nchains
        self.dir = os.path.dirname(self.universe.filename)
        self.n  = 10 # Number of pieces to cut the trajectory into
        self.delta = self.nframes // self.n
        self.cvs = ds.Canvas(plot_width=500, plot_height=500,x_range=(0,360),y_range=(0,360))
        self.partialdfs = []

    def get_resnames(self):
        return [r.resname for r in self.protein.residues]

    def get_all_resnames(self):
        return [r.resname for r in self.all_residues.residues]

    def n_to_residue(self,n):
        """Converts the serial number n (from 1 to nresidues)
        to a string "Res Ch num" e.g. PHE B 221
        """
        chains = list(ascii_uppercase)[:self.nchains]
        chain = chains[(n-1)//self.num_res_per_chain_all]
        resname = self.all_resnames[n-1]
        return u"%s %s %s" %(resname, chain, n%self.num_res_per_chain_all or self.num_res_per_chain_all)

    def slices(self):
        "A list of frame numbers to cut the trajectory"
        m, n = self.nframes, self.n
        l = list(range(0,m,m//n))
        if m%n: l[-1] = m
        else: l += [m]
        return l

    def make_partial_dataframe(self,start,end):
        os.chdir(self.dir)
        J = Janin(self.protein,verbose=True) 
        J.run(start,end)
        angles = J.results['angles'].reshape(end-start,self.nresidues*2)
        df = pd.DataFrame(angles)
        print('Writing out janin_%d_%d.pkl' %(start,end))
        df.to_pickle('janin_%d_%d.pkl' %(start,end))
        self.partialdfs.append('janin_%d_%d.pkl' %(start,end))
        del df

    def make_all_dataframes(self):
        frames = self.slices()
        for i in range(len(frames)-1):
           self.make_partial_dataframe(frames[i],frames[i+1]) 
    
    def concat_dataframes(self):
        df = pd.concat([pd.read_pickle(f) for f in self.partialdfs])
        df.reset_index(inplace=True,drop=True)
        df.astype('float32')
        df.to_pickle('janin_all.pkl')
        self.janin = pd.read_pickle('janin_all.pkl')

    def read_janin(self,file):
        self.janin = pd.read_pickle(file)
            
    def make_janin_images(self,aggregate=False):
        frames = self.slices()
        for i in range(len(frames)-1):
            directory = 'janin' + str(i)
            print(directory)
            if aggregate:
                start, stop = 0, frames[i+1]
            else:
                start, stop = frames[i], frames[i+1]
            dfw = self.janin.iloc[start:stop,:]
            self.create_images(dfw,directory)
            print('Finished directory %s' %directory)
            os.chdir(self.dir)

    def create_images(self,df,d):
        if not os.path.exists(d):
           os.mkdir(d)
        os.chdir(d)
        for n in self.resids:
           self.create_chi1_chi2(df,n)

    def create_chi1_chi2(self,df,n):
        "Creates the chi1_chi2 plot for residue number n during the trajectory"
        nn = self.resids.index(n) + 1
        phi_psi = df[[2*nn-2,2*nn-1]]
        phi_psi.columns = ['Chi1','Chi2']
        agg = self.cvs.points(phi_psi,'Chi1','Chi2')
        img = ds.tf.shade(agg, cmap=viridis)
        img = img.to_pil()
        img.save('%s.png' %(n))

    def create_blanks(self):
        """Create two blank images for the first and last amino acid
        This relies on convert command from ImageMagick to be installed.
        """
        #FIXME: Could be implemented using PIL in pure Python.
        os.chdir(self.dir)
        n = self.nresidues
        os.system('convert -size 500x500 xc:white 1.png')
        os.system('convert -size 500x500 xc:white %d.png' % n)
        for i in self.resids:
            print('Copying blanks to janin%d' %i)           
            copy('1.png','janin%d/1.png'%i) 
            copy('%d.png' %n,'janin%d/%d.png'%(i,n)) 

    def annotate_images(self):
        "Annotate all images with their file numbers"
        #TODO: Could be implemented using PIL in pure Python.
        for i in range(self.n):
            print('Adding labels in dir: janin%d' %i) 
            os.chdir('janin%d' %i)
            for j in self.resids:
                text = self.n_to_residue(j)
                os.system("mogrify -gravity South -annotate 0 '%s' -pointsize 50 %d.png" %(text,j))
            os.chdir(self.dir) 

    def create_gifs(self):
        """Create one gif for each amino acid stacking
        png images from directores df0 ... df9 i.e. times
        """
        if not os.path.exists('janin_gifs'):
            print('Making directory janin_gifs')
            os.mkdir('janin_gifs')
        for i in self.resids:
            print('Making gif: %d' %i)  
            #TODO check if the command below respects the order df0 ... df9 (shell expends this)
            #                                  originally this was ./df{0..9}/
            command = "convert -dispose previous -delay 10 -loop 0 ./janin*/%d.png janin_gifs/%d.gif" %(i,i)
            os.system(command)
        os.chdir(self.dir) 
    
    def montage(self):
        """
        Combine all pngs into a tile for each chain.
        Finally combine all chain montages to a single gif.
        """
        chns = list(ascii_uppercase[:self.nchains]) # ['A','B','C'...,'F']
        N = self.res_per_chain 
        chains = dict(zip(chns,[self.resids[N*i:N*(i+1)] for i in range(self.nchains)]))
        for i in range(self.n):
            print('Making montage in dir: janin%d' %i) 
            os.chdir('janin%d' %i)
            for c,v in chains.items():
                pnglist = ' '.join([str(i)+'.png' for i in v])
                os.system("montage -geometry 100x100 %s %s.png" %(pnglist,c))
                os.system("mogrify -gravity SouthEast -annotate 0 '%s' -pointsize 70 %s.png" % (c,c))
            os.chdir(self.dir) 
        for c in chns:
            flist = ' '.join(['./janin'+str(i)+'/%s.png' %c for i in range(self.n)])
            os.system("convert -dispose previous -delay 10 -loop 0 %s %s.gif" % (flist,c))
        print('All done!')

    def circular_correlation_matrix(self): 
        """After much work finally this function correctly calculates
        correlation coefficient between Ramachandran angles without 
        any need for removig spikes first, it uses original angles
        coming from MDanalysis.
        Uses this formula:
        https://docs.astropy.org/en/stable/api/astropy.stats.circcorrcoef.html
        Also cleverly uses numpy broadcasting to calculate fast!
        """
        l = []
        for i in self.resids:
            l += [2*i-2,2*i-1]
        self.janin.columns = l
        janin = np.radians(self.janin)
        cmean = circmean(janin)
        sinr = np.sin(janin - cmean)
        numerator = np.matmul(sinr.T,np.asarray(sinr))
        sinr2 = np.sum(np.sin(janin - cmean)**2)
        # denominator is np.sqrt od sin2_r_minus_mean[i] * sin2_r_minus_mean[j]
        # we have to get the matrix Aij from array Ai such that Aij = Ai x Aj.
        # Numpy outer product to the rescue:
        denominator = np.sqrt(np.outer(sinr2,sinr2))
        self.cm = numerator/denominator
        self.cm.columns = l

    def truncate_correlation_matrix(self,limit=0.5):
        "Leave only rows and columns which have at least one element with corr > limit or corr < - limit"
        # First set diagonal elements to 0 and the immediate neighbours 
        n = self.cm.shape[0]
        a = np.diag(np.ones(n-1),k=1) + np.diag(np.ones(n)) + np.diag(np.ones(n-1),k=-1) == True
        self.cm.mask(a,0,inplace=True)
        self.cm = self.cm[abs(self.cm)>limit]
        indices = self.cm.any(axis=0) 
        self.cm = self.cm.loc[indices,indices]

    def extract_correlations(self,limit=0.5):
        self.corr = {}
        cm_t = self.cm
        cm_t.columns= list(map(self.num_to_res_num,cm_t.columns))
        cm_t.index = list(map(self.num_to_res_num,cm_t.index))
        for c in cm_t.columns:
            self.corr[c] = dict(cm_t[cm_t[c].notnull()][c])
        pickle.dump(self.corr,open(os.path.join(self.dir,'janin_corr.pkl'),'wb'))

    def num_to_res_num(self,i):
        chi1_chi2 = ('chi1','chi2')[i%2]
        if chi1_chi2 == 'chi1':
            n = (i+2)/2
        else:
            n = (i+1)/2
        return int(n), chi1_chi2

    def run(self):
        self.make_all_dataframes()
        self.concat_dataframes()
        self.make_janin_images()
        #self.create_blanks()
        self.annotate_images()
        self.create_gifs()
        self.montage()
