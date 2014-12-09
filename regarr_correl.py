import numpy as np
import numpy.ma as ma
import iris as iris
import iris.plot as iplt
import iris.quickplot as qplt
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import mycmaps as mc
import scipy.stats as stats
import sys
import troposave as ta
import prettyplotlib as ppl
import pickle

with open('./pickles/regarr.pickle','rb') as f:
    regarr = pickle.load(f)
# for n in xrange(regarr.shape[0]):
#     regarr[n,:] = regarr[n,:]/(regarr[n,:].std())

pval_reg = np.zeros((regarr.shape[1],regarr.shape[1]))
for n in xrange(regarr.shape[1]):
    for m in xrange(regarr.shape[1]):
#         if n<=m:
        pval_reg[n,m] = np.corrcoef(regarr[0:-4,n],regarr[0:-4,m])[0,1]

var = np.array(['Tsfc','smc clim','T700hPa','T300hpa','RH700hPa','RH300hPa','DLWR','DSWR','smc','Cld High','Cld Low','Precip','u_high','u_low','v_high','v_low'])
regs =np.array(['India','MC','TropSthAm','SthSthAm','NthWestAfr','NthEastAfr','TropAfr','SthAfr','Aus','lat10', 'latN20', 'latS20', 'latN30', 'latS30'])

plt.close('all')
fig, axes = plt.subplots(nrows=1) #,ncols=2)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
ppl.pcolormesh(fig, axes, pval_reg, ax_colorbar=cbar_ax, yticklabels=regs, xticklabels=regs,vmin=-1.0,vmax=1.0)
axes.set_xlabel('Regions')
axes.set_ylabel('Regions')
axes.set_title('Correlation between regions for range of variables')
plt.show()
fig.set_size_inches(18,5)
plt.savefig('regresp_cor_manyvar.eps')

