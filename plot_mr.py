# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import pickle
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from pylab import *
from astropy import constants as const

mearth = 5.9721986E27 #grams
rearth = 6.3675E8 #centimeters
msun = 1.988435E33 #grams
rsun = 6.955E10 #centimeters
cmperau=1.496E13
lsun = 3.846E33 #ergs per second
tsun = 5780 #kelvins
ss_names = {"name" : ['Me','V','E','Ma','J','S','U','N']}
ss = pd.DataFrame(ss_names) 
ss['m'] = [0.0553, 0.815, 1.00, 0.107, 317.83, 95.159, 14.536, 17.147]
ss['um'] = ss.m*0.10
ss['r'] = [0.383,0.950,1.0, 0.532, 10.973, 9.140, 3.981, 3.865]
ss['ur'] = ss.r*0.05
ss['a'] = [0.387, 0.723 , 1.0, 1.524, 5.204, 9.582, 19.201, 30.047]
ss['rho'] = (ss.m * mearth)/(4./3 * np.pi * (ss.r * rearth)**3.)
ss['urho'] = np.sqrt((ss.um/ss.m)**2. + (3*ss.ur/ss.r)**2.)*ss.rho
ss['flux'] = 1./ss.a**2
ss['per'] = 365.25 * pow(ss.a, 3./2)
ss['fe'] = [0.0]*8
ss['age'] = [5.0]*8
ss = ss.set_index('name')
ss.head()

df = pickle.load(open('exoplanets.pickle'))
medf = np.median(df.flux)
df_hot = df[df.flux > medf]
df_cold = df[df.flux <= medf]
sp = pickle.load(open('small_exoplanets.pickle'))
#sp = sp.ix[sp.index!='Kepler-78 b',:]
sp['evap'] = sp.flux * sp.age
#sp.um = np.maximum(sp.um,0.1*sp.m)
#sp.ur = np.maximum(sp.ur,0.05*sp.r)
#sp.urho = np.maximum(abs(sp.urho), (np.sqrt((sp.um/sp.m)**2. + (3*sp.ur/sp.r)**2.)*sp.rho))
sp_clean = sp[sp.urho < 6.5]
sp_1sig = sp[sp.m > sp.um]
sp_cold = sp[sp.flux <= np.median(sp.flux)]
sp_hot = sp[sp.flux > np.median(sp.flux)]
sp['ttv'] = [0]*len(sp)
sp.ix[['Kepler-11 b','Kepler-11 c','Kepler-11 d','Kepler-11 f','Kepler-30 b','Kepler-36 b','Kepler-36 c','KOI-152 b','KOI-152 c','KOI-152 e'],'ttv'] = 1
ttvs = sp.ix[['Kepler-11 b','Kepler-11 c','Kepler-11 d','Kepler-11 f','Kepler-30 b','Kepler-36 b','Kepler-36 c','KOI-152 b','KOI-152 c','KOI-152 e'],:]
sp_vol = sp[sp.r >= 1.5]
sp_rocky = sp[sp.r < 1.5].append(ss[ss.r < 1.5])
sp_rocky.ix[:,['r','ur','m','um','rho','urho']]
sp.ix[:,['per','m','um','r','ur','flux','firstref','orbref']]
ttvs.ix[:,['m','r','rho','urho']]



# <codecell>

from scipy.stats import pearsonr
print "Pearson R, Null prob. on M-R correlation:", pearsonr(sp.r, sp.m)

# <codecell>

def plot_all_exoplanets_mrf():
    '''Plot whole transiting exoplanet population in Radius vs. Mass, Density vs. Mass'''
    fig = plt.figure()
    ax1 = fig.add_subplot(1,2,1, xscale="log", yscale="log") 
    ax2 = fig.add_subplot(1,2,2, xscale="log", yscale="log") 
    ax1.errorbar(df_hot.m, df_hot.r, xerr=df_hot.um, yerr=df_hot.ur, 
                 fmt='o', markersize=5, 
                 alpha=0.5, color='r', label='$F > $ %i $F_E$'%medf, 
                 capsize=0);
    ax1.errorbar(df_cold.m, df_cold.r, xerr=df_cold.um, yerr=df_cold.ur, 
             fmt='o', markersize=5, 
             alpha =0.5, color='b', label='$F < $ %i $F_E$'%medf, 
             capsize=0);
    for x in zip(ss.m,ss.r,ss.index):
        ax1.text(x[0], x[1], x[2], fontsize=12, family='monospace', color='k')
    ax2.errorbar(df_hot.m, df_hot.rho, xerr=df_hot.um, yerr=df_hot.urho, 
             fmt='o', markersize=5, 
             alpha=0.5, color='r', label='$F > $ %i $F_E$'%medf, 
             capsize=0);
    ax2.errorbar(df_cold.m, df_cold.rho, xerr=df_cold.um, yerr=df_cold.urho, 
             fmt='o', markersize=5, 
             alpha=0.5, color='b', label='$F < $ %i $F_E$'%medf,
             capsize=0);
    for x in zip(ss.m,ss.rho,ss.index):
        ax2.text(x[0], x[1], x[2], fontsize=12, family='monospace', color='k')
    ax1.set_xlabel('Planet Mass ($\mathrm{M_E}$)', fontsize=16)
    ax2.set_xlabel('Planet Mass ($\mathrm{M_E}$)', fontsize=16)
    ax1.set_ylabel('Planet Radius ($\mathrm{R_E}$)',fontsize=16)
    ax2.set_ylabel('Planet Density (g/cc)',fontsize=16)
    ax1.legend(loc='lower right');
    ax2.legend(loc='upper right');
    fig.set_size_inches(12,6)
    plt.savefig('mrf.eps')


def plot_histograms():
    '''Plot histograms of exoplanet mass, radius, density'''
    fig, axes = plt.subplots(nrows=1, ncols=3);
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    sp['r'].hist(ax=axes[0], bins = np.arange(0,5,0.5));
    axes[0].set_xlabel('Radius ($\mathrm{R_E}$)', fontsize=16); 
    axes[0].set_ylabel('Number', fontsize=16);
    majorLocator   = MultipleLocator(0.5);
    axes[0].xaxis.set_major_locator(majorLocator);
    sp['m'].hist(ax=axes[1], bins = np.arange(-10,20,2)); 
    axes[1].set_xlabel('Mass ($\mathrm{M_E}$)', fontsize=16);
    majorLocator   = MultipleLocator(4);
    axes[1].xaxis.set_major_locator(majorLocator);
    sp['rho'].hist(ax=axes[2], range=(-10,30), bins = np.arange(-10,30,2)); 
    axes[2].set_xlabel('Density (g/cc)', fontsize=16);
    majorLocator   = MultipleLocator(4);
    axes[2].xaxis.set_major_locator(majorLocator);
    fig.set_size_inches(12,6);
    subplots_adjust(wspace=0.1);
    subplots_adjust(hspace=0.5);
    plt.savefig('histograms.png',dpi=100);
    plt.savefig('histograms.eps', format='eps', dpi=100);
    plt.show()

# numerical weighted linear fit
fitfunc = lambda p, x: p[0] + p[1] * x
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err  
from scipy.optimize import leastsq
out = leastsq(errfunc,  [0.2,2.5], args=(sp.r[sp.r >= 1.5], sp.m[sp.r >= 1.5], 
                                         sp.um[sp.r >= 1.5]),full_output=1)
best_pars = out[0]
print "Linear Mass-Radius Fit for R >= 1.5 R_E:"
print "Best parameters:", best_pars
print "Covariance:", out[1]
print "Chi^2:", (errfunc(out[0],sp_vol.r,sp_vol.m,sp_vol.um)**2.).sum()/(len(sp_vol) - len(out[0]))
print "RMS:", np.sqrt(np.mean((sp_vol.m - fitfunc(best_pars,sp_vol.r))**2.))
#calculate standard errors
ssxx = sum((sp_vol.r - np.mean(sp_vol.r))**2.)
ssyy = sum((sp_vol.m - np.mean(sp_vol.m))**2.)
n = len(sp_vol)
s = np.sqrt(sum((sp_vol.m - fitfunc(best_pars,sp_vol.r))**2./(n-2)))
sea = s*np.sqrt(1./n + np.mean(sp.r[sp.r >= 1.5])**2./ssxx)
seb = s/np.sqrt(ssxx)
print "Standard errors:"
print "Intercept:", sea, "Slope:", seb
# generate fine arrays
rarr = linspace(1.5, 4.0, 101)
marr = fitfunc(best_pars,rarr)
darr = (marr*mearth)/(4./3.*pi*(rarr*rearth)**3)
darr[darr < 0] = 100

# <codecell>

# this is an *analytical* weighted polynomial fit to density vs. radius -- the best!
out4 = np.polynomial.polynomial.polyfit(sp_vol.r, sp_vol.rho, 2, 
                                        rcond=None, full=True, w=1./sp_vol.urho**2.)
best_q = out4[0]
poly2 = lambda p, x: p[0] + p[1] * x + p[2] * x**2.
poly3 = lambda p, x: p[0] + p[1] * x + p[2] * x**2. + p[3] * x**3.
poly_d = best_q[0] + best_q[1] * rarr + best_q[2] * rarr**2.# + best_q[3] * rarr**3.
rho_pred = poly2(best_q,sp_vol.r)
print "Best Analytical Weighted Polynomial Fit to Density vs. Radius:"
print "Best parameters:", out4[0]
print "Covariance:", out4[1]
print "Chi^2:", (((sp_vol.rho - rho_pred)/sp_vol.urho)**2.).sum()/(len(sp_vol.rho) - len(out4[0]))
print "RMS:", np.sqrt(np.mean((sp_vol.rho - rho_pred)**2.))
#print "RMS:", np.sqrt(np.mean((sp[sp.r > 1.5].rho - poly2(best_q,sp[sp.r > 1.5].r))**2.))

# <codecell>

powerlaw = lambda p, x: p[0] * x**(p[1])
powerlaw2 = lambda p, x: pow(x,p[0])
errfunc_plaw = lambda p, x, y, err : (y - powerlaw(p, x))/err
errfunc_plaw2 = lambda p, x, y, err : (y - powerlaw2(p, x))/err

mass_plaw = leastsq(errfunc_plaw,  [3.,1.], 
              args=(sp_vol.r, sp_vol.m, sp_vol.um),full_output=1)
mpars = mass_plaw[0]
mplaw = mpars[0] * pow(rarr,mpars[1])
print "Best Power Law Mass-Radius Relation for R >=1.5 R_E:"
print "Best parameters:", mpars
print "Covariance:", mass_plaw[1]
print "Chi^2:", (errfunc_plaw(mpars,sp_vol.r,sp_vol.m,sp_vol.um)**2.).sum()/(len(sp_vol.m) - len(mpars))
print "RMS:", np.sqrt(np.mean((sp_vol.m - powerlaw(mpars,sp_vol.r))**2))

density_plaw = leastsq(errfunc_plaw,  [8.,-1.], 
                       args=(sp_vol.r, sp_vol.rho, sp_vol.urho),full_output=1)
dpars = density_plaw[0]
dplaw = dpars[0] * pow(rarr,dpars[1])
print "********************************"
print "Best Power Law Density-Radius Relation for R >=1.5 R_E:"
print "Best parameters:", dpars
print "Covariance:", density_plaw[1]
print "Chi^2:", (errfunc_plaw(dpars,sp_vol.r,sp_vol.rho,sp_vol.urho)**2.).sum()/(len(sp_vol.m) - len(dpars))
print "RMS:", np.sqrt(np.mean((sp_vol.rho - powerlaw(dpars,sp_vol.r))**2))
#print "RMS:", np.sqrt(np.mean((sp[sp.r > 1.5].rho - powerlaw(dpars,sp[sp.r > 1.5].r))**2))

# <codecell>

# Fits to the rocky planets
rocky_plaw_no_ss = leastsq(errfunc_plaw,  [1.,3.], 
              args=(sp[sp.r < 1.5].r, sp[sp.r < 1.5].m, sp[sp.r < 1.5].um),full_output=1)
rocky_plaw = leastsq(errfunc_plaw,  [1.,3.], 
              args=(sp_rocky.r, sp_rocky.m, sp_rocky.um),full_output=1)
rocky_pars = rocky_plaw[0]
rocky_pars_no_ss = rocky_plaw_no_ss[0]
rocky_r = np.arange(0.01,1.51,0.01)
rockplaw = powerlaw(rocky_pars,rocky_r)
print "********************************"
print "Rocky Mass Power Law Solution:"
print "Best parameters without SS:", rocky_pars_no_ss
print "Best parameters including SS:", rocky_pars
print "Covariance:", rocky_plaw[1]
print "Chi^2 without SS:", (errfunc_plaw(rocky_pars_no_ss,
                              sp[sp.r < 1.5].r,sp[sp.r < 1.5].m,
                              sp[sp.r < 1.5].um)**2.).sum()/(len(sp[sp.r < 1.5].m) - len(mpars))
print "Chi^2 with SS:", (errfunc_plaw(rocky_pars,
                              sp_rocky.r,sp_rocky.m,
                              sp_rocky.um)**2.).sum()/(len(sp_rocky.m) - len(rocky_pars))
print "Chi^2 to mr plaw without SS:", (errfunc_plaw(mpars,
                              sp[sp.r < 1.5].r,sp[sp.r < 1.5].m,
                              sp[sp.r < 1.5].um)**2.).sum()/(len(sp[sp.r < 1.5].m) - len(mpars))
print "Chi^2 to mr plaw with SS:", (errfunc_plaw(mpars,
                              sp_rocky.r,sp_rocky.m,
                              sp_rocky.um)**2.).sum()/(len(sp_rocky.m) - len(mpars))
print "RMS without SS:", np.sqrt(np.mean((sp[sp.r < 1.5].m - powerlaw(mpars,sp[sp.r < 1.5].r))**2))
print "RMS with SS:", np.sqrt(np.mean((sp_rocky.m - powerlaw(rocky_pars,sp_rocky.r))**2))

rocky_dens = leastsq(errfunc_plaw,  [10.,-2.], 
              args=(sp_rocky.r, sp_rocky.rho, sp_rocky.urho),full_output=1)
rocky_d = rocky_dens[0]
rockdlaw = rocky_d[0] * pow(rocky_r,rocky_d[1])
print "********************************"
print "Rocky Density Power Law Solution:"
print "Best parameters:", rocky_d
print "Covariance:", rocky_dens[1]
print "Chi^2:", (errfunc_plaw(rocky_d,
                              sp_rocky.r,sp_rocky.rho,
                              sp_rocky.urho)**2.).sum()/(len(sp_rocky.rho) - len(rocky_d))
print "RMS:", np.sqrt(np.mean((sp_rocky[sp_rocky.rho < 100].rho - powerlaw(rocky_d,sp_rocky[sp_rocky.rho < 100].r))**2))

polytrope = lambda n, x: 4/3.*pi*(x/3.90)**3.*10.55*(1. + (1. - 3./5.*n)*pow((2./3.*pi*(x/3.90)**2.),n))
#polytrope = lambda n, x: 1. + (1. - 3./5.*n)*pow((2./3.*pi*x**2.),n)
#polytrope = lambda k, x: pow(10,k[0]) * pow(x,1/3.)*k[2]*pow(10,k[1]*x)
errfunc_ptrope = lambda n, x, y, err : (y - polytrope(n, x))/err
m_polytrope = leastsq(errfunc_ptrope,  [9/3.], 
              args=(sp_rocky.r, sp_rocky.m, sp_rocky.um),full_output=1)
ptrope = m_polytrope[0]
print "********************************"
print "Rocky Mass Polytrope Solution:"
print "Best parameter:", ptrope
print "Covariance:", m_polytrope[1]
print "Chi^2:", (errfunc_ptrope(ptrope,
                              sp_rocky.r,sp_rocky.m,
                              sp_rocky.um)**2.).sum()/(len(sp_rocky.m) - len(ptrope))
m_ptrope = polytrope(ptrope, rocky_r)

# <codecell>

# find uncertainties in plaw with bootstrap
def draw_data(x,ux):
    x_pbs = np.array([x[i] + random.normal(ux[i]) for i in np.arange(len(x))])
    return x_pbs

def mr_bootstrap(ntrial=100):
    # do 1000 draws and store the power law coefficients
    mr_plaw_pbs = []
    dr_plaw_pbs = []
    dr_poly_pbs = []
    mr_rocky_pbs = []
    dr_rocky_pbs = []
    #ums = [sp.um[i] for i in np.arange(len(sp))]
    #uds = [sp.urho[i] for i in np.arange(len(sp))]
    for i in xrange(0,ntrial):
    #if "x"=="x":
        # draw a mass and density
        ms = draw_data(sp.m, sp.um)
        ds = draw_data(sp.rho, sp.urho)
        ms_rocky = draw_data(sp_rocky.m, sp_rocky.um)
        ds_rocky = draw_data(sp_rocky.rho, sp_rocky.urho)
    #    plt.errorbar(rs,ms,yerr=ums,fmt='o')
        # calculate best fits for these m's and d's
        mr_plaw_pbs.append(leastsq(errfunc_plaw,  [3.,1.], args=(sp.r,ms,sp.um))[0])
        dr_plaw_pbs.append(leastsq(errfunc_plaw,  [10.,-2.], args=(sp.r,ds,sp.urho))[0])
        dr_poly_pbs.append(np.polynomial.polynomial.polyfit(sp.r, ds, 2, rcond=None, w=1./sp.urho**2.))
        mr_rocky_pbs.append(leastsq(errfunc_plaw,  [1.,3.], args=(sp_rocky.r, ms_rocky, sp_rocky.um))[0])
        dr_rocky_pbs.append(leastsq(errfunc_plaw,  [10.,-2.], args=(sp_rocky.r, ds_rocky, sp_rocky.urho))[0])
    print "mr power law coefficients:"
    print "mean:", np.mean(mr_plaw_pbs,axis=0), "\t","\t","std:",np.std(mr_plaw_pbs,axis=0)
    print "****************************************************************"
    #print "dr power law coefficients:"
    #print "mean:", np.mean(dr_plaw_pbs,axis=0),"\t", "\t","std:",np.std(dr_plaw_pbs,axis=0) 
    #print "****************************************************************"
    print "dr polynomial coefficients:"
    print "mean:", np.mean(dr_poly_pbs, axis=0),"\t", "std:",np.std(dr_poly_pbs, axis=0)
    print "****************************************************************"
    print "mr rocky power law coefficients:"
    print "mean:", np.mean(mr_rocky_pbs, axis=0),"\t", "\t","std:", np.std(mr_rocky_pbs,axis=0)
    print "****************************************************************"
    print "dr rocky power law coefficients:"
    print "mean:", np.mean(dr_rocky_pbs,axis=0),"\t","\t", "std:", np.std(dr_rocky_pbs,axis=0)
# <codecell>

#mr_bootstrap

# <codecell>



# <codecell>

# calculated predict masses, densities
#rpred = (sp.m - best_pars[0])/best_pars[1]
vol_mpred = powerlaw(mpars,sp_vol.r)
#dpred = mpred*mearth/(4./3.*(rpred*rearth)**3.)
#r_resid = sp.r - rpred
vol_m_resid = sp_vol.m - vol_mpred
#d_resid = sp.rho - dpred
sp_vol['m_resid'] = vol_m_resid
rocky_mpred = powerlaw(rocky_pars,sp_rocky.r)
rocky_m_resid = sp_rocky.m - rocky_mpred
sp_rocky['m_resid'] = rocky_m_resid

# <codecell>

# bin mass, radius, density
rwm, rwd, urw, mw, umw, rhow, urhow = [], [], [], [], [], [], []
for r in np.arange(0,4.,0.5):
    rbin = sp.loc[r <= sp.r].loc[sp.r < r+0.5]
    rwm.append(np.average(rbin['r'], weights = (1./rbin.um**2.)))
    rwd.append(np.average(rbin['r'], weights = (1./rbin.urho**2.)))
    mw.append(np.average(rbin['m'], weights = 1./rbin.um**2.))
    rhow.append(np.average(rbin['rho'], weights = 1./rbin.urho**2.))
    urw.append(std(rbin.r)/sqrt(len(rbin)))
    umw.append(std(rbin.m)/sqrt(len(rbin)))
    urhow.append(std(rbin.rho)/sqrt(len(rbin)))
means = pd.DataFrame({'r_m':rwm, 'r_d':rwd, 'ur':urw, 'm':mw, 'um':umw, 'rho':rhow, 'urho':urhow})

# <codecell>

# Calculate Earth composition from Seager 2007
# this is good for Ms < 4
m1 = 6.41; r1 = 3.19
k1 = -0.20945; k2 = 0.0804; k3 = 0.394
m_phys = np.logspace(-3, 1, num=50)
ms = m_phys/m1
logrs = k1 + 1/3. * np.log10(ms) - k2 * ms**k3
rs = pow(10.,logrs)
print max(ms), max(rs*r1)
m_seager = ms*m1; r_seager = rs*r1
d_seager = m_seager/(4/3.*pi*r_seager**3.)*mearth/rearth**3.
#print zip(r_seager,d_seager)

# <codecell>

# plot mass, density vs. radius
def plot_mr():
    fig = plt.figure();
    sharex=True;
    ax1 = fig.add_subplot(1,2,1);
    ax1.errorbar(sp.r, sp.m, xerr=sp.ur, yerr=sp.um, fmt='o', markersize=5, 
             alpha =0.3, color='gray', capsize=0, elinewidth=2); # black option
    ax1.plot(rarr, mplaw, '--k', linewidth=4, label=("M = %.2f * R^ %.2f"% (mpars[0], mpars[1])));
    ax1.plot(rarr, marr, '--r', linewidth=4, label=("M = %.2f + %.2f R" % (best_pars[0], best_pars[1])));
    ax1.plot(rocky_r, rockplaw, '--c', linewidth=4, label=("M = %.2f R^ %.2f"% (rocky_pars[0], rocky_pars[1])));
    #ax1.plot(rocky_r, m_ptrope, '--g', linewidth=4, label=("polytrope n= %.1f"% ptrope));
    ax1.errorbar(means.r_m,means.m, xerr=means.ur, yerr=means.um, 
             markersize=10, marker='s',fmt='o', capsize=0, elinewidth=4, alpha = 0.8);
    ax1.errorbar(ttvs.r, ttvs.m, xerr=ttvs.ur, yerr=ttvs.um, fmt='o',capsize=0, elinewidth=2, alpha=0.5, color='orange', label=('TTVs'));
    ax1.set_xlabel('Planet Radius ($\mathrm{R_E}$)',fontsize=24);
    ax1.set_ylabel('Planet Mass ($\mathrm{M_E})$',fontsize=24);
    for x in zip(ss.r,ss.m,ss.index):
        if (x[1] < 20) and (x[1] > -3) and (x[0] < 4.):
            ax1.text(x[0], x[1], x[2], fontsize=12, family='monospace', color='m');
    ax1.legend(loc='upper left');
    plt.axis([0.0,4.1,-3,20]);

    ax2 = fig.add_subplot(1,2,2);
    ax2.plot(r_seager,d_seager, color='g', alpha = 0.5, linewidth=4, label=('Earth Composition (Seager 2007)'));
    ax2.errorbar(sp.r, sp.rho, xerr=sp.ur, yerr=sp.urho, fmt='o', markersize=5, 
             alpha=0.3, color='gray',capsize=0, elinewidth=2);
    #ax2.plot(rarr, poly_d, '--k', linewidth=4, label=("Rho = %.2f + %.2f R + %.2f R^2" % (best_q[0], best_q[1], best_q[2])));
    ax2.plot(rarr, dplaw, '--k', linewidth=4, label=("Rho = %.2f * R^ %.2f"% (dpars[0], dpars[1])));
    ax2.plot(rocky_r, rockdlaw, '--c', linewidth=4, label=("Rho = %.2f * R^ %.2f"% (rocky_d[0], rocky_d[1])));
    ax2.set_xlabel('Planet Radius ($\mathrm{R_E}$)',fontsize=24);
    ax2.set_ylabel('Planet Density ($\mathrm{g/cm^{3}}$)',fontsize=24);
    ax2.errorbar(means.r_d,means.rho, xerr=means.ur, yerr=means.urho, 
             markersize=10, marker='s', fmt='o', capsize=0, elinewidth=4, alpha=0.8);
    ax2.errorbar(ttvs.r, ttvs.rho, xerr=ttvs.ur, yerr=ttvs.urho, fmt='o',capsize=0, elinewidth=2, alpha=0.5, color='orange',label=('TTVs'));
    #ax2.plot(xnew, ynew, '--', c='red',linewidth=4, label=('Univariate Spline (k=3)'));
    #ax2.errorbar(sp.ix['Kepler-78 b','r'], sp.ix['Kepler-78 b','rho'], 
    #             xerr=sp.ix['Kepler-78 b','ur'], yerr=sp.ix['Kepler-78 b','urho'], fmt='o',capsize=0, elinewidth=4, 
    #             alpha=0.5, color='r');
    for x in zip(ss.r,ss.rho,ss.index):
        if (x[1] < 30) and (x[1] > -1) and (x[0] < 4.):
            ax2.text(x[0], x[1], x[2], fontsize=12, family='monospace', color='m');
    plt.axis([0.,4.0,-5,40]);
    ax2.legend();
    fig.set_size_inches(12,6);
    subplots_adjust(wspace=0.3);
    subplots_adjust(hspace=0.2);
    plt.savefig('mr_small.png', dpi=100)
    plt.savefig('mr_small.eps', format='eps', dpi=100)
    plt.show()

# <codecell>

def plot_flux_density():
    fig = plt.figure(figsize=(8,8));
    ax = fig.add_subplot(111,xscale='log');
    rho_resid = sp.rho - rho_pred
    #notnan = numpy.logical_not(numpy.isnan(sp.flux))
    ax.errorbar(sp_clean.flux, sp_clean.rho, 
            xerr=sp_clean.uflux, yerr=sp_clean.urho, 
            fmt='o', alpha = 0.5, capsize=0);
    ax.set_xlabel('Flux (Earth flux)', fontsize=18);
    ax.set_ylabel('Planet Density (g/cc)', fontsize=18);
    plt.axis([0,4000,-5,15])
    fig.show();

#plot residuals
sp_conglom = sp_conglom = sp_rocky.append(sp_vol)
def plot_residuals(sp=sp_conglom):
    '''Plot M-R residuals vs. other physical parameters'''
    fig = plt.figure(figsize=(14,14));
    xarr = [sp.per, sp.a, sp.flux, sp.mstar, sp.rstar, sp.logg, sp.fe, sp.age, sp.vsini]
    xarr_err = [None, sp.ua, sp.uflux, sp.umstar, sp.urstar, 
            None, sp.ufe, None, None]
    xtit = ['Period (d)','a (AU)','Flux (Earth)',
        '$M_\star$ ($\mathrm{M_{Sun}}$)','$R_\star$ ($\mathrm{R_{Sun}}$)','log(g)',
        '[Fe/H]','Age (Gyr)','Vsini (km/s)']
    xscale = ['log','log','log','linear','linear','linear','linear','linear','log']
    fig, axarr = plt.subplots(3, 3, sharey=True);
    ys = sp.m_resid
    yerr = sp.um
    yttv = ys[sp.ttv==1]
    yerrttv = yerr[sp.ttv==1]
    for i in xrange(0,3):
        for j in xrange(0,3):
            xs = xarr[3*i + j]
            xttv = xs[sp.ttv==1]
            xerr = xarr_err[3*i + j]
            axarr[i,j].errorbar(xs,ys, xerr=xerr, yerr=yerr, alpha = 0.5, marker='o', fmt='o')
            if len(xttv) > 0:
                if xerr is not None: xerrttv = xerr[sp.ttv==1] 
                else: xerrttv=None
            axarr[i,j].errorbar(xttv,yttv,xerr=xerrttv,yerr=yerrttv, 
                            alpha = 0.5, marker='o', fmt='o', c='orange')
            if xscale[3*i + j]=='log': axarr[i,j].set_xscale('log')
            if j==0: axarr[i,j].set_ylabel('Residual Mass ($\mathrm{M_E}$)', fontsize=16)
            axarr[i,j].set_xlabel(xtit[3*i+j], fontsize=16)
    fig.set_size_inches(12,12);
    subplots_adjust(wspace=0.1)
    subplots_adjust(hspace=0.5)
    #fig.ylabel('Residual Mass')
    fig.show();
    plt.savefig('mr_resids.png',dpi=100)
    plt.savefig('mr_resids.eps', format='eps', dpi=100)

# <codecell>

def test_ttv_population():
    from scipy.stats import ttest_ind
    ttest_ind(sp.m_resid.loc[sp.ttv==0].loc[sp.r > 1.0], sp.m_resid.loc[sp.ttv==1].loc[sp.r > 1.0])

def test_fe_correlation():
    notnan = numpy.logical_not(numpy.isnan(sp.fe))
    from scipy.stats.stats import pearsonr
    pearsonr(sp[notnan].fe, m_resid[notnan])
    print "r = ",fe_r[0], "p =",1-fe_r[1]
    print "Prob. Fe and resid correlation is real:",(1.-fe_r[1])**9.
plot_histograms()
plot_mr()
plot_residuals()