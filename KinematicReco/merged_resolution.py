#!/usr/bin/env python3
########################
import sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches

# matplotlib.rcParams['pdf.fonttype'] = 42
# matplotlib.rcParams['ps.fonttype'] = 42

plt.rcParams.update({
    "font.family": "sans-serif",
    # "font.sans-serif": ["Liberation Sans"],
    "font.sans-serif": ["DejaVu Sans"],
})

#Output pdf files
pp1 = PdfPages('plots/merged_resolution_1e-2.pdf')

#Marker scale factor
scaleup = 8

#y cut
ylow_cut = 1e-2
yhi_cut = 0.95

#Read data file
# df = pd.read_excel('input/18_275_resolutions.xlsx',engine='openpyxl')
df = pd.read_csv('./18_275_resolutions.csv')
#print(df)

#Scattered Electron
# df_se = df[ (df['SE']>0) & (df['SE']<df['DA']) & (df['SE']<df['ESIG']) & (df['y']>ylow_cut) & (df['y']<yhi_cut) ] #requires all 3 methods to exist
df_se = df[ (df['SE']>0) & (df['SE']<df['DA']) & (df['SE']<df['SIG']) & (df['y']>ylow_cut) & (df['y']<yhi_cut) ] #requires all 3 methods to exist
#df_se = df[ (df['SE']>0) & (df['SE']<df['JB']) & (df['SE']<df['DA']) & (df['SE']<df['SIG'])  ] #requires all 4 methods to exist
#df_se = df[ (df['SE']>0) & ( (df['SE']<df['JB']) | (df['JB']==0) ) & ( (df['SE']<df['DA']) | (df['DA']==0) ) & ( (df['SE']<df['SIG']) | (df['SIG']==0) ) ]
#print(df_se)

x = df_se['x'].tolist()
Q2 = df_se['Q2'].tolist()
resolution = df_se['SE'].tolist()
resolution = [a * scaleup for a in resolution]

plt.scatter(x,Q2,s=resolution)

#JB
#df_jb = df[ (df['JB']>0) & (df['JB']<df['SE']) & (df['JB']<df['DA']) & (df['JB']<df['SIG']) ] #requires all 4 methods to exist
#df_jb = df[ (df['JB']>0) & ( (df['JB']<df['SE']) | (df['SE']==0) ) & ( (df['JB']<df['DA'])  | (df['DA']==0) ) ]
#print(df_jb)

#x = df_jb['x'].tolist()
#Q2 = df_jb['Q2'].tolist()
#resolution = df_jb['JB'].tolist()
#resolution = [a * scaleup for a in resolution]

#plt.scatter(x,Q2,s=resolution)

#DA
# df_da = df[ (df['DA']>0) & (df['DA']<df['SE']) & (df['DA']<df['ESIG']) & (df['y']>ylow_cut) & (df['y']<yhi_cut) ] #requires all 3 methods to exist
df_da = df[ (df['DA']>0) & (df['DA']<df['SE']) & (df['DA']<df['SIG']) & (df['y']>ylow_cut) & (df['y']<yhi_cut) ] #requires all 3 methods to exist
#df_da = df[ (df['DA']>0) & (df['DA']<df['SE']) & (df['DA']<df['JB']) & (df['DA']<df['SIG']) ] #requires all 4 methods to exist
#df_da = df[ (df['DA']>0) & ( (df['DA']<df['SE']) | (df['SE']==0) ) & ( (df['DA']<df['JB']) | (df['DA']==0) ) ]
#print(df_da)

x = df_da['x'].tolist()
Q2 = df_da['Q2'].tolist()
resolution = df_da['DA'].tolist()
resolution = [a * scaleup for a in resolution]

plt.scatter(x,Q2,s=resolution)

#SIG
df_sig = df[ (df['SIG']>0) & (df['SIG']<df['SE']) & (df['SIG']<df['DA']) & (df['y']>ylow_cut) & (df['y']<yhi_cut) ] #requires all 3 methods to exist
#df_sig = df[ (df['SIG']>0) & (df['SIG']<df['SE']) & (df['SIG']<df['JB']) & (df['SIG']<df['DA']) ] #requires all 4 methods to exist
#print(df_sig)

x = df_sig['x'].tolist()
Q2 = df_sig['Q2'].tolist()
resolution = df_sig['DA'].tolist()
resolution = [a * scaleup for a in resolution]

plt.scatter(x,Q2,s=resolution)

#Make plot
plt.subplots_adjust(left=0.125, bottom=0.125, right=0.875, top=0.875)
plt.xlabel('x', loc='right')
plt.ylabel(r'$\mathrm{Q^{2} \thinspace [GeV^{2}]}$', loc='top')
plt.xscale('log')
plt.yscale('log')
axes = plt.gca()
axes.xaxis.label.set_size(22)
axes.yaxis.label.set_size(22)
axes.tick_params(axis='both',labelsize=14)
axes.spines["left"].set_linewidth(2)
axes.spines["right"].set_linewidth(2)
axes.spines["bottom"].set_linewidth(2)
axes.spines["top"].set_linewidth(2)
plt.xlim([1e-5,1])
plt.ylim([8e-1,1e4])

#Add text
# plt.text(1.5e-5,4e3,'18 GeV e$^{-}$ on 275 GeV p',fontsize=15,color='black')
plt.text(1.5e-5,1e3,'Best Reconstruction Method for y',fontsize=12,color='black')
plt.text(1.5e-5,6e2,'Electron Method',fontsize=12,color='blue')
#plt.text(1.5e-5,3.7e2,'JB Method',fontsize=12,color='orange')
#plt.text(1.5e-5,2.25e2,'DA Method',fontsize=12,color='green')
#plt.text(1.5e-5,1.5e2,r'$\Sigma$ Method',fontsize=12,color='red')

plt.text(1.5e-5,2.25e2,'Double-Angle Method',fontsize=12,color='orange')
plt.text(1.5e-5,3.7e2,r'$\Sigma$ Method',fontsize=12,color='green')

plt.text(1.5e-5,6e1,'Resolution on y',fontsize=12,color='black')
yy = [3e0,6e0,1.25e1,3e1]
xx = [2e-5]*len(yy)
s = [1.*scaleup,5.*scaleup,10.*scaleup,25.*scaleup]
plt.scatter(xx,yy,s=s,color='gray')
plt.text(3e-5,2.7,'1 %',fontsize=10,color='black')
plt.text(3e-5,5.5,'5 %',fontsize=10,color='black')
plt.text(3e-5,11,'10 %',fontsize=10,color='black')
plt.text(3e-5,27,'25 %',fontsize=10,color='black')

x1 = [1.5e-5,9e-4]
y1 = [9.25e2,9.25e2]
plt.plot(x1,y1,color='black')

x2 = [1.5e-5,1e-4]
y2 = [55,55]
plt.plot(x2,y2,color='black')

#Add y limits
s_cm = 4.*18.*275.
#s_cm = 4.*5.*41.
xlim = [1.0e-5,1.0]
Q2_lim1 = [s_cm*1.0e-2*x for x in xlim]
Q2_lim2 = [s_cm*0.95*x for x in xlim]
Q2_lim3 = [s_cm*1.0e-3*x for x in xlim]

plt.plot(xlim,Q2_lim1,color='red',linewidth=2)
plt.plot(xlim,Q2_lim2,color='red',linewidth=2)
# plt.plot(xlim,Q2_lim3,color='red',linewidth=2)
plt.text(0.5,1.2e4,'y = 0.95',fontsize=12,color='red',rotation=26)
plt.text(1.1,200,'y = 0.01',fontsize=12,color='red',rotation=26)
# plt.text(1.1,20,'y = 0.001',fontsize=12,color='red',rotation=26)


plt.text(1.5e-5,7.5e3, "ePIC Performance",
    #transform=plt.transAxes,
    fontsize=18,
    fontweight='bold',
    verticalalignment='top',
    horizontalalignment='left') 

plt.text(1.5e-5,4e3, r"$\mathrm{e{+}p,~\sqrt{s}=140~GeV}$",
    #transform=plt.transAxes,
    fontsize=15,
    verticalalignment='top',
    horizontalalignment='left')

plt.text(
    7e-4, 2e4,
    r'$\mathit{Simu\ campaign:\ 10/2025}$',
    ha="right", va="top",
    fontsize=14
)

#Figure Output
#plt.tight_layout()
fig = plt.gcf()
fig.set_size_inches(18.5/2, 10.5/2)
fig.savefig(pp1, format='pdf')

plt.show()

#Close output pdfs
pp1.close()

# end of script
