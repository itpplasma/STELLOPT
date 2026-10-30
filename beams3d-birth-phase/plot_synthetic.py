from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path(__file__).resolve().parents[1]
def load(n):return json.loads((r/(n+'.json')).read_text())
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
b,a=load('beams-before'),load('beams-after');x=np.arange(1,4)
fig,axes=plt.subplots(1,2,figsize=(9,3.8),layout='constrained')
for ax,key,label in zip(axes,['position_error_mm','velocity_error_km_s'],['Birth position error [mm]','Incident velocity error [km/s]']):
 ax.axhline(0,color='black',lw=1,label='Physical birth reference')
 ax.plot(x,b[key],'o--',color='#D55E00',label='Before')
 ax.plot(x,a[key],'x-',color='#0072B2',ms=8,label='After')
 ax.set(xlabel='Synthetic marker ID',ylabel=label,xticks=x);ax.grid(alpha=.2)
axes[0].legend(fontsize=9);fig.savefig(r/'birth-phase-synthetic.png',dpi=180);fig.savefig(r/'birth-phase-synthetic.pdf');plt.close(fig)
b,a=load('mass-before'),load('mass-after')
fig,ax=plt.subplots(figsize=(6,3.8),layout='constrained')
ax.axhline(1,color='black',lw=1,label=r'Relativistic reference $p=\gamma mv$')
ax.plot(x,b['momentum_ratio'],'o--',color='#D55E00',label='Before (serialized amu)')
ax.plot(x,a['momentum_ratio'],'x-',color='#0072B2',ms=8,label='After')
ax.set(xlabel='Synthetic marker ID',ylabel='Momentum / relativistic reference [1]',yscale='log',xticks=x);ax.grid(alpha=.2);ax.legend(fontsize=9)
fig.savefig(r/'particle-mass-synthetic.png',dpi=180);fig.savefig(r/'particle-mass-synthetic.pdf')
