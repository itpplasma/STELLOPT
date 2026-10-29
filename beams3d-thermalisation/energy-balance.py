# Net electron-plus-ion deposition; data are synthetic BEAMS3D benchmark runs.
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root = Path(__file__).parent
rows = json.loads((root/'energy-balance.json').read_text())
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,2,figsize=(8.5,4.8))
fig.subplots_adjust(left=.09, right=.98, top=.82, bottom=.25, wspace=.3)
for ax,model in zip(axes,['Guiding centre','Full orbit']):
 rr=[r for r in rows if r['model']==model];ref=np.mean([r['birth_keV'] for r in rr])
 ax.axhspan(.97*ref,1.03*ref,color='.9',label='Acceptance: ±3%')
 ax.axhline(ref,color='black',linestyle='--',linewidth=1.4,label='Birth energy (no losses)')
 total=[r['electron_keV'] + r['ion_keV'] for r in rr]
 ax.bar(0,total[0],.55,color='#E69F00',hatch='//')
 ax.bar(1,total[1],.55,color='#0072B2')
 for x,r in enumerate(rr):
  ax.text(x,r['electron_keV']+r['ion_keV']+.025*ref,f"{100*r['balance']:.1f}%",ha='center',fontsize=11)
 ax.set(xticks=[0,1],xticklabels=['Before','After'],ylim=(0,1.22*ref),ylabel='Deposited energy per marker [keV]',title=model)
 ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
fig.legend(*axes[0].get_legend_handles_labels(), loc='lower center',
           bbox_to_anchor=(.5,.025), ncol=2, fontsize=9)
fig.suptitle('Energy deposited when markers join the bulk ions')
fig.savefig(root/'energy-balance.png',dpi=180)
