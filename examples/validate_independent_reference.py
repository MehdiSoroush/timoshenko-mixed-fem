"""Independent conventional [w, phi] FEM and continuum frequency validation.

Independent numerical verification with a displacement–rotation finite element reference.
Run: python examples/validate_independent_reference.py
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from timoshenko_mixed import Beam,solve_static,midpoint_deflection,solve_dynamic
from timoshenko_mixed.reference import (solve_reference_static,solve_reference_dynamic,
                                         reference_first_frequencies)
from timoshenko_mixed.analytical import simply_supported_timoshenko_frequency

ROOT=Path(__file__).resolve().parents[1]/'results'; (ROOT/'figures').mkdir(parents=True,exist_ok=True);(ROOT/'data').mkdir(parents=True,exist_ok=True)
beam=Beam(); rows=[]
for n in (9,17,33,65,129):
    m=solve_static(n,beam)
    r=solve_reference_static(n,beam=beam)
    mm=midpoint_deflection(m);rr=float(np.interp(beam.length/2,r.nodes,r.deflection))
    rows.append((n,mm,rr,(mm-rr)/rr))
print('Mesh / midspan mixed / independent reference / relative difference')
for v in rows: print('%3d  %.9f m  %.9f m  %.6f%%'%(v[0],v[1],v[2],v[3]*100))
with open(ROOT/'data'/'independent_static_comparison.csv','w',newline='') as f:
    w=csv.writer(f);w.writerow(['n_nodes','mixed_midspan_m','conventional_fem_midspan_m','signed_relative_difference']);w.writerows(rows)

n=65;m=solve_static(n,beam);r=solve_reference_static(n,beam=beam)
analytic=simply_supported_timoshenko_frequency(beam)
freqs=[]
for count in (9,17,33,65):
    freq=reference_first_frequencies(beam,count,count=1)[0]
    freqs.append((count,freq,analytic,(freq-analytic)/analytic))
with open(ROOT/'data'/'independent_frequency_comparison.csv','w',newline='') as f:
    w=csv.writer(f);w.writerow(['n_nodes','conventional_fem_hz','analytical_timoshenko_hz','relative_error']);w.writerows(freqs)

n_dynamic=41; timestep=.001;finish=.2
base=dict(n_nodes=n_dynamic,beam=beam,boundary='fixed-pinned',load_case='constant',dt=timestep,t_end=finish)
mixed=solve_dynamic(**base,integrator='newmark-average')
ref=solve_reference_dynamic(**base)
wm=mixed.midpoint_deflection();wr=ref.midpoint_deflection()
rms=float(np.sqrt(np.mean((wm-wr)**2)))
print('Independent dynamic reference RMS difference: %.6f mm (%.3f%% of reference peak)'%(rms*1000,100*rms/np.max(abs(wr))))
with open(ROOT/'data'/'independent_dynamic_comparison.csv','w',newline='') as f:
    w=csv.writer(f);w.writerow(['time_s','mixed_newmark_midspan_m','physical_reference_newmark_midspan_m'])
    w.writerows(zip(mixed.times,wm,wr))
fig,axs=plt.subplots(3,1,figsize=(8.6,10.5))
axs[0].plot(m.nodes,m.deflection*1000,label='Mixed [w, M] FEM',lw=2)
axs[0].plot(r.nodes,r.deflection*1000,'--',label='Independent [w, rotation] FEM',lw=1.8)
axs[0].set(ylabel='Deflection (mm)',title='Independent static verification: fixed–pinned beam')
axs[0].legend()
axs[1].semilogy([x[0] for x in rows],[abs(x[3])*100 for x in rows],'o-',label='Mixed/reference midspan difference')
axs[1].set(ylabel='Relative difference (%)',xlabel='Number of beam nodes',title='Both discretizations converge to the same static solution')
axs[1].legend()
axs[2].plot(mixed.times,wm*1000,label='Mixed Newmark',lw=1.7)
axs[2].plot(ref.times,wr*1000,'--',label='Independent physical FE Newmark',lw=1.4)
axs[2].set(xlabel='Time (s)',ylabel='Midspan deflection (mm)',title='Independent dynamic check (41 nodes, 1 ms timestep)')
axs[2].legend()
for ax in axs:ax.grid(alpha=.3)
fig.tight_layout();fig.savefig(ROOT/'figures'/'independent_fem_validation.png',dpi=170);plt.close(fig)
print('Simply-supported continuum Timoshenko f1 = %.6f Hz'%analytic)

fig,ax=plt.subplots(figsize=(7.6,4.5))
ax.plot([n for n,*_ in freqs],[f for _,f,*_ in freqs],'o-',label='Independent physical [w, rotation] FEM')
ax.axhline(analytic,color='black',linestyle='--',label='Exact Timoshenko continuum dispersion')
ax.set(xlabel='Number of nodes',ylabel='First bending frequency (Hz)',title='Independent continuum natural-frequency verification')
ax.grid(alpha=.3);ax.legend();fig.tight_layout()
fig.savefig(ROOT/'figures'/'independent_frequency_verification.png',dpi=180);plt.close(fig)
