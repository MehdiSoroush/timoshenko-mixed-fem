"""Modal time-integration verification against an exact semidiscrete mode.

This checks temporal integration of the FE matrix pencil; it is NOT a
continuum closed-form or Abaqus validation.
"""
from pathlib import Path
import csv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eig
from timoshenko_mixed import Beam, assemble_mass, solve_dynamic
from timoshenko_mixed.model import assemble_stiffness

OUT=Path('results')
(OUT/'figures').mkdir(parents=True,exist_ok=True)
(OUT/'data').mkdir(parents=True,exist_ok=True)
beam=Beam()
n_nodes=9
x,K=assemble_stiffness(beam,n_nodes)
_,M=assemble_mass(beam,n_nodes)
fixed=[0,1,2*n_nodes-2,2*n_nodes-1]
free=np.array([i for i in range(2*n_nodes) if i not in fixed])
lam, phi=eig(K.toarray()[np.ix_(free,free)],M.toarray()[np.ix_(free,free)])
positive=np.where(np.isfinite(lam) & (lam.real>0) & (np.abs(lam.imag)<1e-5*lam.real))[0]
idx=positive[np.argmin(lam[positive].real)]
omega=np.sqrt(lam[idx].real)
initial=np.zeros(2*n_nodes)
initial[free]=phi[:,idx].real
initial*=.001 / np.max(np.abs(initial[::2]))

sol1=solve_dynamic(n_nodes,boundary='pinned-pinned',point_force=0,distributed_intensity=0,
                   initial_dofs=initial,integrator='source-three-level',dt=.001,t_end=.08)
sol2=solve_dynamic(n_nodes,boundary='pinned-pinned',point_force=0,distributed_intensity=0,
                   initial_dofs=initial,integrator='newmark-average',dt=.001,t_end=.08)
exact=initial[2*(n_nodes//2)]*np.cos(omega*sol1.times)
fig,ax=plt.subplots(figsize=(9,4.8))
ax.plot(sol1.times,1000*exact,label=f'Exact semidiscrete mode ({omega/(2*np.pi):.2f} Hz)')
ax.plot(sol1.times,1000*sol1.midpoint_deflection(),label='2017 weighted-residual scheme')
ax.plot(sol2.times,1000*sol2.midpoint_deflection(),label='Newmark average acceleration',ls='--')
ax.grid(alpha=.3);ax.legend()
ax.set(xlabel='Time (s)',ylabel='Midspan deflection (mm)',
       title='Free vibration of the first discretized mode | numerical damping audit')
fig.tight_layout()
fig.savefig(OUT/'figures'/'modal_time_integration_verification.png',dpi=180)
plt.close(fig)

freq=[]
for n in (5,7,9,13,17,25):
 x,K=assemble_stiffness(beam,n)
 _,M=assemble_mass(beam,n)
 fixed=[0,1,2*n-2,2*n-1]
 free=np.array([i for i in range(2*n) if i not in fixed])
 ev=eig(K.toarray()[np.ix_(free,free)], M.toarray()[np.ix_(free,free)],right=False)
 v=ev[np.isfinite(ev)&(ev.real>0)&(np.abs(ev.imag)<1e-5*ev.real)]
 freq.append((n,float(np.sqrt(np.min(v.real))/(2*np.pi))))
ref_freq=(np.pi**2 / beam.length**2) * np.sqrt(beam.EI/(beam.density*beam.area)) / (2*np.pi)
fig,ax=plt.subplots(figsize=(8,4.5))
ax.plot([e[0] for e in freq],[e[1] for e in freq],'-o',label='Mixed FEM, first mode')
ax.axhline(ref_freq,ls='--',label='Euler-Bernoulli beam (slender-beam reference)')
ax.set(xlabel='Number of FE nodes',ylabel='Frequency (Hz)',title='Simply-supported beam | first natural-frequency convergence')
ax.grid(alpha=.3);ax.legend();fig.tight_layout()
fig.savefig(OUT/'figures'/'first_frequency_mesh_convergence.png',dpi=180)
plt.close(fig)
with (OUT/'data'/'natural_frequency_mesh.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['n_nodes','first_fem_frequency_hz','EB_reference_hz'])
 for n,fr in freq:w.writerow([n,f'{fr:.10f}',f'{ref_freq:.10f}'])
print(f'Mixed FE first modal frequency (n=9): {omega/(2*np.pi):.6f} Hz')
print(f'Euler-Bernoulli slender beam estimate: {ref_freq:.6f} Hz')
print('Wrote figures and frequency CSV.')
