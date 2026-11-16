"""Generate auditable time-history plots and CSV from the 2017 beam problem.

Run from repository root: python examples/reproduce_dynamic_cases.py
Computed figures and datasets are produced by the Python implementation.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from timoshenko_mixed import solve_dynamic

OUT = Path('results')
(OUT / 'figures').mkdir(parents=True, exist_ok=True)
(OUT / 'data').mkdir(parents=True, exist_ok=True)

cases = ['constant', 'sine', 'source-ramp']
sols = {name: solve_dynamic(21, t_end=2.0, dt=.002, load_case=name) for name in cases}

fig, ax = plt.subplots(figsize=(9.0, 4.8))
for name, sol in sols.items():
    ax.plot(sol.times, 1000*sol.midpoint_deflection(),label=name.replace('-', ' ').title())
ax.axvline(1,ls='--',lw=1,label='Distributed load ends (t = 1 s)')
ax.set(xlabel='Time (s)', ylabel='Midspan deflection (mm)',
       title='2017 source-derived mixed beam model | 3 distributed load profiles')
ax.grid(alpha=.3)
ax.legend()
fig.tight_layout()
fig.savefig(OUT / 'figures' / 'dynamic_three_loading_cases.png',dpi=180)
plt.close(fig)

sample = sols['constant']
newmark = solve_dynamic(21,t_end=2,dt=.002,load_case='constant',integrator='newmark-average')
fig, ax = plt.subplots(figsize=(9.0, 4.8))
ax.plot(sample.times,1000*sample.midpoint_deflection(),label='2017 weighted-residual 3-level method')
ax.plot(newmark.times,1000*newmark.midpoint_deflection(),label='Newmark average acceleration (new comparison)',alpha=.85)
ax.axvline(1,ls='--',lw=1)
ax.set(xlabel='Time (s)',ylabel='Midspan deflection (mm)',
       title='Same semidiscrete mixed model | different time integration')
ax.grid(alpha=.3);ax.legend()
fig.tight_layout()
fig.savefig(OUT/'figures'/'time_integrator_comparison.png',dpi=180)
plt.close(fig)

with (OUT/'data'/'dynamic_midspan_time_histories.csv').open('w',newline='') as f:
    writer=csv.writer(f)
    writer.writerow(['time_s','constant_mm','sine_mm','source_ramp_mm','constant_newmark_mm'])
    for i,t in enumerate(sample.times):
        writer.writerow([f'{t:.6g}',
            *[f'{1000*sols[c].midpoint_deflection()[i]:.9g}' for c in cases],
            f'{1000*newmark.midpoint_deflection()[i]:.9g}'])

summary=[]
for name, sol in list(sols.items())+[('constant-newmark',newmark)]:
    mid=sol.midpoint_deflection()
    summary.append([name,sol.integrator,sol.deflection.shape[1],sample.times[1],
                   float(mid.max()*1000),float(mid.min()*1000),
                   sol.startup_residual,sol.maximum_step_residual])
with (OUT/'data'/'dynamic_case_summary.csv').open('w',newline='') as f:
    writer=csv.writer(f)
    writer.writerow(['load_case','integrator','n_nodes','dt_s','midspan_max_mm','midspan_min_mm',
                     'startup_equation_residual_inf','step_equation_residual_inf'])
    writer.writerows(summary)

print('Wrote dynamic figures and CSV results under results/.')
for item in summary:
    print(f'{item[0]:18} {item[1]:20} peak={item[4]:.5f} mm')
