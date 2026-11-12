"""Physically labeled bending/shear stress recovery with original-source audit.

Signed stress values follow the mixed nodal moment convention.  The outer
fiber normal stress is M*h/(2*I).  Original source h/(4*I) underestimates
that outer-fiber value by 50%.  Shear has several distinct conventions.
"""
from pathlib import Path
import csv
import matplotlib.pyplot as plt
import numpy as np
from timoshenko_mixed import Beam,solve_static,solve_dynamic
from timoshenko_mixed.stress import recover_element_stress

OUT=Path(__file__).resolve().parents[1]/'results';(OUT/'figures').mkdir(parents=True,exist_ok=True);(OUT/'data').mkdir(parents=True,exist_ok=True)
b=Beam();stat=solve_static(101,beam=b)
e=recover_element_stress(stat.nodes,stat.bending_moment,b)
fig,axs=plt.subplots(2,1,figsize=(8.2,7.2),sharex=True)
axs[0].plot(e.x_mid,e.sigma_outer/1e6,label='Outer fiber M h/(2 I)',lw=2)
axs[0].plot(e.x_mid,e.sigma_source_2017/1e6,'--',label='2017 source M h/(4 I)',lw=1.6)
axs[0].set(ylabel='Signed normal stress (MPa)',title='Static fixed–pinned beam: bending stress (different evaluation locations)')
axs[1].plot(e.x_mid,e.tau_section_average/1e6,label='Section average V/A',lw=2)
axs[1].plot(e.x_mid,e.tau_rectangular_max/1e6,label='Rectangular section peak 1.5 V/A',lw=1.5)
axs[1].plot(e.x_mid,e.tau_effective_timoshenko/1e6,'--',label='2017 effective V/(ks A)',lw=1.5)
axs[1].set(xlabel='Beam position x (m)',ylabel='Signed shear stress (MPa)',title='Shear stress: average, physical peak, and effective Timoshenko')
for ax in axs: ax.grid(alpha=.3);ax.legend(fontsize=9)
fig.tight_layout();fig.savefig(OUT/'figures'/'verified_stress_recovery.png',dpi=180);plt.close(fig)
with open(OUT/'data'/'verified_static_stress.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['x_mid_m','moment_Nm','shear_N','sigma_outer_Pa','sigma_2017_Pa','tau_avg_Pa','tau_rectangular_peak_Pa','tau_effective_source_Pa']);w.writerows(zip(e.x_mid,e.moment_mid,e.shear_force,e.sigma_outer,e.sigma_source_2017,e.tau_section_average,e.tau_rectangular_max,e.tau_effective_timoshenko))

# Dynamic element stresses are also computed and summarized (Newmark comparator).
dyn=solve_dynamic(21,beam=b,dt=.001,t_end=.2,load_case='constant',integrator='newmark-average')
ds=recover_element_stress(dyn.nodes,dyn.bending_moment,b)
with open(OUT/'data'/'dynamic_stress_summary.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['time_s','peak_abs_outer_normal_Pa','peak_abs_section_average_shear_Pa','peak_abs_rectangular_shear_Pa'])
 for t,sig,tavg,tp in zip(dyn.times,ds.sigma_outer,ds.tau_section_average,ds.tau_rectangular_max):
  w.writerow([t,max(abs(sig)),max(abs(tavg)),max(abs(tp))])
print('Static maximum |outer-fiber sigma| = %.3f MPa'%(max(abs(e.sigma_outer))/1e6))
print('Static maximum |source sigma|      = %.3f MPa'%(max(abs(e.sigma_source_2017))/1e6))
print('Static maximum |section avg tau|  = %.3f MPa'%(max(abs(e.tau_section_average))/1e6))
print('Dynamic peak |outer sigma|        = %.3f MPa'%(np.max(abs(ds.sigma_outer))/1e6))
