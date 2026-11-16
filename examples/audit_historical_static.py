"""Compare the archived stiffness-matrix sign with the verified coefficient.

This deliberately demonstrates why the printed '-a' variant cannot be used
as a convergent physical solver. Neither line is a recovered runnable source.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from timoshenko_mixed import Beam,solve_static,midpoint_deflection
from timoshenko_mixed.reference import solve_reference_static

ROOT=Path(__file__).resolve().parents[1]/'results';(ROOT/'figures').mkdir(exist_ok=True,parents=True);(ROOT/'data').mkdir(exist_ok=True,parents=True)
b=Beam();records=[]
for n in (5,9,15,25,35,55,81,101):
  le=b.length/(n-1)
  a=le**2/(3*b.EI)+1/b.kGA; bb=le**2/(6*b.EI)-1/b.kGA
  # Verbatim 2017 printed *static* last element term: -a (unlike the dynamic +a).
  ke=np.array([[0,1,0,-1],[-1,a,1,bb],[0,-1,0,1],[1,bb,-1,-a]])/le
  K=lil_matrix((2*n,2*n))
  for j in range(n-1):
    ids=np.arange(2*j,2*j+4)
    for i in range(4):
      for k in range(4):K[ids[i],ids[k]]+=ke[i,k]
  # Explicitly reproduce original static load assembly, with odd n.
  f=np.zeros(2*n)
  for j in range((n-1)//2):
    f[2*j]+=4800*le/2;f[2*(j+1)]+=4800*le/2
  f[n-1]+=3000  # original 1-based Ft(n), Python zero-based
  free=np.ones(2*n,dtype=bool);free[[0,2*n-2,2*n-1]]=False
  u=np.zeros(2*n);u[free]=spsolve(K.tocsr()[free][:,free],f[free])
  old=float(u[n-1]); fixed=midpoint_deflection(solve_static(n,beam=b));reference=solve_reference_static(n,beam=b).deflection[n//2]
  records.append((n,old,fixed,reference))
  print('n=%3d printed -a: %+9.5f mm | corrected +a: %+9.5f mm | reference: %+9.5f mm'%(n,old*1e3,fixed*1e3,reference*1e3))
with open(ROOT/'data'/'static_matrix_sign_audit.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['n_nodes','literal_printed_minus_a_midspan_m','corrected_plus_a_midspan_m','independent_reference_midspan_m']);w.writerows(records)
fig,ax=plt.subplots(figsize=(8,4.9))
for j,l,sty in [(1,'Literal printed -a (non-convergent)','o--'),(2,'Corrected +a (mixed)','s-'),(3,'Independent physical FEM','^-')]:
 ax.plot([v[0] for v in records],[v[j]*1000 for v in records],sty,label=l)
ax.set(xlabel='Number of nodes',ylabel='Signed midpoint displacement (mm)',title='2017 source inconsistency: printed static matrix term')
ax.grid(alpha=.3);ax.legend(fontsize=9);fig.tight_layout();fig.savefig(ROOT/'figures'/'static_matrix_sign_audit.png',dpi=180);plt.close(fig)
