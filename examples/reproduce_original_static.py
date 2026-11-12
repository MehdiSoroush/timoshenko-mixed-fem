"""Run the static case as specified by the original coursework parameters.

The plots are new computations, NOT imported figures from the PDF.
"""
from pathlib import Path
import csv
import matplotlib.pyplot as plt
from timoshenko_mixed import Beam, solve_static, midpoint_deflection

OUT = Path(__file__).resolve().parents[1] / "results"
(OUT / "figures").mkdir(parents=True, exist_ok=True)
(OUT / "data").mkdir(parents=True, exist_ok=True)

beam = Beam()
n_values = (5, 15, 25, 35, 55, 101)
records = []
fig, ax = plt.subplots(figsize=(8.2, 4.8))
for n in n_values:
    solution = solve_static(n_nodes=n, beam=beam)
    mid = midpoint_deflection(solution)
    peak = float(solution.deflection.max())
    records.append((n, mid, peak, solution.free_residual_l2))
    ax.plot(solution.nodes, solution.deflection, label=f"{n} nodes", linewidth=1.4)
    print(f"n={n:3} | w(L/2)={mid:.9f} m | max(w)={peak:.9f} m | free residual={solution.free_residual_l2:.2e}")
ax.set_xlabel("Beam position x (m)")
ax.set_ylabel("Beam deflection w (m)")
ax.set_title("Static fixed-pinned mixed FEM: Q at center, q on left half")
ax.grid(alpha=.3)
ax.legend()
fig.tight_layout()
fig.savefig(OUT / "figures" / "original_case_mesh_convergence.png", dpi=180)
plt.close(fig)
with open(OUT / "data" / "original_case_convergence.csv", "w", newline="", encoding="utf-8") as fh:
    writer = csv.writer(fh)
    writer.writerow(("n_nodes", "midspan_displacement_m", "max_displacement_m", "free_residual_l2"))
    writer.writerows(records)
