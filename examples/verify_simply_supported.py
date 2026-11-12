"""Independent analytical validation under a uniform distributed load.

This is NOT the original left-half load configuration; it is a benchmark.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from timoshenko_mixed import Beam, solve_static, midpoint_deflection, simply_supported_uniform_midspan

OUT = Path(__file__).resolve().parents[1] / "results"
(OUT / "data").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
beam = Beam()
q = 4800.0
exact = simply_supported_uniform_midspan(beam, q)
values = []
for n in (3, 5, 9, 17, 33, 65, 129):
    sol = solve_static(n_nodes=n, beam=beam, boundary="pinned-pinned",
                       distributed_intensity=q, distributed_end=beam.length,
                       point_force=0.0)
    numerical = midpoint_deflection(sol)
    relative_error = abs(numerical - exact) / abs(exact)
    values.append((n, numerical, exact, relative_error))
    print(f"n={n:3} | numerical={numerical:.9f} m | exact={exact:.9f} m | rel. error={relative_error:.5%}")
fig, ax = plt.subplots(figsize=(7.5, 4.5))
ax.loglog([row[0] - 1 for row in values], [row[3] for row in values], "o-")
ax.set_xlabel("Number of elements")
ax.set_ylabel("Relative midspan deflection error")
ax.set_title("Analytical Timoshenko benchmark (uniform load)")
ax.grid(alpha=.35, which="both")
fig.tight_layout()
fig.savefig(OUT / "figures" / "analytical_validation_convergence.png", dpi=180)
plt.close(fig)
with open(OUT / "data" / "analytical_uniform_benchmark.csv", "w", newline="", encoding="utf-8") as fh:
    writer = csv.writer(fh)
    writer.writerow(("n_nodes", "numerical_midspan_m", "analytical_midspan_m", "relative_error"))
    writer.writerows(values)
