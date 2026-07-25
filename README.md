# SPH Lid-Driven Cavity Simulation in C

A 2D fluid dynamics simulation of the classic **lid-driven cavity** problem, implemented from scratch in C using the **Smoothed Particle Hydrodynamics (SPH)** method. The fluid is discretized into Lagrangian particles enclosed in a square cavity whose top wall (lid) moves at constant velocity, driving a recirculating vortex.

Neighbor search uses a **cell linked list**, reducing the cost from $O(N^2)$ to $O(N)$.

📼 Result: [`sph_animation.mp4`](sph_animation.mp4)

> Each section below is collapsed. Click a heading to expand it.

---

<details open>
<summary><b>🚀 Build and Run</b></summary>

<br>

```bash
# Compile (link the math library explicitly with -lm)
gcc -O2 -Wall -o sph sph.c -lm

# The state files are written to ./output — create it first
mkdir -p output

# Run for N time steps (e.g. 2000 steps = 0.1 s of simulated time)
./sph 2000
```

Tested with `gcc` on Arch Linux; any C99-compliant compiler should work since only the standard library and `libm` are used.

Then animate the result:

```bash
python3 visualize_sph.py                 # writes sph_animation.mp4
python3 visualize_sph.py --color density # or color by density / pressure
python3 visualize_sph.py --step 500      # show a single step on screen
gnuplot -e "step=500" visualize_sph.gp   # single-step view with gnuplot
```

<details>
<summary>Reproducing the validation against the brute-force reference</summary>

<br>

```bash
gcc -O2 -Wall -o sph_bruteforce sph_bruteforce.c -lm
```

Run both binaries in separate directories with the same step count, then compare the final state files column by column:

```bash
paste a/output/state_0100 b/output/state_0100 | awk '
  {for(c=2;c<=12;c++){d=$c-$(c+12); if(d<0)d=-d;
   s=($c<0?-$c:$c); if(s<1e-30)s=1e-30;
   r=d/s; if(r>m)m=r}}
  END{printf "max rel diff = %.3e\n", m}'
```

Expected output: `max rel diff ≈ 1e-13`.

</details>

</details>

---

<details>
<summary><b>📐 Theoretical Background</b></summary>

<br>

<details>
<summary>The SPH Method</summary>

<br>

SPH is a mesh-free, Lagrangian numerical method in which the fluid is represented by a finite set of particles carrying mass, velocity, and thermodynamic properties. Any field $A(\mathbf{r})$ is approximated by a kernel interpolation over neighboring particles:

$$
A(\mathbf{r}_i) \approx \sum_{j} \frac{m_j}{\rho_j}\, A_j \, W(|\mathbf{r}_i - \mathbf{r}_j|, h)
$$

where $m_j$ and $\rho_j$ are the mass and density of particle $j$, $W$ is the smoothing kernel, and $h$ is the smoothing length that defines the interaction range. Spatial derivatives of the field are transferred onto the kernel, which is known analytically:

$$
\nabla A(\mathbf{r}_i) \approx \sum_{j} \frac{m_j}{\rho_j}\, A_j \, \nabla_i W_{ij}, \qquad W_{ij} \equiv W(|\mathbf{r}_i - \mathbf{r}_j|, h)
$$

This makes SPH naturally suited to problems with free surfaces, moving boundaries, and large deformations, since no mesh needs to be maintained.

</details>

<details>
<summary>Smoothing Kernel</summary>

<br>

This implementation uses the **cubic spline kernel** (Monaghan & Lattanzio, 1985) in two dimensions. With $R = r/h$:

$$
W(R, h) = \alpha_D
\begin{cases}
\dfrac{2}{3} - R^2 + \dfrac{1}{2}R^3, & 0 \le R < 1 \\[8pt]
\dfrac{1}{6}(2 - R)^3, & 1 \le R \le 2 \\[8pt]
0, & R > 2
\end{cases}
$$

where the 2D normalization constant is

$$
\alpha_D = \frac{15}{7\pi h^2}
$$

so that $\int W \, dA = 1$. The kernel has **compact support** of radius $2h$: only particles with $r_{ij} \le \kappa h_{ij}$ (with $\kappa = 2$) contribute, which is exactly the criterion used by the neighbor search — and the property the cell grid is built on. The gradient is computed component-wise as

$$
\frac{\partial W}{\partial x_i} = \frac{\partial W}{\partial R}\frac{x_{ij}}{h\, r_{ij}}
$$

with a guard against division by zero when two particles coincide ($r_{ij} < 10^{-12}$).

</details>

<details>
<summary>Density Estimation</summary>

<br>

The density of each fluid particle is computed by **kernel summation**, including the self-contribution $m_i W(0, h_i)$:

$$
\rho_i = \sum_{j} m_j W_{ij}
$$

Near boundaries the kernel support is truncated and the raw summation underestimates the density. To correct this, a **Shepard filter** (zeroth-order renormalization) is applied:

$$
\rho_i^{\text{corr}} = \frac{\sum_j m_j W_{ij}}{\sum_j \dfrac{m_j}{\rho_j} W_{ij}}
$$

The denominator equals 1 in the bulk of the fluid (partition of unity) and drops below 1 near walls, restoring a consistent density there.

</details>

<details>
<summary>Equation of State</summary>

<br>

The flow is treated as **weakly compressible**: instead of solving a pressure Poisson equation, pressure is obtained algebraically from density through a linear (isothermal-like) equation of state,

$$
p_i = c^2 \rho_i
$$

with an artificial speed of sound $c = 0.01\ \text{m/s}$. The value of $c$ is chosen to be roughly an order of magnitude smaller than a physical sound speed but large compared to the lid velocity ($v_{\text{lid}} = 1.5\times10^{-2}\ \text{m/s}$ gives a Mach number of order 1 here — see **Limitations**), keeping density fluctuations bounded while allowing a reasonable time step.

</details>

<details>
<summary>Momentum Equation</summary>

<br>

The Navier–Stokes momentum equation in Lagrangian form,

$$
\frac{d\mathbf{v}}{dt} = -\frac{1}{\rho}\nabla p + \boldsymbol{\Theta}_{\text{visc}} + \mathbf{f}_{\text{boundary}}
$$

is discretized with the **symmetric pressure-gradient form**, which conserves linear momentum exactly by construction (pairwise antisymmetric forces):

$$
\frac{d\mathbf{v}_i}{dt} = -\sum_j m_j \left( \frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2} \right) \nabla_i W_{ij}
$$

The associated rate of change of internal energy per unit mass is also accumulated:

$$
\frac{du_i}{dt} = \frac{1}{2}\sum_j m_j \left( \frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2} \right) \mathbf{v}_{ij} \cdot \nabla_i W_{ij}
$$

with $\mathbf{v}_{ij} = \mathbf{v}_i - \mathbf{v}_j$.

</details>

<details>
<summary>Artificial Viscosity</summary>

<br>

Physical viscosity is modeled with the standard **Monaghan artificial viscosity** (Monaghan, 1992), which provides shear and bulk dissipation and stabilizes the scheme. It acts only on approaching particle pairs ($\mathbf{v}_{ij}\cdot\mathbf{r}_{ij} < 0$):

$$
\Pi_{ij} =
\begin{cases}
\dfrac{-\alpha \bar{c}_{ij} \phi_{ij} + \beta \phi_{ij}^2}{\bar{\rho}_{ij}}, & \mathbf{v}_{ij}\cdot\mathbf{r}_{ij} < 0 \\[8pt]
0, & \text{otherwise}
\end{cases}
\qquad
\phi_{ij} = \frac{\bar{h}_{ij}\, \mathbf{v}_{ij}\cdot\mathbf{r}_{ij}}{|\mathbf{r}_{ij}|^2 + \varepsilon^2}
$$

where barred quantities are arithmetic means between particles $i$ and $j$, $\alpha = \beta = 1$, and $\varepsilon^2 = 0.01\,(\Delta x)^2$ prevents singularities when particles get very close. The viscous term enters the momentum equation as

$$
\left.\frac{d\mathbf{v}_i}{dt}\right|_{\text{visc}} = -\sum_j m_j \Pi_{ij} \nabla_i W_{ij}
$$

</details>

<details>
<summary>Boundary Treatment</summary>

<br>

The walls are discretized with **virtual (boundary) particles** of type `-1`, placed with half the fluid spacing along the four edges of the cavity. They contribute to the density and pressure sums of the fluid (they participate in the neighbor search), and additionally exert a short-range **Lennard-Jones-type repulsive force** (Monaghan, 1994) that prevents fluid penetration:

$$
\mathbf{f}_{ij} = D\left[ \left(\frac{r_0}{r_{ij}}\right)^{n_1} - \left(\frac{r_0}{r_{ij}}\right)^{n_2} \right] \frac{\mathbf{r}_{ij}}{r_{ij}^2}, \qquad r_{ij} < r_0
$$

with $n_1 = 12$, $n_2 = 4$, cut-off radius $r_0 = \Delta x / 2$, and strength $D = 0.01$. The force vanishes for $r_{ij} \ge r_0$ and diverges as $r_{ij} \to 0$, acting as a soft wall.

The **top lid** is a row of boundary particles with constant prescribed velocity $v_x = 1.5\times10^{-2}\ \text{m/s}$. Their velocity enters the viscous and XSPH terms of nearby fluid particles, which is the mechanism that drags the fluid and generates the cavity vortex.

</details>

<details>
<summary>XSPH Velocity Correction</summary>

<br>

To keep the particle distribution orderly and prevent particle interpenetration, the **XSPH correction** (Monaghan, 1989) moves each particle with a velocity smoothed over its neighbors:

$$
\hat{\mathbf{v}}_i = \mathbf{v}_i - \epsilon \sum_j \frac{m_j}{\bar{\rho}_{ij}} \, \mathbf{v}_{ij} \, W_{ij}
$$

with $\epsilon = 0.3$. In this implementation the correction is applied directly to the stored velocity after the force computation.

</details>

<details>
<summary>Time Integration</summary>

<br>

Time stepping uses a **leap-frog scheme in drift-kick-drift (DKD) form**, which is second-order accurate and symplectic-like:

$$
\begin{aligned}
\mathbf{r}_i^{\,n+1/2} &= \mathbf{r}_i^{\,n} + \frac{\Delta t}{2}\,\mathbf{v}_i^{\,n} & \text{(drift)} \\
\mathbf{v}_i^{\,n+1} &= \mathbf{v}_i^{\,n} + \Delta t\,\mathbf{a}_i^{\,n+1/2} & \text{(kick)} \\
\mathbf{r}_i^{\,n+1} &= \mathbf{r}_i^{\,n+1/2} + \frac{\Delta t}{2}\,\mathbf{v}_i^{\,n+1} & \text{(drift)}
\end{aligned}
$$

The internal energy $u$ is drifted alongside the position using $du/dt$. The fixed time step $\Delta t = 5\times10^{-5}\ \text{s}$ must satisfy the CFL-type condition $\Delta t \lesssim 0.25\, h / c$ for stability.

</details>

</details>

---

<details>
<summary><b>🔍 The Cell-Linked-List Neighbor Search</b></summary>

<br>

Since the cubic spline kernel vanishes for $r > \kappa h$, a particle can only interact with particles inside a disc of radius $\kappa h$. Partitioning the domain into square cells of side exactly $\kappa h$ means every neighbor is guaranteed to lie in the particle's own cell or in one of the 8 adjacent ones, so the search visits a bounded number of candidates per particle instead of all $N$ of them:

$$
\text{cost} : O(N^2) \;\longrightarrow\; O(N)
$$

The grid is held in two global arrays (Hockney & Eastwood):

```c
int *cellHead;   /* first particle of each cell, -1 if empty        */
int *cellNext;   /* next particle in the same cell, -1 at the end   */
```

This costs `nCells + nPart` integers total and imposes no per-cell particle limit.

<details>
<summary>Two details that matter for correctness</summary>

<br>

- **Cell size.** Smaller cells would force scanning more than a $3\times3$ block; larger cells would waste distance checks. $\kappa h_{\max}$ is the exact optimum, and $h_{\max}$ (not $h_i$) must be used because the criterion averages $\bar{h}_{ij} = (h_i + h_j)/2$.
- **Escaping particles.** A fluid particle may drift outside $[0, L]$. Its cell index is clamped into range rather than rejected. This stays correct because clamping is monotone: any particle within $\kappa h$ of a clamped particle is itself clamped to the same edge strip, so no pair is lost, and the distance test discards the false candidates the clamp introduces.

</details>

<details>
<summary>Validation and measured speedup</summary>

<br>

| Check | Result |
|---|---|
| `gcc -O2 -Wall -Wextra` | no warnings |
| 100 steps vs. `sph_bruteforce.c` | max relative difference $\approx 10^{-13}$ |
| `-fsanitize=address,undefined`, 20 steps | clean |
| 600 steps, wall clock | 12.3 s → **5.2 s** |

The $10^{-13}$ agreement is floating-point summation order, not physics: a single missed neighbor would show up as an $O(1)$ discrepancy within a few steps. The timing includes writing one state file per step, so the speedup of the search alone is considerably larger.

</details>

</details>

---

<details>
<summary><b>🧩 Code Structure</b></summary>

<br>

The simulation lives in a single translation unit, `sph.c`. A second file, `sph_bruteforce.c`, contains the same physics with a direct $O(N^2)$ neighbor search and is kept only as an independent reference to validate `sph.c` against.

<details>
<summary>Data layout</summary>

<br>

```c
typedef struct
{
  int    id;          /* global particle index                       */
  double pos[2];      /* position (x, y)                             */
  double vel[2];      /* velocity                                    */
  double accel[2];    /* acceleration accumulated each step          */
  double mass, rho;   /* mass and density                            */
  double h;           /* smoothing length                            */
  double p, c;        /* pressure and speed of sound                 */
  double u, du;       /* internal energy and its time derivative     */
  int    *nn;         /* neighbor indices (dynamic array)            */
  int    nNeighbors;  /* neighbor count                              */
  int    nnCap;       /* capacity reserved in the neighbor arrays    */
  double *dx, *dy;    /* per-neighbor separations x_i - x_j, y_i-y_j */
  double *r;          /* per-neighbor distances                      */
  double *W, *dWx, *dWy; /* per-neighbor kernel and gradient values  */
  int    type;        /* 1 = fluid, -1 = boundary                    */
} Particles;
```

Each particle caches its neighbor list **and** the kernel values evaluated for each neighbor, so that `density()`, `navierStokes()`, `viscosity()` and `meanVelocity()` reuse the same geometric data without recomputing kernels. The arrays are grown by **capacity doubling** (`nnCap`) and reused across steps: since the neighbor count of a particle barely changes from one step to the next, the `realloc` calls stop almost immediately and the search becomes allocation-free in steady state.

</details>

<details>
<summary>Function-by-function overview</summary>

<br>

| Function | Role |
|---|---|
| `main` | Parses the number of steps, allocates the fluid particles, calls `ics` and `buildCellGrid`, then runs the DKD loop: `buildCellList` → `NN` → `density` → `drift` → `acceleration` → `kick` → `drift` → `printState`. |
| `ics` | Builds the initial conditions: a `40 × 40` lattice of fluid particles in a `1 mm × 1 mm` box, then appends the four walls of boundary particles via `realloc`, spaced at $\Delta x / 2$. The top wall carries the lid velocity. Dumps every wall and the full initial state to `.output` files. |
| `W`, `dW` | Cubic spline kernel and the Cartesian components of its gradient, with the 2D normalization $15 / (7\pi h^2)$ and a division-by-zero guard in `dW`. |
| `testKernel` | Samples $W$ and $dW$ over $r \in [-3, 3]$ into `kernel_test.output` for plotting/validation. |
| `buildCellGrid` | Allocates the cell grid once, with cell side equal to the support radius $\kappa h_{\max} = 2h$, giving a $20 \times 20$ grid for the default parameters. |
| `buildCellList` | Rebuilds the `cellHead`/`cellNext` linked lists each step from the current positions, in $O(N)$. Must run before `NN`, since the particles moved in the previous step. |
| `NN` | Neighbor search restricted to the $3 \times 3$ block of cells around the particle's own cell, applying the compact-support criterion $r_{ij} \le \kappa \bar{h}_{ij}$ with $\kappa = 2$ (plus a cheap bounding-box reject before the square root). Scans both fluid and boundary particles and stores index, separations, distance, kernel and gradient per neighbor. |
| `test_NN` | Picks 20 random fluid particles and dumps them with their neighbors to `NN_test.output`, allowing a visual check of the search (run only on the first step). |
| `density` | Kernel-summation density including the self term $m_i W(0, h_i)$, followed by the Shepard renormalization. |
| `eos` | Sets $c = 0.01$ and $p = c^2 \rho$ for every particle, boundary included, so wall particles exert pressure back on the fluid. |
| `navierStokes` | Resets accelerations, calls `eos`, then accumulates the symmetric pressure-gradient acceleration and the corresponding $du/dt$. |
| `viscosity` | Adds the Monaghan artificial-viscosity acceleration and its energy contribution for approaching pairs only. |
| `boundaryInteraction` | Adds the Lennard-Jones-type repulsive force from boundary neighbors closer than $r_0 = \Delta x / 2$. |
| `meanVelocity` | Applies the XSPH correction with $\epsilon = 0.3$. |
| `acceleration` | Orchestrates the four contributions above in sequence. |
| `drift`, `kick` | The two halves of the leap-frog integrator: position/energy half-updates and the full velocity update. |
| `printState` | Writes the full particle state (id, position, velocity, acceleration, density, mass, pressure, sound speed, energy) to `./output/state_XXXX`, one file per step, and aborts with a clear message if `./output` does not exist. |

</details>

</details>

---

<details>
<summary><b>⚙️ Simulation Parameters</b></summary>

<br>

| Parameter | Symbol | Value |
|---|---|---|
| Domain size | $L_x \times L_y$ | $10^{-3} \times 10^{-3}\ \text{m}$ |
| Fluid particles | $n_x \times n_y$ | $40 \times 40 = 1600$ |
| Boundary particles | — | ~320 (four walls, spacing $\Delta x/2$) |
| Particle spacing | $\Delta x$ | $2.5\times10^{-5}\ \text{m}$ |
| Smoothing length | $h$ | $\Delta x$ |
| Support radius | $\kappa h$ | $2h$ |
| Cell grid | $n_{cx} \times n_{cy}$ | $20 \times 20$ cells of side $\kappa h$ |
| Reference density | $\rho_0$ | $1000\ \text{kg/m}^3$ (water) |
| Lid velocity | $v_{\text{lid}}$ | $1.5\times10^{-2}\ \text{m/s}$ |
| Artificial sound speed | $c$ | $0.01\ \text{m/s}$ |
| Viscosity coefficients | $\alpha,\ \beta$ | $1.0,\ 1.0$ |
| XSPH coefficient | $\epsilon$ | $0.3$ |
| Boundary force | $D,\ n_1,\ n_2,\ r_0$ | $0.01,\ 12,\ 4,\ \Delta x/2$ |
| Time step | $\Delta t$ | $5\times10^{-5}\ \text{s}$ |

</details>

---

<details>
<summary><b>📄 Output Files</b></summary>

<br>

| File | Content |
|---|---|
| `output/state_XXXX` | Full particle state at step `XXXX`: `id x y vx vy ax ay rho m p c u`. |
| `fluid_ics.output` | Initial state of all particles (fluid + boundary). |
| `bottom/right/top/left_border.output` | Initial state of each wall, for plotting the geometry. |
| `kernel_test.output` | Sampled kernel and gradient, for validating $W$ and $\nabla W$. |
| `NN_test.output` | 20 random particles with their neighbor lists, for validating the search. |

The `state_XXXX` files are plain ASCII and can be animated directly with `gnuplot`, or with Python (`numpy.loadtxt` + `matplotlib` + `ffmpeg`). All of them are regenerated by a run, so they are excluded from version control.

</details>

---

<details>
<summary><b>🚧 Limitations and Possible Improvements</b></summary>

<br>

- **Single-threaded.** The neighbor search and the force loops are now $O(N)$ and trivially parallel over particles; OpenMP over the `i` loop would be the next speedup.
- **Grid rebuilt every step.** `buildCellList` is $O(N)$ and cheap, but particles move far less than a cell per step, so an incremental update would avoid most of the work.
- **Sound speed vs. lid velocity.** With $c = 0.01$ and $v_{\text{lid}} = 0.015\ \text{m/s}$ the Mach number is $\sim 1.5$, violating the weakly compressible assumption ($\text{Ma} \lesssim 0.1$); increasing $c$ (and reducing $\Delta t$ accordingly) would yield a more incompressible flow.
- **Artificial viscosity only.** A physical laminar viscosity term (e.g. Morris et al., 1997) would allow matching a target Reynolds number and comparing against the reference solution of Ghia et al. (1982).
- **Fixed $\Delta t$.** An adaptive CFL-based time step would improve robustness.

</details>

---

<details>
<summary><b>📚 References</b></summary>

<br>

- Monaghan, J. J. (1992). *Smoothed Particle Hydrodynamics*. Annual Review of Astronomy and Astrophysics, 30, 543–574.
- Monaghan, J. J. (1994). *Simulating Free Surface Flows with SPH*. Journal of Computational Physics, 110(2), 399–406.
- Monaghan, J. J., & Lattanzio, J. C. (1985). *A refined particle method for astrophysical problems*. Astronomy and Astrophysics, 149, 135–143.
- Liu, G. R., & Liu, M. B. (2003). *Smoothed Particle Hydrodynamics: A Meshfree Particle Method*. World Scientific.
- Ghia, U., Ghia, K. N., & Shin, C. T. (1982). *High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method*. Journal of Computational Physics, 48(3), 387–411.
- Hockney, R. W., & Eastwood, J. W. (1988). *Computer Simulation Using Particles*. Adam Hilger. (cell linked lists)

</details>
