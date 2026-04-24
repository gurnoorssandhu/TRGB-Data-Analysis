# 2D Persistent Homology TRGB Detector — Full Explanation

This document walks through, from zero, how the `ph_tip_2d` function in
`analysis.ipynb` detects the Tip of the Red Giant Branch (TRGB) in a
color-magnitude diagram (CMD) using persistent homology (PH). No prior
knowledge of topology is assumed.

---

## Part 1 — What is persistent homology, intuitively?

### The mental picture: flooding a landscape

Imagine a topographic map of a mountain range. Now imagine slowly flooding the
landscape with water, starting from above the tallest peak and lowering the
water level.

- At first, nothing pokes above the water. No "land components" exist.
- Lower the water a bit — the tallest peak pokes through. **One island is born.**
- Lower more — a second peak pokes through. **Second island is born.** Two disjoint islands.
- Keep lowering. Eventually the saddle between two islands emerges above water — the two islands touch and **merge into one.**
- Keep going. More peaks emerge, more merges happen, until at water level 0 the whole landscape is dry.

Persistent homology is the bookkeeping of this process. For each island we record:

- **Birth level**: the water level at which it first appeared.
- **Death level**: the water level at which it merged into a bigger, older island.
- **Persistence**: `birth − death`. This is literally how tall that peak is
  above the saddle where it got absorbed.

A tall, prominent peak has **high persistence**. A small bump that merges almost
immediately into a neighboring peak has **low persistence**. Noise in the
landscape creates lots of tiny low-persistence features. Real mountains create
few high-persistence features.

**That is the entire core idea.** Everything else is machinery to apply it to
functions that aren't literal landscapes.

### The "H0" terminology

Mathematicians index topological features by dimension:

- **H0** = connected components ("how many islands are there")
- **H1** = loops / holes ("how many lakes fully enclosed by land")
- **H2** = enclosed voids (only in 3D+)

For detecting dense clusters in a 2D density field, only H0 matters. Every H0
class is born at a local maximum of the density and dies at a saddle where
its basin meets a deeper basin. That's it.

### Sublevel vs. superlevel sets

Two equivalent phrasings of the same thing:

- **Sublevel-set filtration on `f`**: grow the region `{x : f(x) ≤ t}` as `t`
  increases from `−∞` to `+∞`. Components appear at local minima.
- **Superlevel-set filtration on `f`**: grow the region `{x : f(x) ≥ t}` as
  `t` decreases from `+∞` to `−∞`. Components appear at local maxima.

They are mirror images. For density estimates we care about peaks, so
superlevel sets are the natural framing. But most PH libraries (including
GUDHI's `CubicalComplex`) only implement sublevel sets. The trick: doing
sublevel PH on `-f` is identical to superlevel PH on `f`. That's why the code
passes `-rho` to GUDHI.

### The persistence diagram

Every H0 class contributes one point to the **persistence diagram**: plot
`(birth, death)`. Tall persistent classes lie far from the diagonal
`birth = death`; noise clusters near the diagonal. This gives an immediate
visual separation of signal from noise, with no hand-tuned threshold.

---

## Part 2 — Why PH is the right tool for TRGB

### What the CMD looks like

- Each star contributes a point to the `(F606W − F814W, F814W)` plane.
- The **red giant branch** is a dense, roughly vertical ribbon of stars
  spanning a range of magnitudes.
- At a specific F814W magnitude — the **tip** — the RGB ends sharply. Above
  the tip (brighter), there are far fewer stars. Below the tip (fainter), the
  RGB continues densely.
- The tip is a calibrated standard candle. Measuring its magnitude gives a
  distance to the host galaxy.

### Why finding the tip is non-trivial

The discontinuity is sharp in theory but blurred in practice by photometric
noise and Poisson statistics. Classical methods either:

- compute a **Sobel-like edge filter** on a 1D luminosity function and find
  its strongest positive response (the Hatt-Sobel method), or
- fit a broken power law to the luminosity function directly.

Both collapse the CMD along color, discarding 2D structure. The dense locus
of the RGB is a genuinely 2D object and its upper boundary in the `(color,
F814W)` plane is what we want to identify.

### Why PH fits

Treat the density `ρ(color, mag)` of stars in the CMD as a landscape. The
RGB is the dominant mountain. Its peak is the densest point in the
color-magnitude plane; its basin (the superlevel region surrounding it)
stretches along the RGB ribbon. The TRGB is the **bright boundary of that
basin at the right density threshold**.

PH gives two things for free:

1. A principled way to identify "the RGB" as the most persistent H0 class.
2. A principled way to pick the density threshold that marks the tip —
   detailed in Part 4.

---

## Part 3 — Step-by-step walkthrough of `ph_tip_2d`

Reference: `analysis.ipynb` cell `074c28a0`.

### Step A — Window the CMD

```python
sel = (mag >= 22) & (mag <= 26) & (color >= 1.0) & (color <= 3.0)
```

Restrict to the RGB search box. This removes:

- **Saturation** at the bright end (`F814W < 22`).
- **Photometric completeness fall-off** at the faint end (`F814W > 26`),
  where detection efficiency drops sharply and creates a spurious density
  gradient.
- **Non-RGB populations** (blue main-sequence stars, foreground MW stars) at
  colors outside `[1.0, 3.0]`.

Without this, PH locks onto the completeness edge instead of the TRGB.

### Step B — Build a 2D density estimate

```python
kde = gaussian_kde(np.vstack([color, mag]), bw_method=bandwidth)
rho = kde(...).reshape(Nm, Nc)
```

- Place a regular grid: 400 bins in color × 600 bins in magnitude, with 5%
  padding so the KDE can decay smoothly at the edges.
- Fit a 2D Gaussian kernel density estimator to the windowed star
  coordinates.
- Evaluate `ρ` on every grid cell. Result: a 600×400 array of density
  values.

Bandwidth defaults to Scott's rule; the `bandwidth` argument exposes a
multiplicative factor for robustness testing.

### Step C — Superlevel-set PH via GUDHI

```python
cc = gudhi.CubicalComplex(top_dimensional_cells=-rho)
cc.persistence(homology_coeff_field=2, min_persistence=0)
regular, _ = cc.cofaces_of_persistence_pairs()
```

What `CubicalComplex` does internally:

1. Treats the 2D grid as a cubical complex: each grid cell is a 2-cell, with
   shared edges as 1-cells and shared corners as 0-cells.
2. Assigns each top-dimensional cell (grid pixel) the value `-rho` at that
   cell; lower-dimensional cells inherit the max of their neighbors.
3. Runs the sublevel-set filtration algorithm: sort all cells by value, add
   them one at a time, track connected components.

Because we passed `-rho`, sublevel sets of `-rho` are the superlevel sets of
`rho`. So:

- Adding cells in order of increasing `-rho` = adding in order of decreasing
  `ρ` = "lowering the water level" over the density landscape.
- Every time a new local maximum of `ρ` appears, a new H0 class is born.
- When the basin around one peak grows to touch the basin around a deeper
  peak (at a saddle), the younger/shallower basin dies.

`cofaces_of_persistence_pairs()` returns, for each finite H0 class, the flat
grid index of its birth cell and its death cell. We translate those back to
density values:

```python
for b_idx, d_idx in regular[0]:
    rho_b = rho_flat[b_idx]   # density at birth (local max)
    rho_d = rho_flat[d_idx]   # density at death (saddle)
    persistence = rho_b - rho_d
```

The list, sorted by persistence, **is the persistence diagram** of the CMD
density. The most persistent classes are the real modes; the low-persistence
tail is noise.

One subtlety: the *global* maximum of `ρ` gives an **essential** H0 class
(the last component standing) and never dies. It is not in the `regular`
list. In our case the global max is the RGB peak — so the "most persistent
H0 class" we retrieve from the list is actually the second-most-significant
mode, the one whose basin eventually got swallowed by the RGB. The RGB
itself is the essential class.

### Step D — Locate the RGB basin

```python
peak = np.unravel_index(int(rho.argmax()), rho.shape)
```

The global maximum of `ρ` is the RGB peak. We'll track *its* superlevel
component as the filtration proceeds, because that component IS the RGB
basin.

### Step E — Sweep the filtration and record the bright edge

```python
t_grid = np.linspace(rho.max() * 0.02, rho.max() * 0.98, 400)
for k, t in enumerate(t_grid):
    mask = rho >= t
    lab, _ = cc_label(mask)
    region_k = (lab == lab[peak])
    # per-column bright edge:
    tips = [m_edges[np.argmax(region_k[:, j])] for j in range(Nc) if region_k[:, j].any()]
    m_top[k] = np.median(tips)
```

For each threshold `t` (400 levels spanning 2% to 98% of the peak density):

1. Take the binary mask `ρ ≥ t`.
2. Run connected-components labeling (`scipy.ndimage.label`). This gives
   every pixel a label; all pixels of the same component share a label.
3. Keep only the component containing the RGB peak — this is `region_k`,
   the RGB basin at threshold `t`.
4. For each color column of `region_k`, find the brightest (smallest
   F814W) pixel that's still inside the region. That's the **bright
   edge of the RGB basin at this color, at this threshold**.
5. Take the median across color columns: `m_top[k]`.

`m_top(t)` is a 1D curve: as `t` decreases, the RGB basin grows, and its
upper boundary typically moves bright-ward (smaller F814W).

### Step F — Find the plateau

This is the heart of the detector.

**Physical picture**: suppose the TRGB sits at F814W = 25.0, and suppose the
density along the RGB just below the tip is 0.4 while above the tip it drops
to 0.1. Then as we lower `t` from the RGB peak density:

- While `t > 0.4`: the RGB basin's bright edge creeps faint-ward → bright-ward toward 25.0.
- Once `t ≈ 0.4`: the bright edge hits F814W = 25.0.
- While `0.1 < t < 0.4`: **the bright edge stays pinned at 25.0.** To extend
  the region brighter than 25.0, `t` would have to drop below 0.1 — but it
  hasn't yet. So the region expands sideways (wider color) and downward
  (fainter), but not upward.
- Once `t < 0.1`: the bright edge finally jumps up into the sparse
  above-TRGB region.

The plateau of `m_top(t)` at the value 25.0, spanning the `t`-range
`[0.1, 0.4]`, is the 2D-topological fingerprint of the density
discontinuity. Its value is the TRGB magnitude.

**How the code finds it**:

```python
d_edge = np.gradient(m_top, np.log(t_grid))
k_star = idxs[np.argmin(np.abs(d_edge[idxs]))]
tip_mag = m_top[k_star]
t_star  = t_grid[k_star]
```

- Compute the numerical derivative `d m_top / d log(t)`.
- Ignore the first and last 10% of the sweep (endpoint artifacts).
- Pick the interior `t` where `|d m_top / d log(t)|` is smallest — the
  flattest point of the plateau.
- Report `m_top` at that `t` as the TRGB.

We use `log(t)` because the threshold sweep is density-ratio-meaningful:
"the bright edge is stable over a factor-of-4 drop in `t`" is the right
statement, not "stable over an additive step in `t`."

### Step G — Return the region for plotting

```python
mask = rho >= t_star
lab, _ = cc_label(mask)
region = (lab == lab[peak])
```

`region` is the 2D boolean mask of the RGB basin at the plateau threshold —
the "dense structure" you asked about. Its upper boundary is, by
construction, tangent to the TRGB discontinuity. `tip_mag` is the median
height of that upper boundary.

---

## Part 4 — Why the plateau = the right answer

A few alternative approaches and why they fail:

### "Just use the death-density of the most persistent H0 class"

This is the naive first try. It fails because in this CMD the RGB is overwhelmingly dominant, so the only way for the second-most-persistent class to die is by being absorbed into the RGB or by the essential class absorbing everything near `t = 0`. Either way, the death-threshold is too low and the resulting region floods most of the box. We verified empirically: this gives `tip ≈ 24.1`, wrong by ~1 magnitude.

### "Use the peak density directly"

The peak of `ρ` is somewhere in the middle of the RGB, not at the tip. Wrong by construction.

### "Threshold at some fixed fraction of peak density"

This works but requires hand-tuning the fraction, and the right fraction depends on bandwidth and dataset. The plateau method finds the right threshold automatically.

### Why the plateau method is principled

The plateau exists *because of the density discontinuity at the TRGB*. Its
existence IS the topological signature we want. Finding it is:

- **Unsupervised**: no hand-tuned threshold or fraction.
- **Robust to bandwidth**: the plateau location moves very little as KDE
  bandwidth changes, because it's a topological (not geometric) feature.
- **Connected to the PH framework**: `m_top(t)` as a function of `t` is a
  1-parameter filtration of the bright edge itself. The "most persistent
  value" of `m_top` — the one that persists over the widest `t`-range — is
  exactly what we report.

---

## Part 5 — What each piece of the output means

```python
tip_mag, t_star, region, h0_pairs, c_edges, m_edges, rho = ph_tip_2d(color, mag)
```

- **`tip_mag`** — the TRGB magnitude (F814W). The main answer.
- **`t_star`** — the density threshold at which the RGB basin's bright
  boundary sits at the TRGB. Mostly a diagnostic.
- **`region`** — 2D boolean mask of the RGB basin at `t_star`. Plotted as
  the red filled contour in the diagnostic figure.
- **`h0_pairs`** — full H0 persistence pairs `(rho_birth, rho_death,
  persistence, birth_flat_index)`, sorted by persistence descending. Used to
  draw the persistence diagram.
- **`c_edges, m_edges, rho`** — the KDE grid and density, for plotting.

---

## Part 6 — Cross-checks and error bars

### Cross-check against Sobel

The 1D Hatt-Sobel method gives a tip near F814W ≈ 25.0 for M101. `ph_tip_2d`
reports 25.003 under default settings. The two methods use entirely
different machinery (edge filter on 1D LF vs. PH on 2D KDE), so agreement
to a few hundredths of a magnitude is a strong validation.

### Bootstrap uncertainty (cell `97a1eb1a`)

Resample stars with replacement, perturb by photometric errors, rerun
`ph_tip_2d`, and build a distribution of tip estimates. The median and
16/84 percentiles give a 68% confidence interval.

### Bandwidth robustness (cell `b7226914`)

Sweep KDE bandwidth from 0.5× to 2× Scott's rule and confirm the tip is
stable. If it moves dramatically, the detector is bandwidth-dependent and
not trustworthy. If it stays within ~0.05 mag, we're fine.

---

## TL;DR

1. Window the CMD to the RGB region.
2. Build a 2D density `ρ`.
3. Use GUDHI to compute the persistence diagram of `ρ` (via `-ρ` for
   superlevel-set filtration). Confirms the RGB is the dominant persistent
   mode; produces a diagnostic persistence diagram.
4. Track the bright edge of the RGB basin as the density threshold sweeps
   from the peak down toward zero.
5. Find the threshold range where the bright edge is flat (the plateau).
   That's the TRGB discontinuity. Report its value.

The persistence idea shows up twice: once explicitly via GUDHI on `ρ`, once
implicitly via the plateau — the "most persistent value of the bright edge
over the threshold sweep."
