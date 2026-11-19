# Important Equations: Critical Current of a Josephson Junction in a Parallel Field

## 1. Setup and the sine-Gordon equation

Take a Josephson junction with current flowing in the $x$ direction. The barrier occupies $-d/2 < x < d/2$, and the phase varies along the barrier in the $y$ direction. A magnetic field is applied in the $z$ direction:

$$
\mathbf{B} = B_z(y)\,\hat{\mathbf{z}}, \qquad
B_z = (\nabla\times\mathbf{A})_z = \frac{\partial A_y}{\partial x} - \frac{\partial A_x}{\partial y}.
$$

Here $q$ is the charge of the superconducting carriers (Cooper pairs, $q = 2e$), and the flux quantum is $\phi_0 = h/q$.

**Supercurrent (Ginzburg–Landau).**

$$
\mathbf{J}_s = \frac{q|\psi|^2}{m}\left(\hbar\nabla\phi - q\mathbf{A}\right).
$$

**Josephson relation and gauge-invariant phase difference.**

$$
J_x(y) = J_{\mathrm{DJ}} \sin\Delta\phi(y),
\qquad
\Delta\phi(y) = \phi(+d/2, y) - \phi(-d/2, y) - \frac{q}{\hbar}\int_{-d/2}^{+d/2} A_x(x,y)\,dx .
$$

$J_{\mathrm{DJ}}$ is the amplitude of the Josephson current density. $d$ is the **effective magnetic thickness** of the junction, $d = t + 2\lambda_L$, where $t$ is the barrier thickness and $\lambda_L$ is the London penetration depth. See the note on screening below.

### Phase gradient along the barrier

Differentiate $\Delta\phi$ with respect to $y$ (suppressing the $y$ argument) and use $\partial_y A_x = \partial_x A_y - B_z$:

$$
\begin{aligned}
\frac{\partial(\Delta\phi)}{\partial y}
&= \frac{\partial\phi(d/2)}{\partial y} - \frac{\partial\phi(-d/2)}{\partial y}
 - \frac{q}{\hbar}\int_{-d/2}^{d/2}\left(\frac{\partial A_y}{\partial x} - B_z\right)dx \\[4pt]
&= \left[\frac{\partial\phi(d/2)}{\partial y} - \frac{q}{\hbar}A_y(d/2)\right]
 - \left[\frac{\partial\phi(-d/2)}{\partial y} - \frac{q}{\hbar}A_y(-d/2)\right]
 + \frac{q}{\hbar}\int_{-d/2}^{d/2} B_z\,dx .
\end{aligned}
$$

Assume the field is uniform across the barrier thickness, so the last term becomes $\frac{qd}{\hbar}B_z$.

The brackets are proportional to the supercurrent flowing *parallel* to the barrier, $J_{s,y} \propto \hbar\,\partial_y\phi - qA_y$. Neglecting this current relative to the current across the barrier sets each bracket to zero:

$$
\hbar\frac{\partial\phi}{\partial y} = qA_y .
$$

This gives the phase gradient in terms of the local field:

$$
\boxed{\frac{\partial(\Delta\phi)}{\partial y} = \frac{qd}{\hbar}\,B_z(y) = \frac{2\pi d}{\phi_0}B_z(y)}
$$

> **Note on screening.** Screening currents do flow parallel to the barrier within roughly $\lambda_L$ of the interfaces. Using $d = t + 2\lambda_L$ rather than the bare barrier thickness absorbs their effect, since the field penetrates this far into each electrode. The approximation is then that the field is uniform over $d$ and that no further parallel current flows.

### Relation between vector potential $A_x(y)$ and phase difference $\Delta\phi(y)$

In the London gauge with negligible variation of $A_y$ across the barrier ($\partial_x A_y \approx 0$), the magnetic field is related to $A_x$ simply by

$$
B_z(y) = -\frac{\partial A_x}{\partial y} .
$$

Combining this with the phase-gradient relation $\partial_y(\Delta\phi) = \frac{qd}{\hbar}B_z(y)$ gives

$$
\frac{\partial A_x}{\partial y} = -\frac{\hbar}{qd}\frac{\partial(\Delta\phi)}{\partial y} = -\frac{\phi_0}{2\pi d}\frac{\partial(\Delta\phi)}{\partial y} .
$$

Integrating from the center of the barrier ($y = 0$) with the antisymmetric gauge choice $A_x(0) = 0$ yields the $x$-component of the vector potential along the barrier directly in terms of $\Delta\phi(y)$:

$$
\boxed{A_x(y) = -\int_0^y B_z(y')\,dy' = -\frac{\hbar}{qd}\bigl(\Delta\phi(y) - \Delta\phi(0)\bigr) = -\frac{\phi_0}{2\pi d}\bigl(\Delta\phi(y) - \Delta\phi(0)\bigr)}
$$

### The sine-Gordon equation

Maxwell–Ampère gives $(\nabla\times\mathbf{B})_x = \partial_y B_z = \mu_0 J_x$. Differentiating the boxed result and substituting the Josephson relation:

$$
\frac{\partial^2(\Delta\phi)}{\partial y^2} = \frac{qd\mu_0}{\hbar}J_{\mathrm{DJ}}\sin\Delta\phi = \frac{1}{\lambda_J^2}\sin\Delta\phi,
$$

where the **Josephson penetration depth** is

$$
\lambda_J = \sqrt{\frac{\hbar}{q\,d\,\mu_0 J_{\mathrm{DJ}}}} = \sqrt{\frac{\phi_0}{2\pi d\mu_0 J_{\mathrm{DJ}}}} .
$$

---

## 2. Narrow junctions ($w \ll \lambda_J$)

Define the critical current density as the maximum, over the free phase offset, of the average current density across the junction width:

$$
J_c = \max_{\Delta\phi_0}\;\frac{1}{w}\int_{-w/2}^{w/2} J_{\mathrm{DJ}}\sin\Delta\phi(y)\,dy .
$$

For $w \ll \lambda_J$, the junction's own (self-)field and screening currents are negligible, and the field is the uniform applied field $B_z$. The phase gradient is then constant:

$$
\Delta\phi(y) = \frac{qdB_z}{\hbar}\,y + \Delta\phi_0 = \frac{2\pi d B_z}{\phi_0}\,y + \Delta\phi_0,
$$

where $\Delta\phi_0$ is an integration constant.

Write $k = 2\pi dB_z/\phi_0$. Expanding the sine:

$$
\sin(ky + \Delta\phi_0) = \sin(ky)\cos\Delta\phi_0 + \cos(ky)\sin\Delta\phi_0 .
$$

The $\sin(ky)$ term is odd in $y$ and integrates to zero over the symmetric interval. The average is therefore $\sin\Delta\phi_0$ times the integral of $\cos(ky)$, and the maximum over $\Delta\phi_0$ occurs at $\Delta\phi_0 = \pm\pi/2$ (taking the positive value):

$$
J_c = J_{\mathrm{DJ}}\left|\frac{1}{w}\int_{-w/2}^{w/2}\cos(ky)\,dy\right|
= J_{\mathrm{DJ}}\left|\frac{\sin(kw/2)}{kw/2}\right| .
$$

With $kw/2 = \pi w d B_z/\phi_0$:

$$
\boxed{J_c(B_z) = J_{\mathrm{DJ}}\left|\frac{\sin\!\left(\pi w d B_z/\phi_0\right)}{\pi w d B_z/\phi_0}\right|}
$$

This is the classic **Fraunhofer diffraction pattern**. The first zero occurs when the flux through the junction area, $B_z w d$, equals one flux quantum $\phi_0$.

### Spatial profiles $J_x(y)$ and $A_x(y)$ in the narrow limit

At the critical-current phase offset $\Delta\phi_0 = \pi/2$, substituting $\Delta\phi(y) = ky + \pi/2$ into the Josephson relation gives the spatial supercurrent profile along the barrier:

$$
\boxed{J_x^{\mathrm{narrow}}(y) = J_{\mathrm{DJ}}\sin\!\left(ky + \frac{\pi}{2}\right) = J_{\mathrm{DJ}}\cos(ky) = J_{\mathrm{DJ}}\cos\!\left(\frac{2\pi d B_z}{\phi_0}\,y\right)}
$$

Because self-field screening is negligible ($B_z(y) = B_{\mathrm{app}}$ uniformly across the width), integrating $\partial_y A_x = -B_{\mathrm{app}}$ with $A_x(0) = 0$ gives the linear vector potential profile:

$$
\boxed{A_x^{\mathrm{narrow}}(y) = -B_{\mathrm{app}}\,y}
$$

---

## 3. Wide junctions ($w \gg \lambda_J$)

In a wide junction, the interior is in the Meissner state. Current flows only within about $\lambda_J$ of each edge, so the total current is set by the boundary conditions.

### Edge field of a semi-infinite junction

Consider a semi-infinite junction occupying $y > 0$, with field $B_0$ applied at the edge $y = 0$. Far inside, the field decays to zero and the phase difference approaches a multiple of $2\pi$, which we take to be $0$ without loss of generality:

$$
\left.\frac{\partial\Delta\phi}{\partial y}\right|_{y=0} = \frac{2\pi d}{\phi_0}B_0,
\qquad
\left.\frac{\partial\Delta\phi}{\partial y}\right|_{y\to\infty} = 0,
\qquad
\Delta\phi(y\to\infty) \to 0 .
$$

Multiply the sine-Gordon equation by $\partial_y\Delta\phi$ and integrate once:

$$
\frac{1}{2}\left(\frac{\partial\Delta\phi}{\partial y}\right)^2 = -\frac{1}{\lambda_J^2}\cos\Delta\phi + C .
$$

Fix $C$ using the **far-field ($y\to\infty$) conditions**, where $\Delta\phi\to 0$ and $\partial_y\Delta\phi\to 0$. This gives $C = 1/\lambda_J^2$, so

$$
\frac{1}{2}\left(\frac{\partial\Delta\phi}{\partial y}\right)^2 = \frac{1}{\lambda_J^2}\bigl(1-\cos\Delta\phi\bigr)
\;\;\Rightarrow\;\;
\left|\frac{\partial\Delta\phi}{\partial y}\right| = \frac{2}{\lambda_J}\left|\sin\frac{\Delta\phi}{2}\right| .
$$

Evaluate at the edge $y = 0$ using the boundary condition (the sign of the gradient follows the sign of the applied field):

$$
\frac{2\pi d}{\phi_0}|B_0| = \frac{2}{\lambda_J}\left|\sin\frac{\Delta\phi(0)}{2}\right| .
$$

Since $|\sin(\Delta\phi(0)/2)| \le 1$, there is a **maximum edge field** that the Meissner state can sustain, reached at $\Delta\phi(0) = \pi$:

$$
\boxed{B_0 = \frac{\phi_0}{\pi d\lambda_J} = 2\mu_0 J_{\mathrm{DJ}}\lambda_J}
$$

The second form follows from the definition of $\lambda_J$.

### Single-edge spatial profile

To obtain the spatial dependence $\Delta\phi(\xi)$ as a function of distance $\xi \ge 0$ inward from an edge into the Meissner interior ($\Delta\phi(\xi\to\infty) \to 0$), take the decaying branch of the first integral:

$$
\frac{\partial\Delta\phi}{\partial\xi} = -\frac{2}{\lambda_J}\sin\frac{\Delta\phi}{2} .
$$

Separating variables and integrating from $\xi = 0$ (where $\Delta\phi(0) = \phi_{\mathrm{edge}}$) to $\xi$:

$$
\int_{\phi_{\mathrm{edge}}}^{\Delta\phi(\xi)} \frac{d(\Delta\phi')}{2\sin(\Delta\phi'/2)} = -\frac{\xi}{\lambda_J}
\;\;\Rightarrow\;\;
\left[\ln\left|\tan\frac{\Delta\phi'}{4}\right|\right]_{\phi_{\mathrm{edge}}}^{\Delta\phi(\xi)} = -\frac{\xi}{\lambda_J} .
$$

Exponentiating gives the single-edge sine-Gordon profile:

$$
\tan\frac{\Delta\phi(\xi)}{4} = \tan\frac{\phi_{\mathrm{edge}}}{4}\,e^{-\xi/\lambda_J}
\;\;\Rightarrow\;\;
\boxed{\Delta\phi(\xi) = 4\arctan\!\left(\tan\frac{\phi_{\mathrm{edge}}}{4}\,e^{-\xi/\lambda_J}\right)}
$$

### Critical current of a finite wide junction

Now take a junction of width $w$ ($y \in [-w/2, w/2]$) with the Meissner state in the interior. Define $J$ as the average current density, i.e. the total current per unit length (along $z$) divided by $w$. Integrating Ampère's law across the width:

$$
B(w/2) - B(-w/2) = \mu_0 J w .
$$

The applied field fixes the mean of the two edge fields:

$$
B(w/2) + B(-w/2) = 2B_{\mathrm{app}} .
$$

Hence the edge fields are the applied field plus or minus the transport-current (self-field) contribution:

$$
B(\pm w/2) = B_{\mathrm{app}} \pm \frac{\mu_0 J w}{2} .
$$

The critical current is reached when the field at the more heavily loaded edge reaches the maximum $B_0$:

$$
B_{\mathrm{app}} + \frac{\mu_0 J_c w}{2} = B_0
\;\;\Rightarrow\;\;
J_c = \frac{2}{w\mu_0}\bigl(B_0 - B_{\mathrm{app}}\bigr).
$$

Substituting $B_0 = 2\mu_0 J_{\mathrm{DJ}}\lambda_J$:

$$
\boxed{J_c = \frac{2}{w}\left(2\lambda_J J_{\mathrm{DJ}} - \frac{|B_{\mathrm{app}}|}{\mu_0}\right),
\qquad 0 \le |B_{\mathrm{app}}| \le B_0}
$$

**Validity.** This result applies in the Meissner state, i.e. $|B_{\mathrm{app}}| \le B_0$. The opposite edge, $B(-w/2) = 2|B_{\mathrm{app}}| - B_0$, then stays within $\pm B_0$. The expression is symmetric in the direction of the field, hence $|B_{\mathrm{app}}|$. Note that vortices become thermodynamically favourable above the lower field $B_{c1,J} = 2B_0/\pi$. The linear regime up to $B_0$ therefore describes a metastable Meissner state, assuming no vortex entry.

**Zero-field limit.**

$$
J_c(0) = \frac{4\lambda_J J_{\mathrm{DJ}}}{w} = \frac{4}{w}\sqrt{\frac{J_{\mathrm{DJ}}\,\phi_0}{2\pi d\mu_0}} .
$$

The corresponding critical current per unit length is $J_c w = 4\lambda_J J_{\mathrm{DJ}}$, independent of width. This is the signature of edge-dominated current flow.

### Spatial profiles $J_x(y)$ and $A_x(y)$ in a finite wide junction

For $w \gg \lambda_J$, the two edge tails decay exponentially into the interior with negligible overlap, so the full phase profile across $y \in [-w/2, w/2]$ is obtained by superposing the single-edge solutions from $y = +w/2$ (inward distance $\xi_+ = w/2 - y$) and $y = -w/2$ (inward distance $\xi_- = y + w/2$):

$$
\Delta\phi_{\mathrm{wide}}(y) = \Delta\phi_+(y) + \Delta\phi_-(y),
$$

with

$$
\Delta\phi_+(y) = 4\arctan\!\left(\tan\frac{\phi_+}{4}\,e^{-(w/2 - y)/\lambda_J}\right),
\qquad
\Delta\phi_-(y) = 4\arctan\!\left(\tan\frac{\phi_-}{4}\,e^{-(y + w/2)/\lambda_J}\right) .
$$

Differentiating each tail with respect to $y$ gives

$$
\frac{\partial\Delta\phi_+}{\partial y} = +\frac{2}{\lambda_J}\sin\frac{\Delta\phi_+(y)}{2},
\qquad
\frac{\partial\Delta\phi_-}{\partial y} = -\frac{2}{\lambda_J}\sin\frac{\Delta\phi_-(y)}{2} .
$$

Matching $\partial_y\Delta\phi = \frac{2\pi d}{\phi_0}B_z$ at the two boundaries $y = \pm w/2$ (using $\frac{2\pi d}{\phi_0}B_0 = \frac{2}{\lambda_J}$) fixes the edge phase values $\phi_\pm = \Delta\phi_\pm(\pm w/2)$ in terms of the edge fields $B(\pm w/2) = B_{\mathrm{app}} \pm \frac{\mu_0 J w}{2}$. Evaluated at the wide-junction low-field linear critical current $J = J_c = \max\!\left(0,\,\frac{2}{w\mu_0}(B_0 - |B_{\mathrm{app}}|)\right)$, the more heavily loaded edge reaches $\pm B_0$ while the opposite edge stays within $\pm B_0$:

$$
\sin\frac{\phi_+}{2} = \frac{B(+w/2)}{B_0},
\qquad
\sin\frac{\phi_-}{2} = -\frac{B(-w/2)}{B_0},
$$

or explicitly (within the Meissner regime $|B(\pm w/2)| \le B_0$):

$$
\phi_+ = 2\arcsin\!\left(\frac{B(+w/2)}{B_0}\right),
\qquad
\phi_- = 2\arcsin\!\left(-\frac{B(-w/2)}{B_0}\right) .
$$

> **Sign check.** In zero applied field ($B_{\mathrm{app}} = 0$) with positive transport current $J > 0$, the edge fields are $B(\pm w/2) = \pm\mu_0 J w/2$, giving $\phi_+ = \phi_- > 0$ so that both edges carry positive transport current $J_x > 0$. Conversely, for $J = 0$ and $B_{\mathrm{app}} > 0$, $B(+w/2) = B(-w/2) = B_{\mathrm{app}}$, giving $\phi_- = -\phi_+$ so that equal and opposite Meissner screening currents circulate on the two edges.

Substituting $\Delta\phi_{\mathrm{wide}}(y)$ into the Josephson relation and the vector potential relation gives the wide-junction spatial profiles:

$$
\boxed{J_x^{\mathrm{wide}}(y) = J_{\mathrm{DJ}}\sin\bigl(\Delta\phi_{\mathrm{wide}}(y)\bigr) = J_{\mathrm{DJ}}\sin\bigl(\Delta\phi_+(y) + \Delta\phi_-(y)\bigr)}
$$

$$
\boxed{A_x^{\mathrm{wide}}(y) = -\frac{\phi_0}{2\pi d}\bigl(\Delta\phi_{\mathrm{wide}}(y) - \Delta\phi_{\mathrm{wide}}(0)\bigr)}
$$

---

## 4. Dimensionless Ginzburg–Landau (TDGL) simulation units

In the dimensionless TDGL units used by `MulTDGL`:

- Lengths are measured in units of the coherence length $\xi$, with $\hbar = 1$ and $q = 1$, so the dimensionless flux quantum is $\phi_0 = 2\pi$ and $\frac{2\pi d}{\phi_0} = \frac{qd}{\hbar} = d$.
- The London penetration depth is $\lambda_L = \kappa$ (the Ginzburg–Landau parameter), and the effective magnetic permeability is $\mu_0 = 1/\kappa^2$ (from Ampère's law $\kappa^2 \nabla\times\mathbf{B} = \mathbf{J}$).
- The effective magnetic thickness is $d = t + 2\kappa$, where $t$ is the dimensionless barrier thickness.
- The Josephson penetration depth and maximum Meissner edge field simplify to

$$
\lambda_J = \frac{\kappa}{\sqrt{d\,J_{\mathrm{DJ}}}},
\qquad
B_0 = \frac{2}{d\,\lambda_J} = \frac{2}{\kappa^2}\,J_{\mathrm{DJ}}\lambda_J .
$$

### 1D Ginzburg–Landau depairing current density

As a reference scale for $J_{\mathrm{DJ}}$, we use the 1D bulk Ginzburg–Landau depairing current density. For a uniform superconductor ($\alpha = 1$, $\beta = 1$, $m = 1$) carrying superfluid velocity $v_s = \partial_x\phi - A_x$, the stationary GL equation is

$$
(1 - |\psi|^2 - v_s^2)\psi = 0
\;\;\Rightarrow\;\;
|\psi|^2 = 1 - v_s^2,
$$

and the supercurrent density is

$$
J_s(v_s) = |\psi|^2 v_s = v_s(1 - v_s^2) .
$$

Maximizing with respect to $v_s$ gives $\frac{dJ_s}{dv_s} = 1 - 3v_s^2 = 0 \Rightarrow v_s = \frac{1}{\sqrt{3}}$, yielding the dimensionless 1D depairing current density:

$$
\boxed{J_{\mathrm{DJ}} = J_{\mathrm{dp}} = \frac{2}{3\sqrt{3}} \approx 0.3849}
$$


### 1D Junction Depairing Current Density

Blair derived an expression for the maximum current that can flow across a junction at a given point, extending Fink's zero-field solution to low-field.
Blair begins by taking the contour integral of the GI phase. By applying Stokes' theorem, he obtains

$$
\oint \nabla\gamma \cdot d\mathbf{l}
=
\oint \nabla\theta \cdot d\mathbf{l}
-
\int \tilde{\mathbf{B}} \cdot d\mathbf{S}.
\tag{2.5.1}
$$

To ensure that the order parameter is single-valued, the contour integral of $\theta$ (first term on the RHS) is equal to $2\pi n$, where $n$ is the number of vortices within the contour. In the very narrow limit, one can assume there are no fluxons in the system ($n=0$), and so this term has no contribution. He assumes that the order parameter magnitude is symmetric about both the $y$ and $x$-axes, and that the screening currents, and hence $\partial_y\gamma$, are antisymmetric about the $y$ and $x$-axes. To first order, the transport current is uniform in the $y$-direction; there is a laminar flow of transport current in the $x$-direction. Since the supercurrent in the $y$-direction is assumed to be zero, the $x$-component of the supercurrent is,

$$
\tilde{J}_{s:x}
=
\frac{1}{\tilde{m}_n}\,|\psi|^2\,\partial_x\gamma(\bar{y}=0).
\tag{2.5.2}
$$

Inserting this back into Eq. (2.5.1), applying Stokes' theorem for a second time, and taking the length of the contours parallel to the $x$-axis to be sufficiently short leads to the result

$$
\partial_x\gamma(\bar{y})
=
\frac{\tilde{J}_{s:x}}{|\psi|^2}
+
\tilde{B}_{\mathrm{app}}\,\bar{y}.
\tag{2.5.3}
$$

Since the system is uniform in the $y$-direction, Blair utilises the mean value theorem and takes an average over the width of the system, replacing quantities by their averaged values:

$$
\tilde{F}
=
\frac{1}{\tilde{w}_s}
\int_{-\tilde{w}_s/2}^{\tilde{w}_s/2}
|\psi|\,d\bar{y},
\tag{2.5.4}
$$

$$
\langle \tilde{J}_{s:x} \rangle_y
=
\frac{1}{\tilde{w}_s}
\int_{-\tilde{w}_s/2}^{\tilde{w}_s/2}
\tilde{J}_{s:x}\,d\bar{y}
\equiv
\langle \tilde{J}_{\mathrm{app}} \rangle_y.
\tag{2.5.5}
$$

Similarly, the average of the square of the result found in Eq. (2.5.3) is given as

$$
\left\langle
\left( \partial_x\gamma(\bar{y}) \right)^2
\right\rangle
=
\left(
\frac{\tilde{B}_{\mathrm{app}}\tilde{w}_s}{\sqrt{12}}
\right)^2
+
\frac{\tilde{m}_n^2
\left( \langle \tilde{J}_{\mathrm{app}} \rangle_y \right)^2}
{\tilde{F}^4}.
\tag{2.5.6}
$$

Inserting this result back into Eq. (2.3.42) yields a single equation,

$$
\frac{1}{\tilde{m}_n}
\frac{\partial^2 \tilde{F}}{\partial \tilde{x}^2}
+
\left[
\left( \tilde{\alpha}_n - \frac{q^2}{\tilde{m}_n} \right)
-
\tilde{\beta}_n \tilde{F}^2
-
\frac{\tilde{m}_n
\left( \langle \tilde{J}_{\mathrm{app}} \rangle_y \right)^2}
{\tilde{F}^4}
\right]
\tilde{F}
= 0,
\tag{2.5.7}
$$

where $q^2$ has been defined as

$$
q^2
=
\left(
\frac{\tilde{B}_{\mathrm{app}}\tilde{w}_s}{\sqrt{12}}
\right)^2.
\tag{2.5.8}
$$

Eq. (2.5.7) represents a field dependent generalization of Fink's result, which is recovered in the limit where $q^2=0$. The expression for $q^2$ is related to the parallel critical field for sufficiently thin films. In effect, Blair's formalism treats the very narrow system as being within the overlapping sheath regions, with an enhanced upper critical field.

### General Mathematical Solution

In order to extract an analytic solution, it is convenient to rescale the variables. Blair made the following substitutions:

$$
\tilde{F}
=
\tilde{\mathcal{F}}
\sqrt{1 - \frac{q^2}{\tilde{m}_n \tilde{\alpha}_n}},
\qquad
\tilde{x}
=
\frac{\tilde{X}}
{\sqrt{1 - \frac{q^2}{\tilde{m}_n \tilde{\alpha}_n}}},
\qquad
\tilde{J}
=
\tilde{j}
\left(
1 - \frac{q^2}{\tilde{m}_n \tilde{\alpha}_n}
\right)^{3/2}.
\tag{2.5.9}
$$

Inserting these into Eq. (2.5.7) leads to the result

$$
\frac{1}{\tilde{m}_n}
\frac{\partial^2 \tilde{\mathcal{F}}}{\partial \tilde{X}^2}
+
\left[
\tilde{\alpha}_n
-
\tilde{\beta}_n \tilde{\mathcal{F}}^2
-
\frac{\tilde{m}_n
\left( \langle \tilde{J}_{\mathrm{app}} \rangle_y \right)^2}
{\tilde{\mathcal{F}}^4}
\right]
\tilde{\mathcal{F}}
= 0.
\tag{2.5.10}
$$

This is the same equation that Fink considered, allowing Blair to use the same machinery to obtain an exact analytic solution. In effect, the substitutions rescaled the variables to remove the field dependence, and provide the in-field observables $\tilde{F}$, $\tilde{x}$ and $\tilde{J}$.

To ensure that the critical current density vanishes at the upper critical field $B_{c2}^*$, $q^2$ is replaced in practice by the reduced field $B_{\mathrm{app}} / B_{c2}^*$ (which equals $\tilde{B}_{\mathrm{app}}$ in dimensionless GL units where $B_{c2}^* = \phi_0 / (2\pi\xi_s^2) = 1$). This identification follows from replacing the junction width with the inter-fluxon spacing $a_J \approx 1.84\,\phi_0 / (B_{\mathrm{app}} D)$ as the relevant transverse length scale in the definition of $q^2$, together with the Ginzburg–Landau relation $B_{c2}^* = \phi_0 / (2\pi\xi_s^2)$.

The critical current density of a narrow junction is then given by

$$
J_c = J_{\mathrm{DJ}} \left|\frac{\phi_0}{\pi B D (d_{\mathrm{n}} + 2\lambda)} \sin\left(\frac{\pi B D (d_{\mathrm{n}} + 2\lambda)}{\phi_0}\right)\right|,
$$

where $\phi_0 = h/q$ is the magnetic flux quantum and $\lambda$ is the penetration depth. The junction depairing current density is given by Blair (2022),

$$
\begin{aligned}
J_{\mathrm{DJ}} &= 4J_0\left(1 - \tilde{B}^*\right)^{\frac{3}{2}} \frac{1 - \sqrt{1-\tilde{s}\tilde{f}^2_{d/2}}}{\tilde{s}\tilde{v}}\exp\left(-\frac{\tilde{d}_{\mathrm{n}}}{\tilde{\xi}_{\mathrm{n}}}\right), \\[6pt]
\tilde{f}_{d/2}^2 &= \frac{\tilde{v}^2 + 1 - \sqrt{\tilde{v}^2(2-\tilde{s}) + 1}}{\tilde{v}^2 + \tilde{s}}, \\[6pt]
\tilde{v} &= \tilde{m}_{\mathrm{n}}\tilde{\xi}_{\mathrm{n}}\sqrt{1-\tilde{B}^*}, \qquad
\tilde{s} = \frac{\tilde{\beta}_{\mathrm{n}}\left(1-\tilde{B}^*\right)}{\tilde{\alpha}_{\mathrm{n}} - \tilde{B}^*\tilde{m}_{\mathrm{n}}^{-1}}, \\[6pt]
\tilde{\xi}_{\mathrm{n}} &= \sqrt{\frac{1}{-\tilde{m}_{\mathrm{n}}\tilde{\alpha}_{\mathrm{n}} + \tilde{B}^*}}, \qquad
\tilde{B}^* = \frac{B}{B_{c2}^*} \qquad\mathrm{and}\qquad
J_0 = \frac{\phi_0}{2\pi\mu_0\xi_s^3\kappa_s^2}.
\end{aligned}
$$

---

## 5. Transverse profiles $|\psi(x, 0)|^2$ and $B_z(x, 0)$ across the junction

To compare the simulation with the narrow and wide analytic models along the transverse line $y = 0$ (perpendicular to the barrier, centered at $x = 0$), we evaluate the spatial variation of both the magnetic field $B_z(x, 0)$ and the Cooper pair density $|\psi(x, 0)|^2$.

### Transverse magnetic field profile $B_z(x, 0)$

Inside the normal barrier ($|x| \le t/2$, where $t = \tilde{d}_{\mathrm{n}}$ is the bare barrier thickness), the Josephson supercurrent is uniform in $x$ and no London screening occurs, so the magnetic field is constant across the barrier thickness: $B_z(x, 0) = B_{\mathrm{center}}$. Inside the superconducting electrodes ($|x| > t/2$), Meissner screening causes the magnetic field to decay exponentially with the London penetration depth $\lambda_L = \kappa$:

$$
\boxed{B_z(x, 0) = \begin{cases} B_{\mathrm{center}}, & |x| \le t/2, \\[4pt] B_{\mathrm{center}}\,\exp\!\left(-\dfrac{|x| - t/2}{\kappa}\right), & |x| > t/2. \end{cases}}
$$

Note that integrating $B_z(x, 0)$ over $x \in (-\infty, \infty)$ gives $\int_{-\infty}^{\infty} B_z(x, 0)\,dx = (t + 2\kappa)B_{\mathrm{center}} = d\,B_{\mathrm{center}}$, directly recovering the effective magnetic thickness $d = t + 2\kappa$.

The barrier-center field $B_{\mathrm{center}} = B_z(0, 0)$ differs between the two regimes:
- **Narrow regime ($w \ll \lambda_J$):** Self-field screening along $y$ is negligible, so the field at the center of the barrier equals the applied field:

$$
B_{\mathrm{center}}^{\mathrm{narrow}} = B_{\mathrm{app}}.
$$

- **Wide regime ($w \gg \lambda_J$):** From the phase-gradient relation $B_z(0, y) = \frac{1}{d}\frac{\partial\Delta\phi_{\mathrm{wide}}}{\partial y}$ (in dimensionless units where $\phi_0 = 2\pi$) and the derivatives $\partial_y\Delta\phi_\pm = \pm\frac{2}{\lambda_J}\sin(\Delta\phi_\pm/2)$ with $B_0 = \frac{2}{d\lambda_J}$, the field at $y = 0$ is the superposition of the two exponentially decaying edge-field tails:

$$
B_{\mathrm{center}}^{\mathrm{wide}} = B_0\left(\sin\frac{\Delta\phi_+(0)}{2} - \sin\frac{\Delta\phi_-(0)}{2}\right).
$$

### Transverse order parameter profile $|\psi(x, 0)|^2$

Across an SNS junction carrying local supercurrent density $J_{\mathrm{local}} = J_x(y = 0)$ in local field $B_{\mathrm{local}} = |B_z(0, 0)|$, the 1D transverse order parameter magnitude $f(x) = |\psi(x, 0)|$ satisfies the stationary Ginzburg–Landau equation with effective field pair-breaking parameter $B_{\mathrm{local}}$:

$$
\frac{1}{m(x)}\frac{d^2 f}{dx^2} + \left[\alpha(x) - \frac{B_{\mathrm{local}}}{m(x)} - \beta(x)f^2 - \frac{m(x)J_{\mathrm{local}}^2}{f^4}\right]f = 0.
$$

Rescaling separately in each region removes the explicit material and field parameters:

1. **Superconducting electrodes ($|x| > t/2$, where $\alpha_s = \beta_s = m_s = 1$):** With $\tilde{x} = x\sqrt{1 - B_{\mathrm{local}}}$, $\tilde{f}_s = f / \sqrt{1 - B_{\mathrm{local}}}$, and $\tilde{j}_x = J_{\mathrm{local}}(1 - B_{\mathrm{local}})^{-3/2}$, the GL equation becomes

$$
\frac{d^2\tilde{f}_s}{d\tilde{x}^2} + \left(1 - \tilde{f}_s^2 - \frac{\tilde{j}_x^2}{\tilde{f}_s^4}\right)\tilde{f}_s = 0.
$$

Far inside the electrodes ($|x| \to \infty$), $d\tilde{f}_s/d\tilde{x} \to 0$ and $\tilde{f}_s^2 \to \tilde{f}_\infty^2$, where $\tilde{f}_\infty^2 \in [2/3, 1]$ is the largest root of the bulk current-carrying equilibrium condition $\tilde{f}_\infty^4(1 - \tilde{f}_\infty^2) = \tilde{j}_x^2$. Multiplying by $2\,d\tilde{f}_s/d\tilde{x}$ and integrating from $\tilde{x}$ to $\infty$ yields the first integral in the electrodes:

$$
\left(\frac{d\tilde{f}_s}{d\tilde{x}}\right)^2 = (\tilde{f}_\infty^2 - \tilde{f}_s^2) - \frac{1}{2}(\tilde{f}_\infty^4 - \tilde{f}_s^4) - \tilde{j}_x^2\left(\frac{1}{\tilde{f}_s^2} - \frac{1}{\tilde{f}_\infty^2}\right).
$$

2. **Normal barrier ($|x| \le t/2$, with parameters $\alpha_n < 0$, $\beta_n$, $m_n$):** With $\tilde{u} = x / \tilde{\xi}_n$, $\tilde{\xi}_n = \bigl(-m_n\alpha_n + B_{\mathrm{local}}\bigr)^{-1/2}$, $\tilde{f}_n = f\sqrt{\beta_n / (-\alpha_n + B_{\mathrm{local}}/m_n)}$, and $\tilde{j}_u = J_{\mathrm{local}}\beta_n\sqrt{m_n}\,\tilde{\xi}_n^3$, the GL equation becomes

$$
\frac{d^2\tilde{f}_n}{d\tilde{u}^2} - \left(1 + \tilde{f}_n^2 - \frac{\tilde{j}_u^2}{\tilde{f}_n^4}\right)\tilde{f}_n = 0,
$$

with symmetry condition $\left.d\tilde{f}_n/d\tilde{u}\right|_{\tilde{u}=0} = 0$ at the center of the barrier.

3. **S/N interface matching ($|x| = t/2$):** Continuity of the order parameter $f$ and its flux $m^{-1}df/dx$ across $|x| = t/2$ links the rescaled barrier and electrode variables via $\tilde{s}$ and $\tilde{v}$ (defined in Section 4):

$$
\tilde{f}_n^2(t/2) = (-\tilde{s})\,\tilde{f}_s^2(t/2),
\qquad
\left.\frac{d\tilde{f}_n}{d\tilde{u}}\right|_{t/2}^2 = (-\tilde{s})\,\tilde{v}^2\left.\frac{d\tilde{f}_s}{d\tilde{x}}\right|_{t/2}^2.
$$

Shooting from $\tilde{u} = 0$ with trial center value $\tilde{f}_n^2(0) = \tilde{f}_0^2$ to $\tilde{u} = t / (2\tilde{\xi}_n)$ fixes $\tilde{f}_0^2$ uniquely by satisfying this interface flux-matching condition, after which integrating outward in $\tilde{x}$ gives the complete profile $|\psi(x, 0)|^2 = f^2(x)$.

The two analytic models evaluate this 1D boundary-value solution at their respective mid-barrier conditions $(J_{\mathrm{local}}, B_{\mathrm{local}})$:
- **Narrow regime:** $J_{\mathrm{local}} = J_x^{\mathrm{narrow}}(0) = J_{\mathrm{DJ}}$ and $B_{\mathrm{local}} = B_{\mathrm{app}}$.
- **Wide regime:** $J_{\mathrm{local}} = J_x^{\mathrm{wide}}(0) = J_{\mathrm{DJ}}\sin\bigl(\Delta\phi_{\mathrm{wide}}(0)\bigr)$ and $B_{\mathrm{local}} = |B_{\mathrm{center}}^{\mathrm{wide}}|$.

---

## Summary of regimes

| Regime | Condition | $J_c$ vs. applied field | $J_x(y)$ profile | $A_x(y)$ profile | $B_z(x, 0)$ center field $B_{\mathrm{center}}$ | $|\psi(x, 0)|^2$ local inputs $(J_{\mathrm{local}}, B_{\mathrm{local}})$ |
|---|---|---|---|---|---|---|
| Narrow | $w \ll \lambda_J$ | Fraunhofer: $J_{\mathrm{DJ}}\left\lvert\dfrac{\sin(\pi w d B/\phi_0)}{\pi w d B/\phi_0}\right\rvert$ | $J_{\mathrm{DJ}}\cos(ky)$ | $-B_{\mathrm{app}}\,y$ | $B_{\mathrm{app}}$ | $\bigl(J_{\mathrm{DJ}},\, B_{\mathrm{app}}\bigr)$ |
| Wide | $w \gg \lambda_J$ | Linear decrease: $\dfrac{2}{w}\left(2\lambda_J J_{\mathrm{DJ}} - \dfrac{\lvert B_{\mathrm{app}}\rvert}{\mu_0}\right)$ | $J_{\mathrm{DJ}}\sin\bigl(\Delta\phi_+(y) + \Delta\phi_-(y)\bigr)$ | $-\dfrac{\phi_0}{2\pi d}\bigl(\Delta\phi_{\mathrm{wide}}(y) - \Delta\phi_{\mathrm{wide}}(0)\bigr)$ | $B_0\left(\sin\dfrac{\Delta\phi_+(0)}{2} - \sin\dfrac{\Delta\phi_-(0)}{2}\right)$ | $\bigl(J_{\mathrm{DJ}}\sin\Delta\phi_{\mathrm{wide}}(0),\, \lvert B_{\mathrm{center}}^{\mathrm{wide}}\rvert\bigr)$ |
