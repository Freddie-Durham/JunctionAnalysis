# Important Equations: Critical Current of a Josephson Junction in a Parallel Field

## 1. Setup and the sine-Gordon equation

Take a Josephson junction with current flowing in the $x$ direction. The barrier occupies $-d/2 < x < d/2$, and the phase varies along the barrier in the $y$ direction. A magnetic field is applied in the $z$ direction:

$$
\mathbf{B} = B_z(y)\,\hat{\mathbf{z}}, \qquad
B_z = (\nabla\times\mathbf{A})_z = \frac{\partial A_y}{\partial x} - \frac{\partial A_x}{\partial y}.
$$

Here $q$ is the charge of the superconducting carriers (Cooper pairs, $q = 2e$), and the flux quantum is $\phi_0 = h/q$.

**Supercurrent (Ginzburg–Landau).** We use the standard sign convention:

$$
\mathbf{J}_s = \frac{q|\psi|^2}{m}\left(\hbar\nabla\phi - q\mathbf{A}\right).
$$

**Josephson relation and gauge-invariant phase difference.** Current flows toward increasing phase, consistent with $\mathbf{J}_s \propto \nabla\phi$:

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

### Critical current of a finite wide junction

Now take a junction of width $w$ with the Meissner state in the interior. Define $J$ as the average current density, i.e. the total current per unit length (along $z$) divided by $w$. Integrating Ampère's law across the width:

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

---

## Summary of regimes

| Regime | Condition | $J_c$ vs. applied field |
|---|---|---|
| Narrow | $w \ll \lambda_J$ | Fraunhofer: $J_{\mathrm{DJ}}\left\lvert\dfrac{\sin(\pi w d B/\phi_0)}{\pi w d B/\phi_0}\right\rvert$ |
| Wide | $w \gg \lambda_J$ | Linear decrease: $\dfrac{2}{w}\left(2\lambda_J J_{\mathrm{DJ}} - \dfrac{\lvert B_{\mathrm{app}}\rvert}{\mu_0}\right)$, with $J_c \propto 1/w$ |
