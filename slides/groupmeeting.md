---
marp: true
author: Casper Pranger
paginate: true
math: katex
footer: SIO/LMU CompEQteam group meeting, March 12, 2024
header: Pranger -- Numerical time integration of multi-state rock rheology
style: |
    section {
        font-size: 25px;
        justify-content: flex-start;
    }


---
$$
    \\[5em]
$$

# <!-- fit --> Numerical time integration of multi-state  rock rheology <br/> $\hspace{11em} \hookrightarrow$ multi-rate
_Casper Pranger, Dave May, Alice Gabriel_


---
# Motivation

- Damage-Breakage rheology of Lyakhovsky, Ben-Zion, _et al_.
    - **slow:** damage
    - **fast:** breakage
    - **intermittent:** elastic wave physics
    - includes bulk rate-and-state (Pranger et al., 2022, JGR:SE) <br>  as critically damaged end-member
- Grain-size assisted thermal runaway of Thielmann, _et al_.
    - **slow:** grain size
    - **fast:** temperature
    - **intermittent:** elastic wave physics
- Others: multi-phase, thermo-chemical, etc.


---
# Damage-Breakage rheology (DBR)
Lyakhovksy, Ben-Zion, _et al_.

- Cauchy's linear momentum balance law:
$$
    \partial_t v = \nabla \cdot s + g,
$$
-   - $\partial_t$: (partial) derivative with respect to time $t$
    - $v = v(t,x)$: mass flux at time $t$ at material coordinate $x = x(t)$
    - $s = s(t,x)$: Cauchy's stress tensor
    - $g = g(t,x)$: (gravitational) body forces;  $g = 0$.


---
# Damage-Breakage rheology (DBR)

- Isotropic nonlinear constitutive equation with elastic and structural stresses:
$$
	\begin{align*}
        s(e, \alpha) &= s_\mathrm{e}(e, \alpha) + s_\mathrm{s}(\nabla \alpha) \\[.5em]
        s_\mathrm{e}(e, \alpha) &= \lambda(x,\alpha, \jmath e) \,\delta \;\mathrm{tr}\; e + 2 \mu(x,\alpha, \jmath e) e \\[.5em]
        s_\mathrm{s}(\nabla \alpha) &= - \vartheta(x) (\nabla \alpha) \otimes (\nabla \alpha)
    \end{align*}
$$
-   - $e = e(t,x)$: elastic strain field
    - $\alpha = \alpha(t,x)$: damage field
    - $\lambda$ and $\mu$: Lamé parameters of isotropic elastic solid; $\vartheta$: modulus of structural stress
    - $\jmath e$: scalar invariants $(J_1, J_2, J_3)$ of $e$, indicating isotropy
    - $\delta$: three-dimensional identity tensor


---
# Damage-Breakage rheology (DBR)

- Additive elastic-anelastic strain decomposition (total = elastic + anelastic), combined with Koiter's rule:
$$
	\partial_t e = \nabla^\mathrm{s} (r v) - \partial_t \bar{e} = \nabla^\mathrm{s} (r v) - \gamma(\beta) \partial_s G(\jmath s(e,\alpha)),
$$
-   - $\bar{e} = \bar{e}(t,x)$: _anelastic_ strain field
    - $r v = rv(t,x)$: material velocity (density $r$ $\times$ mass flux $v$)
    - $\gamma(\beta)$: breakage rheology, with $\beta = \beta(t, x)$ the breakage field
    - $\partial_s G(\jmath s(e,\alpha))$: 'Schmidt tensor', with plastic potential $G$; $\jmath s(e,\alpha)$ the  scalar invariants of $s$

---
# Damage-Breakage rheology (DBR)

- Evolution equations for damage field $\alpha$ and breakage field $\beta$:
$$
	\begin{align*}
        \partial_t \alpha &= \nabla \cdot \left[ D_\alpha(\alpha)\nabla \alpha \right] + R_\alpha(\jmath e, \alpha, \beta) \\[.7em]
        \partial_t \beta &= \nabla \cdot \left[ D_\beta(\beta) \nabla \beta \right] + R_\beta(\jmath e, \alpha, \beta)
    \end{align*}
$$
-   - $D_\alpha(\alpha)$, $D_\beta(\beta)$: nonlinear diffusivities
    - $R_\alpha(\jmath e, \alpha, \beta)$, $R_\beta(\jmath e, \alpha, \beta)$: reaction, interaction, forcing terms


---
# Damage-Breakage rheology (DBR)

- Summary:
$$
	\begin{align*}
        \partial_t v &= \nabla \cdot s_\mathrm{e}(e, \alpha) + R_v(\alpha)\\[.7em]
        \partial_t e &= \nabla^\mathrm{s} (r v) + R_e(\jmath e, \alpha, \beta)  \\[.3em] \hline\\[-.8em]
        \partial_t \alpha &= \nabla \cdot \left[ D_\alpha(\alpha)\nabla \alpha \right] + R_\alpha(\jmath e, \alpha, \beta)  \\[.3em] \hline\\[-.8em]
        \partial_t \beta &= \nabla \cdot \left[ D_\beta(\beta) \nabla \beta \right] + R_\beta(\jmath e, \alpha, \beta)
    \end{align*}
$$
-   - $R_v(\alpha)$: contribution to momentum balance due to the divergence of structural stresses
    - $R_e(\jmath e, \alpha, \beta)$: _negative_ contribution to elastic strain accumulation by anelastic processes


---
# Example: Damage-Breakage rheology (DBR)

- Summary in quasi-Jacobian form: $\quad \partial_t U = \left(J_\mathrm{PDE} + J_\mathrm{interact} \right) U, \quad U = [v, e, \alpha, \beta]^\mathrm{T}\\[1em]$

$$
    \begin{align*}
	    J_\mathrm{PDE} &= \left[\begin{array}{cc|c|c}
            0 & \nabla \cdot s_\mathrm{e}(\bullet, \alpha) & 0 & 0 \\
            \nabla^\mathrm{s} (r\, \bullet) & 0 & 0 & 0 \\[.3em] \hline\\[-.8em]
            0 & 0 & \nabla \cdot \left[ D_\alpha(\bullet)\nabla \bullet \right] & 0 \\[.3em] \hline\\[-.8em]
            0 & 0 & 0 & \nabla \cdot \left[ D_\beta(\bullet)\nabla \bullet \right]
        \end{array}\right] \\[3em]

        J_\mathrm{interact} &=  \left[\begin{array}{c|ccc}
            0 & 0 & \nabla \cdot s_\mathrm{e}(e, \bullet) + R_v(\bullet) & 0 \\[.3em] \hline\\[-.8em]
            0 & R_e(\jmath\, \bullet, \alpha, \beta) & R_e(\jmath e, \bullet, \beta) & R_e(\jmath e, \alpha, \bullet) \\
            0 & R_\alpha(\jmath\, \bullet, \alpha, \beta) & R_\alpha(\jmath e, \bullet, \beta) & R_\alpha(\jmath e, \alpha, \bullet) \\
            0 & R_\beta(\jmath\, \bullet, \alpha, \beta) & R_\beta(\jmath e, \bullet, \beta) & R_\beta(\jmath e, \alpha, \bullet)
        \end{array}\right]
    \end{align*}
$$

---
# Damage-Brakage rheology (DBR)

- The subsystem $\partial_t [v, e]^\mathrm{T}$ is an inhomogeneous, non-linear elastic wave equation in coupled first-order form
- The subsystems $\partial_t \alpha$ and $\partial_t \beta$ are non-linear reaction-diffusion equations
- "Stick-Slip" episodes of seismogenesis:
    - **Stick**: interaction terms $R_{v,e,\alpha,\beta} \ll 1$; $\quad\partial_t v = \nabla \cdot s_\mathrm{e}(e, \alpha) \approx 0$
    - **Slip**: interaction terms $R_{v,e,\alpha,\beta} \gg 0$.
- Explicit time integration subject to "CFL condition":
    - Hyperbolic (wave equation): $\Delta t \propto C^{-1} \Delta x$
    - Parabolic (diffusion equation): $\Delta t \propto D^{-1} \Delta x^2$
- Implicit time integration requires nonlinear (iterative) solvers
    - places constraints on the structure of $J_\mathrm{PDE}$,  $J_\mathrm{interact}$
    - no stability constraints on $\Delta t$.


---
# Additive IMEX Runge-Kutta methods

- We need to solve the ODE $\quad \partial_t U = F(U; t) = F^\mathrm{im}(U; t) + F^\mathrm{ex}(U; t),$
-   - $F^\mathrm{im}(U; t)$ contains the terms to be integrated implicitly (PDE terms),
    - $F^\mathrm{ex}(U; t)$ contains the terms to be integrated explicitly (interaction terms).

- We use the three-stage additive Runge-Kutta (ARK) scheme:
$$
    \begin{align*}
	    t_s = t^n + &\Delta t\,c_s, \qquad s = 1,2,3 \\
        W_s = U^n + &\Delta t \sum_{j = 1}^{s} a^\mathrm{im}_{sj} F^\mathrm{im}(W_j; t_j) + \Delta t \sum_{j = 1}^{s-1} a^\mathrm{ex}_{sj} F^\mathrm{ex}(W_j; t_j), \qquad s = 1,2,3 \\
        U^{n+1} = U^n + &\Delta t \sum_{j = 1}^{3} b_{j}\left(F^\mathrm{im}(W_j; t_j) + F^\mathrm{ex}(W_j; t_j)\right)
    \end{align*}
$$


---
# Additive IMEX Runge-Kutta methods

- We need to solve the ODE $\quad \partial_t U = F(U; t) = F^\mathrm{im}(U; t) + F^\mathrm{ex}(U; t),$
-   - $F^\mathrm{im}(U; t)$ contains the terms to be integrated implicitly (PDE terms),
    - $F^\mathrm{ex}(U; t)$ contains the terms to be integrated explicitly (interaction terms).

- We use the three-stage additive Runge-Kutta (ARK) scheme:
$$
    \begin{align*}
	    \bar{T} = t^n \bar{\mathrm{1}} + &\Delta t\, \bar{C}, \\
        \bar{W} = U^n \bar{\mathrm{1}} + &\Delta t \left[\bar{\mathbf{A}}^\mathrm{im} \bar{F}^\mathrm{im}(\bar{W}; \bar{T}) + \bar{\mathbf{A}}^\mathrm{ex} \bar{F}^\mathrm{ex}(\bar{W}; \bar{T})\right] \\
        U^{n+1} = U^n + &\Delta t\, \bar{B} \left[\bar{F}^\mathrm{im}(\bar{W}; \bar{T}) + \bar{F}^\mathrm{ex}(\bar{W}; \bar{T}) \right]
    \end{align*}
$$
- with $\bar{\mathrm{1}} = [1,1,1]^\mathrm{T}$, $\bar{C} = [c_1, c_2, c_3]^\mathrm{T}$, $\bar{B} = [b_1, b_2, b_3]$ and likewise $\bar{\mathbf{A}}^\mathrm{im}$ and $\bar{\mathbf{A}}^\mathrm{ex}$ matrices of coefficients $a^\mathrm{im}$ and $a^\mathrm{ex}$.

---
# Additive IMEX Runge-Kutta methods

- We use the three-stage additive Runge-Kutta (ARK) scheme:
$$
    \begin{align*}
	    \bar{T} = t^n \bar{\mathrm{1}} + &\Delta t\, \bar{C}, \\
        \bar{W} = U^n \bar{\mathrm{1}} + &\Delta t \left[\bar{\mathbf{A}}^\mathrm{im} \bar{F}^\mathrm{im}(\bar{W}; \bar{T}) + \bar{\mathbf{A}}^\mathrm{ex} \bar{F}^\mathrm{ex}(\bar{W}; \bar{T})\right] \\
        U^{n+1} = U^n + &\Delta t\, \bar{B} \left[\bar{F}^\mathrm{im}(\bar{W}; \bar{T}) + \bar{F}^\mathrm{ex}(\bar{W}; \bar{T}) \right]
    \end{align*}
$$
- with $\bar{\mathrm{1}} = [1,1,1]^\mathrm{T}$, $\bar{C} = [c_1, c_2, c_3]^\mathrm{T}$, $\bar{B} = [b_1, b_2, b_3]$ and likewise $\bar{\mathbf{A}}^\mathrm{im}$ and $\bar{\mathbf{A}}^\mathrm{ex}$ matrices of coefficients $a^\mathrm{im}$ and $a^\mathrm{ex}$.
- Additionally,
    - $\bar{B} \cdot \bar{1} = 1$, $\quad\bar{C} = \bar{\mathbf{A}}^\mathrm{im}\bar{1} = \bar{\mathbf{A}}^\mathrm{ex} \bar{1}$,
    - $\bar{\mathbf{A}}^\mathrm{im}$ is lower-triangular (diagonally implicit),
    - $\bar{\mathbf{A}}^\mathrm{ex}$ is _strictly_ lower-triangular (explicit).


---
# Additive IMEX Runge-Kutta methods

- Populate $\bar{\mathbf{A}}^\mathrm{im}$ with coefficients from the TR-BDF2 time integration scheme [e.g. Giraldo 2013]
    1) an explicit first stage $W_1 = U^{n}$ at $T_1 = t^n$,
    2) a trapezoidal second stage at $t_2 = t^n + c_2 \Delta t$ given by </br> $W_2 = U^{n} + (c_2/2) F^\mathrm{im}(W_1; T_1) + (c_2/2) F^\mathrm{im}(W_2; T_2))$,
    3) a third stage at $T_3 = t^n + \Delta t$ given by the 2nd-order backwards difference formula (BDF2)
    4) a finishing rule compatible with the implicit third stage; i.e.</br> $U^{n+1} = W_3 \iff \bar{B}_2 = \bar{\mathbf{A}}^\mathrm{im} [0, 0, 1]^\mathrm{T}$.
- Keep open the values of $\bar{B}$ in the hope of generating different schemes that re-use $\bar{W}$.

---
# Additive IMEX Runge-Kutta methods

- The resulting scheme is summarized as

$$
    \begin{align*}
        \bar{\mathbf{A}}^\mathrm{im} &= \begin{bmatrix}
            0 & 0 & 0 \\
            \frac{1}{2}c_2 & \frac{1}{2}c_2 & 0  \\
            \frac{1}{2}(2-c_2)^{-1} & \frac{1}{2}(2-c_2)^{-1} & (1-c_2)(2-c_2)^{-1}
        \end{bmatrix} \\[2em]
        \bar{B}_2 &= \left[\begin{matrix}
            \frac{1}{2}(2-c_2)^{-1} & \frac{1}{2}(2-c_2)^{-1} & (1-c_2)(2-c_2)^{-1}
        \end{matrix}\right] \\[.5em]
        \bar{\mathbf{A}}^\mathrm{ex} &= \begin{bmatrix}
            0 & 0 & 0 \\
            c_2 & 0 & 0 \\
            1-a^{ex}_{32} & a^{ex}_{32} & 0
        \end{bmatrix}, \quad \bar{C} = \begin{bmatrix}
            0 \\
            c_2  \\
            1
        \end{bmatrix} \\[2em]
        \bar{B}_3 &= [b_1, b_2, 1 - b_1 - b_2]
    \end{align*}
$$

- Four tunable parameters: $b_1$, $b_2$, $c_2$ and $a^{ex}_{32}$.


---
# Additive IMEX Runge-Kutta methods

- IMEX version of Dahlquist's problem [Dahlquist, 1963]: $\quad\Delta t\, \partial_t U = \zeta^\mathrm{im} U + \zeta^\mathrm{ex} U$
- Analytical solution over one step: $\tilde{U}^{n+1} = \mathrm{exp}(\zeta^\mathrm{im} + \zeta^\mathrm{ex}) U^n$
- Interpret $\zeta^\mathrm{im}$ and $\zeta^\mathrm{ex}$ as 'eigenvalues' of $\Delta t\, F^\mathrm{im}(\bullet; t)$ and $\Delta t\, F^\mathrm{ex}(\bullet; t)$.
- The Additive Runge-Kutta scheme becomes:
$$
	U^{n+1} = P(\zeta^\mathrm{im}, \zeta^\mathrm{ex}) U^n,
$$

$$
	P(\zeta^\mathrm{im}, \zeta^\mathrm{ex}) = 1 + (\zeta^\mathrm{im} + \zeta^\mathrm{ex})\, \bar{B} \left[\bar{\mathbf{I}} - \zeta^\mathrm{im} \bar{\mathbf{A}}^\mathrm{im} - \zeta^\mathrm{ex} \bar{\mathbf{A}}^\mathrm{ex} \right]^{-1} \bar{\mathrm{1}}.
$$
- The consistency of the method is determined by the error
$$
    \tilde{\epsilon}^{n+1} = \tilde{U}^{n+1} - U^{n+1} = \left[\mathrm{exp}(\zeta^\mathrm{im} + \zeta^\mathrm{ex}) - P(\zeta^\mathrm{im}, \zeta^\mathrm{ex})\right] U^n
$$

- Coefficients $b_1$, $b_2$, $c_2$ and/or $a^{ex}_{32}$ are obtained by Taylor series expansion at $\zeta^\mathrm{im}, \zeta^\mathrm{ex} = 0$ 


---
# Additive IMEX Runge-Kutta methods

- The Taylor series expansion looks like this:
$$
    \begin{align*}
        \mathrm{exp}(\zeta^\mathrm{im} + \zeta^\mathrm{ex}) - P(\zeta^\mathrm{im}, \zeta^\mathrm{ex}) &= \mathcal{O}(\zeta_\mathrm{im} + \zeta_\mathrm{ex})^3\\
            &= \mathcal{O}(\zeta_\mathrm{im}^3) + \mathcal{O}(\zeta_\mathrm{im}^2 \zeta_\mathrm{ex}) + \mathcal{O}(\zeta_\mathrm{im} \zeta_\mathrm{ex}^2) + \mathcal{O}(\zeta_\mathrm{ex}^3)
    \end{align*}
$$

- Second-order consistency (accuracy) is given by the TR-BDF2 scheme with $\bar{B} = \bar{B}_2$

- Big-O notation implies proportionality; coefficients expressed in terms of tunable ARK parameters; eliminated by
    - $b_1 = \frac{1}{2} + b_2 (c_2 - 1)$
    - $b_2 = \frac{1}{6}(c_2 - c_2^2)^{-1}$
    - $a_{32}^\mathrm{ex} = (c_2 - 1)(3 c_2^2 - 2 c_2)^{-1}$
    - $\Rightarrow \mathrm{exp}(\zeta^\mathrm{im} + \zeta^\mathrm{ex}) - P(\zeta^\mathrm{im}, \zeta^\mathrm{ex}) = \mathcal{O}(\zeta_\mathrm{im} + \zeta_\mathrm{ex})^4$


---
# Additive IMEX Runge-Kutta methods

- Choice: $c_2 = 1 - \frac{1}{3}\sqrt{3} \approx 0.42265$
$$
    \begin{align*}
        \bar{B}_2 &= \left[\begin{matrix}
            \frac{1}{2}(2-c_2)^{-1} &
            \frac{1}{2}(2-c_2)^{-1} &
            (1-c_2)(2-c_2)^{-1}
        \end{matrix}\right] \\[.5em]
        \bar{B}_3 &= \left[\begin{matrix}
            \frac{1}{2}-\frac{1}{6}c_2^{-1} &
            \frac{1}{6}(c_2 - c_2^2)^{-1} &
            \frac{1}{2} + \frac{1}{6}(c_2-1)^{-1}
        \end{matrix}\right]
    \end{align*}
$$

- _Embedded_ ARK method:
$$
    \begin{align*}
	    \bar{T} = t^n \bar{\mathrm{1}} + &\Delta t\, \bar{C}, \\
        \bar{W} = U^n \bar{\mathrm{1}} + &\Delta t \left[\bar{\mathbf{A}}^\mathrm{im} \bar{F}^\mathrm{im}(\bar{W}; \bar{T}) + \bar{\mathbf{A}}^\mathrm{ex} \bar{F}^\mathrm{ex}(\bar{W}; \bar{T})\right] \\
        U^{n+1} = U^n + &\Delta t\, \bar{B}_2 \left[\bar{F}^\mathrm{im}(\bar{W}; \bar{T}) + \bar{F}^\mathrm{ex}(\bar{W}; \bar{T}) \right] \\
        E^{n+1} = &\Delta t\, (\bar{B}_3 - \bar{B}_2) \left[\bar{F}^\mathrm{im}(\bar{W}; \bar{T}) + \bar{F}^\mathrm{ex}(\bar{W}; \bar{T}) \right]
    \end{align*}
$$

- Dimensionless error: $\epsilon_{n+1} = \lVert E^{n+1}\rVert_2 \left(\tau_r \lVert U^{n+1} \rVert_2 + \tau_a\right)^{-1}$, $\quad \Delta t^{n+1} = \Delta t^n \epsilon_{n+1}^{-1/3}$

---
# Additive IMEX Runge-Kutta methods

- _Stability_ region $S = \{ \zeta^\mathrm{im}, \zeta^\mathrm{ex} \in \mathbb{C} : \lvert P(\zeta^\mathrm{im}, \zeta^\mathrm{ex}) \rvert < 1 \}$

$\qquad\qquad\qquad$ ![h:380px](https://media.githubusercontent.com/media/cpranger/DamageBreakage/main/paper/figures/stability.png)

- **L**-stable at $\zeta^\mathrm{ex} = 0$ when $\lim\limits_{|\zeta^\mathrm{im}| \to \infty} \lvert P(\zeta^\mathrm{im}, 0) \rvert  = 0$ (TR-BDF2 is L-stable)


---
# Example: Damage-Breakage rheology (DBR)

- Return to 'CFL' conditions: $\Delta t$ constrained by the 'spectral radius' $\rho$ of $J^\mathrm{ex}$

$$
    \begin{align*}
	    J_\mathrm{PDE} &= \left[\begin{array}{cc|c|c}
            0 & \nabla \cdot s_\mathrm{e}(\bullet, \alpha) & 0 & 0 \\
            \nabla^\mathrm{s} (r\, \bullet) & 0 & 0 & 0 \\[.3em] \hline\\[-.8em]
            0 & 0 & \nabla \cdot \left[ D_\alpha(\bullet)\nabla \bullet \right] & 0 \\[.3em] \hline\\[-.8em]
            0 & 0 & 0 & \nabla \cdot \left[ D_\beta(\bullet)\nabla \bullet \right]
        \end{array}\right] \\[3em]

        J_\mathrm{interact} &=  \left[\begin{array}{c|ccc}
            0 & 0 & \nabla \cdot s_\mathrm{e}(e, \bullet) + R_v(\bullet) & 0 \\[.3em] \hline\\[-.8em]
            0 & R_e(\jmath\, \bullet, \alpha, \beta) & R_e(\jmath e, \bullet, \beta) & R_e(\jmath e, \alpha, \bullet) \\
            0 & R_\alpha(\jmath\, \bullet, \alpha, \beta) & R_\alpha(\jmath e, \bullet, \beta) & R_\alpha(\jmath e, \alpha, \bullet) \\
            0 & R_\beta(\jmath\, \bullet, \alpha, \beta) & R_\beta(\jmath e, \bullet, \beta) & R_\beta(\jmath e, \alpha, \bullet)
        \end{array}\right]
    \end{align*}
$$


---
# Example: Damage-Breakage rheology (DBR)

- Return to 'CFL' conditions: $\Delta t$ constrained by the 'spectral radius' $\rho$ of $J^\mathrm{ex}$

$$
    \begin{align*}
	    J_\mathrm{interact} &=  \left[\begin{array}{c|ccc}
            0 & 0 & \nabla \cdot s_\mathrm{e}(e, \bullet) + R_v(\bullet) & 0 \\[.3em] \hline\\[-.8em]
            0 & R_e(\jmath\, \bullet, \alpha, \beta) & R_e(\jmath e, \bullet, \beta) & R_e(\jmath e, \alpha, \bullet) \\
            0 & R_\alpha(\jmath\, \bullet, \alpha, \beta) & R_\alpha(\jmath e, \bullet, \beta) & R_\alpha(\jmath e, \alpha, \bullet) \\
            0 & R_\beta(\jmath\, \bullet, \alpha, \beta) & R_\beta(\jmath e, \bullet, \beta) & R_\beta(\jmath e, \alpha, \bullet)
        \end{array}\right]
    \end{align*}
$$

- By the Schur determinant theorem, $\rho(J_\mathrm{interact})$ is independent of the entry $\nabla \cdot s_\mathrm{e}(e, \bullet) + R_v(\bullet)$!

- Integrate PDE terms implicitly (L-stable, independent of grid)
- Integrate interaction terms explicitly!
- https://github.com/cpranger/DamageBreakage
- THE END :-)