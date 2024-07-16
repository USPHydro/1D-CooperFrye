# (1+1)-D Cooper-Frye
Code for (1+1)D Cooper-Frye (not only for SPH)

**_dndeta.cpp:_**

This code receives an emitted Freeze-out (1+1)D hypersurface and calculates the pseudo-rapidity distribution. It implements the gross distribution calculation and a smeared calculation. The input file must be organized into columns as follows (space between the columns):

$\tau$    $\eta_s$    $u_\tau$    $u_{\eta_s}$    $\mathrm{d}s_\tau$    $\mathrm{d}s_\eta$    $s$    $T$    $p$    $\epsilon$

where $s$ is entropy density, $T$ is temperature, $p$ is pressure and $\epsilon$ is energy density. $\mathrm{d}s_\mu$ is the covariant hypersurface element. The output files is _dndeta.dat_ (gross distribution calculation) and _sdndeta.dat_ (smeared calculation). The unit of Freeze-out temperature must be GeV.


**_Calculation Method:_**

Since the invariant one-particle distribution is

$$ \mathcal N_1(p) = p^0 \dfrac{\mathrm{d}^3N}{\mathrm{d}^3p} = \dfrac{1}{m_T}\dfrac{\mathrm{d}^3N}{\mathrm{d}m_T\mathrm{d}\phi \mathrm{d}y}$$

[where $p$ is parameterized as $p=(m_T \cosh y, p_T \cos \phi, p_T \sin \phi, m_T \sinh y)$] and 

$$ \mathrm{d} y = \dfrac{|\vec p|}{p^0} \mathrm{d}\eta $$ 

we get

$$ \mathcal N_1(p) = \dfrac{p^0}{|\vec p| m_T}\dfrac{\mathrm{d}^3N}{\mathrm{d}m_T\mathrm{d}\phi \mathrm{d}\eta} $$

and this leads for

$$\dfrac{\mathrm{d} N}{\mathrm{d}\eta} = \int\limits_0^{2\pi} \mathrm{d} \phi \int\limits_m^{\infty} \mathrm{d} m_T \frac{|\vec p| m_ T}{p^0} \mathcal N_1(p). $$ 

The Cooper-Frye description for the invariant one-particle distribution is

$$ \mathcal N_1(p) = \int \mathrm{d} \sigma \cdot p f(u \cdot p) $$

and we discretize this in hypersurface elements, as

$$ \mathcal N_1(p) = \sum_i \mathrm{d} \sigma_i \cdot p f(u_i \cdot p). $$

At this point, the SPH method provides at Freeze-out $\mathrm{d}\sigma_i, u_i$ and termodynamical quantities such that temperature, energy density, etc. Then we must only construt $p$ to calculate \mathcal N_1(p) and therefore ${\mathrm{d} N}/{\mathrm{d}\eta}$


