# (1+1)-D Cooper-Frye
Code for (1+1)D Cooper-Frye (not only for SPH)

**_dndeta.cpp:_**

This code receives an emitted Freeze-out (1+1)D hypersurface and calculates the pseudo-rapidity distribution. It implements the gross distribution calculation and a smeared calculation. The input file must be organized into columns as follows (space between the columns):

$\tau$    $\eta_s$    $u_\tau$    $u_{\eta_s}$    $\mathrm{d}s_\tau$    $\mathrm{d}s_\eta$    $s$    $T$    $p$    $\epsilon$

where $s$ is entropy density, $T$ is temperature, $p$ is pressure and $\epsilon$ is energy density. $\mathrm{d}s_\mu$ is the covariant hypersurface element. The output files is _dndeta.dat_ (gross distribution calculation) and _sdndeta.dat_ (smeared calculation). The unit of Freeze-out temperature must be GeV.

Since 

$$ p^0 \dfrac{\mathrm{d}^3N}{\mathrm{d}^3p} = \dfrac{1}{m_T}\dfrac{\mathrm{d}^3N}{\mathrm{d}m_T\mathrm{d}\phi \mathrm{d}y}$$

[where $p$ is parametrized as $p=(m_T \cosh y, p_T \cos \phi, p_T \sin \phi, m_T \sinh y)$] and 

$$ \mathrm{d} y = \dfrac{|\vec p|}{E} \mathrm{d}\eta $$ 

we get

$$ p^0 \dfrac{\mathrm{d}^3N}{\mathrm{d}^3p} = \dfrac{|\vec p|}{p^0 m_T}\dfrac{\mathrm{d}^3N}{\mathrm{d}m_T\mathrm{d}\phi \mathrm{d}\eta} $$

and this leads for

$$ \frac{\mathrm{d} N}{\eta} = \int\limits_0^{2\pi} \mathrm{d} \phi \int\limits_m^{\inf} \mathrm{d} m_T \frac{E m_ T}{|\vec p|} \left(p^0 \dfrac{\mathrm{d}^3N}{\mathrm{d}^3p}\right)$$.


