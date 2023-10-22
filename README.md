# -1-1-D-Cooper-Frye-SPH
Code for (1+1)D Cooper-Frye (not only for SPH)

**_dndeta.cpp:_**

This code receives an emitted Freeze-out (1+1)D hypersurface and calculates the pseudo-rapidity distribution. It implements the gross distribution calculation and a smeared calculation. The input file must be organized into columns as follows (space between the columns):

$\tau$    $\eta_s$    $u_\tau$    $u_{\eta_s}$    $\mathrm{d}s_\tau$    $\mathrm{d}s_\eta$    $s$    $T$    $p$    $\epsilon$

where $s$ is entropy density, $T$ is temperature, $p$ is pressure and $\epsilon$ is energy density. $\mathrm{d}s_\mu$ is the covariant hypersurface element. The output files is _dndeta.dat_ (gross distribution calculation) and _sdndeta.dat_ (smeared calculation). 
