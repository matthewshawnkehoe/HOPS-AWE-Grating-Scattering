"""Global switches.

ALPHA_FIX (default True): correct three places where the MATLAB code (src.zip / GitHub) drops
terms that only matter for oblique incidence (alpha != 0).  With alpha = 0 every switch is a
no-op, so all alpha = 0 results (paper Figs. 2-10, the MATLAB reference data) are unchanged.

 1. T_dno is called with p instead of alpha_p = alpha + p in the field solvers.  The first
    frequency coefficient of gamma_p(delta) is (k^2 - alpha*alpha_p)/gamma_p, so MATLAB is off by
    alpha^2/gamma_p (the manufactured solution, which uses gamma_exp, is correct).
 2. setup_zeta_psi_n_m keeps only alpha_bar of alpha(delta) = alpha_bar (1 + delta) in the
    f_x * i*alpha term of psi (the loop variable ell = m after the loop), dropping the delta part.
 3. The HOPS unknowns carry the Bloch phase e^{i alpha x} removed, and the DNOs G, J are computed
    with d_x only, so the true normal derivatives are  G + i alpha g_x U  and  J - i alpha g_x W.
    The interface equation therefore needs  R -= i alpha(delta) g_x (U - tau^2 W).  MATLAB's
    two_layer_solve.m tries this with shifted indices (U_{n-1,m-1} instead of U_{n-1,m}) and
    two_layer_solve_fast.m has it commented out.
 4. The TFE right-hand side keeps 2 i alpha S d_x' u but drops the part of 2 i alpha d_x that the
    change of variables z' = a (z - g)/(a - g) sends to d_z':  d_x = d_x' + (dz'/dx) d_z', and
    (a-g)^2/a^2 * dz'/dx = eps A1_xz + eps^2 A2_xz (same coefficients as the A_xz terms already in
    the recursion; likewise for the lower layer).  Missing terms: 2 i alpha(delta) A_xz d_z' u.
    Verified by the manufactured solution with alpha != 0 (tests::test_mms_oblique_incidence).
Set ALPHA_FIX = False to reproduce the MATLAB code exactly.
"""
ALPHA_FIX = True
