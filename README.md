This repo contains the implementation of similarity-restoring randomized-Arnoldi process proposed in [1], and scripts to reproduce the numerical experiments.

To produce numerical results in Section 5.1 about eigenvalue problems, run code in the eigenvalue folder. 

  main_ks_performance.m for Figure 5.1

  main_dMax.m for Figure 5.2 and Figure 5.3

To produce numerical results in Section 5.2 about matrix functions, run code in the matfun folder.
  
  test_arnoldi_matfun_conv.m for Figure 5.4
  
  test_arnoldi_matfun_ns_times.m for Figure 5.5
  
  test_arnoldi_matfun_graph.m for Figure 5.6 and Figure 5.7

[1]: Laura Grigori, Daniel Kressner, Nian Shao, Igor Simunec. ``Restoring similarity in randomized Krylov methods with applications to eigenvalue problems and matrix functions.'' 
