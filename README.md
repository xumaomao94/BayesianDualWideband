# Bayesian modeling for dual-wideband channel estimation

This is the source code for our TSP paper "Overcoming Beam Squint in mmWave MIMO Channel Estimation: A Bayesian Multi-Band Sparsity Approach" [1]. You can start by running the demo code. Apart from the proposed Variational EM algorithm, my own implementations of OMP and TensorEsprit [2] are also included.

Remember to cite our paper [1] if you use the code.

## Functions

- Channel_build.m

  Channel setup

- f_em_doa_mimo/

  Codes for the proposed method

  - em_offgrid_dualwideband.m
 
  Use this method to test the proposed algorithm

- f_somp/

  OMP for solving dual-wideband channel estimation

- f_tensor_esprit/

  Tensor Esprit for solving the CE problem

## Reference

[1] Xu, Le, et al. "Overcoming Beam Squint in mmWave MIMO Channel Estimation: A Bayesian Multi-Band Sparsity Approach." IEEE Transactions on Signal Processing 72 (2024): 1219-1234.

[2] Y. Lin, S. Jin, M. Matthaiou, and X. You, “Tensor-based channel estimation for millimeter wave MIMO-OFDM with dual-wideband effects,” IEEE Transactions on Communications, vol. 68, no. 7, pp. 4218–4232, 2020.
