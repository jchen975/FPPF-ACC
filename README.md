# FPPF-ACC
MATLAB code for the ACC2023 paper on fixed-point power flow: https://ieeexplore.ieee.org/document/10156226

If you find the results interesting/helpful, please cite the paper as follows: 

*L. Chen and J. W. Simpson-Porco, "A Fixed-Point Algorithm for the AC Power Flow Problem," 2023 American Control Conference (ACC), San Diego, CA, USA, 2023, pp. 4449-4456, doi: 10.23919/ACC55779.2023.10156226.*

Requirements: 
1. MATPOWER library and test cases: https://matpower.org/
2. LUQ decomposition: https://www.mathworks.com/matlabcentral/fileexchange/11120-null-space-of-a-sparse-matrix (in `sparse_null` folder)

To run: 
1. `full_fppf` runs FPPF with flat start conditions
2. `full_fppf_IC` runs the robustness test for FPPF, NR and FDLF
3. `full_fppf_RX` runs the convergence test for the network R/X ratio upper limit
4. `full_fppf_update_order` investigates the effect of update order on the convergence of FPPF

Please contact me at liangjie.chen@mail.utoronto.ca if you have any questions or bug report -- this repository is not being actively maintained. 
