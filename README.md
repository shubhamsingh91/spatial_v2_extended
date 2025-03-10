# spatial_v2_extended

This package directly builds upon Roy Featherstone's [spatial_v2 library](http://royfeatherstone.org/spatial/v2/) and his closely associated [book](https://link.springer.com/book/10.1007/978-1-4899-7560-7).

New algorithms include:
* Methods to compute the Coriolis matrix ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/blob/main/v3/dynamics/CoriolisMatrix.m)) [all systems] and Christoffel symbols ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/blob/main/v3/dynamics/Christoffel.m)) ([paper](http://dx.doi.org/10.1115/1.4051169))
* Methods for assessing identifiability ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/identifiability), [paper](https://arxiv.org/abs/1711.03896)) 
* Methods to calculate second-order partial derivatives of Inverse Dynamics for multi-DoF joints ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/blob/main/v3/derivatives/ID_SO_derivatives.m)) 
* Methods to calculate second-order partial derivatives of Forward Dynamics for multi-DoF joints ([link](https://github.com/shubhamsingh91/spatial_v2_extended/blob/main/v3/derivatives/FD_SO_derivatives.m)) ([paper](https://arxiv.org/abs/2302.06001))
* Methods to calculate First/Second-order partial derivatives of KKT Forward Dynamics and Impact Dynamics for multi-DoF joints ([link](https://github.com/shubhamsingh91/spatial_v2_extended/blob/main/v3/unit_tests/UnitTest_Derivatives_KKT.m)) 

New features include:
* Extensions of most algorithms (RNEA, ABA, CRBA, etc.) to address dynamic effects from motor rotors ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/dynamics))
* Updated most algorithms (RNEA, ABA, CRBA, etc.) to handle multi-DoF joints (e.g., spherical or floating base) without needing specialized versions of the algorithms ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/dynamics))
* Regressor calculation algorithms ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/regressor))
* Variety of tools for converting between different representations of orientation ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/orientation_tools))
* Methods for computing the partial derivatives of Inverse Dynamics ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/blob/main/v3/derivatives/ID_derivatives.m)) 
* Partial compatibility for complex-valued input arguments toward support of complex-step derivative comptuations in unit tests (including the complex step on matrix Lie groups [link](https://ieeexplore-ieee-org.proxy.library.nd.edu/abstract/document/8957301)).
