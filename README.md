**#CanteraTransportFoam**

Cantera software has been coupled with OpenFOAM, and a new mixture-averaged diffusion model with Soret effect has been developed, where both mass and thermal diffusion coefficients were calculated inside Cantera. Additionally more robust correction velocity approach has been deployed to balance the imbalance in mass fluxes. The included report has an overview and implementation of the existing thermophysical transport models in OpenFOAM, which at least to the author’s knowledge, are not listed openly anywhere.

Originally developed and published by Aleksi Rintanen as a project "Coupling OpenFOAM with Cantera for advanced diffusion modeling", during the Chalmers University course [CFD with OpenSource Software (2024)], held by Prof. Håkan Nilsson [https://www.tfd.chalmers.se/~hani/kurser/OS_CFD/].

The present version is an **updated version**. Some bugs from the earlier version are fixed and a zero mass flux boundary condition [zeroMassFlux] at walls for species is implemented. This prevents species from exiting the computational domain through the walls due to the Soret effect.

**Publications using the CanteraTransportFoam:**

S. Karimkashi, P. Tamadonfar, A. Rintanen, T. Zirwes, A. Stagni, V. Vuorinen, O. Kaario, Humidification effects on near-wall NO and N2O formation in NH3/H2 flames, Proceedings of the Combustion Institute, 2026, https://doi.org/10.1016/j.proci.2026.106707

V.-P. Salomaa, P. Tamadonfar, M. Gadalla, V. Vuorinen, O. Kaario, A numerical study of ammonia-hydrogen flames: Propagation and emission formation in crevices, Proceedings of the Combustion Institute, 2026, https://doi.org/10.1016/j.proci.2026.106112

G. Mehdi, E. Almusa, V.-P. Salomaa, M.G. De Giorgi, O. Kaario, Multiscale numerical study of Plasma-Assisted combustion in ammonia and Ammonia–Hydrogen blends, Fuel, 2027, https://doi.org/10.1016/j.fuel.2026.141454
