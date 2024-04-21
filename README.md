# TPMconvex
Thermophysical Model (TPM) for convex asteroid shapes

Code to calculate surface temperatures on asteroids and thermally emitted radiation.
For a detailed description, see the PhD thesis of the developer, https://arxiv.org/abs/1208.3993 (defended 2007).

The code allows for (nearly) arbitrary asteroid shapes and spin axes. Thermal conduction into the subsoil is modeled explicitly, as is multiple scattering of visible and IR photons within small concavities ("surface roughness").  This code was used in a variety of scientific studies. It has been validated extensively against observational data and against similar codes by other authors.  

The code is organized in several subfolders of folder 'source'. The most immediately useful project to start exploring the possibilities of the TPM is probably source/ThermalLC that calculates and outputs thermal lightcurves for a given shape, spin state, observational circumstances, and assumed physical properties.
