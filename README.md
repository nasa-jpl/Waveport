# Waveport
Waveport Scattering Library - Version 1

## Contact
Mark Haynes

Jet Propulsion Laboratory, California Institute of Technology

mark.s.haynes@jpl.nasa.gov

## Intro
The Waveport Scattering Library is a collection of equations, derivations, and Matlab codes on selected topics in electromagnetic scattering. This work started with a personal collection of notes and codes with the idea of creating a consistent scattering code library with which analysis and new routines could be developed and checked more quickly. The documentation was designed to serve as part quick reference, part code explanation, and part teaching guide.

## Topics
The topics of the library follow closely the textbooks of Chew, Tsang, Kong, and Ulaby: 
- Green’s functions
- Spherical wave functions
- Spherical wave function rotation and translation
- S-matrix and T-matrix scattering
- Facet scattering under the Kirchhoff approximation
- Fast Multipole Method operators
- Generation of rough surfaces and random objects
- Geometric (ray-tracing) scattering solutions
- Coordinate transforms and special functions

## Document
The main document is Waveport.pdf. It is written in the style of Numerical Recipes. Equations, derivations, explanations and code are combined inline in the text. All the Latex source files are in Tex/.

## Code 
Matlab routines are in directory Code/. Most routines are building blocks for analysis rather than complete heavy-hitting numerical codes. Example scripts for running the routines are included in the repository under each topic directory.

## License
The code is released under the Apache 2.0 Open Source license. It is free to run and modify for non-commercial use.

## Citation
This work can be cited generally as:

Haynes, M. (2021). *Waveport Scattering Library*. Jet Propulsion Laboratory - California Institute of Technology. https://doi.org/10.48588/JPL.HE9D-BA55

Specific sections and subsections can be cited as, for example:

Haynes, M. (2021). "Spherical Harmonics." *Waveport Scattering Library*. Jet Propulsion Laboratory - California Institute of Technology. https://doi.org/10.48588/JPL.HE9D-BA55


## Copyright
Copyright (c) 2020-21, by the California Institute of Technology. ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology. JPL document clearance CL#21-3866.

## Export
This software may be subject to U.S. export control laws. By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. User has the responsibility to obtain export licenses, or other export authority as may be required before exporting such information to foreign countries or providing access to foreign persons.

