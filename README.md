# Nonlinear Spring and Axial-Angular Spring Restraints for OpenFOAM

## Related article

These directories and files are related to the following article:

Poorya Poozesh, Antonio J. Álvarez, Arturo N. Fontán, and Félix Nieto, *Implementation and verification of nonlinear spring-damper restraints in OpenFOAM's six-DOF rigid-body motion framework*, OpenFOAM Journal, 2026.

## Copyright and licensing

Copyright © 2026 Poorya Poozesh, Antonio J. Álvarez, Arturo N. Fontán, and Félix Nieto.

The article and its accompanying documentation are made available under the Creative Commons Attribution-ShareAlike 4.0 International license (CC BY-SA 4.0). The source code is distributed under the GNU General Public License, version 3 or later, as stated in the individual source files. Existing third-party copyright and license notices remain in effect and take precedence for the files to which they apply.

## Overview

This repository contains two custom six-degree-of-freedom (six-DOF) restraint libraries for OpenFOAM:

- `nonLinearSpring`: nonlinear translational spring-damper restraint.
- `nonLinearAxialAngularSpring`: nonlinear axial/angular spring-damper restraint.

Both classes derive from `sixDoFRigidBodyMotionRestraint` and are available through the OpenFOAM runtime-selection table.

Supported spring models are:

- `linear`: `F = A*x`, where `x` is displacement.
- `Duffing`: `F = A*x + B*x^3`.
- `stiffnessList`: a tabulated displacement-force or rotation-moment relation.

Supported damper models are:

- `linear`: `F = c1*v`, where `v` is velocity.
- `Rayleigh`: `F = c1*v + c2*v^3`.
- `Mann`: `F = c1*v + c2*x^2*v`.
- `dampingList`: a tabulated displacement-force or velocity-force relation and its rotational equivalent.

Tabulated stiffness and damping data can use linear or cubic-spline interpolation. The repository also includes two libraries that provide linear extrapolation outside the tabulated range:

- `libinterExtrapolationTable`: linear interpolation with linear extrapolation.
- `libinterExtrapolateSplineXY`: cubic-spline interpolation with linear extrapolation.

## Compatibility and requirements

- OpenFOAM ESI must be loaded before compiling or running a case. The supplied code has been tested with OpenFOAM v2212 and v2406.
- The tutorial cases are configured for eight subdomains in `system/decomposeParDict`.
- The plotting scripts require Python 3, NumPy, and Matplotlib. The `os`, `sys`, and `copy` modules used by the scripts are part of the Python standard library.
- Run the plotting scripts in a normal Python environment after the OpenFOAM run. A Python installation loaded inside the OpenFOAM environment can conflict with OpenFOAM's MPI libraries.

## Repository layout

- `src/`: restraint and interpolation/extrapolation source code, plus `Allwmake` and `Allwclean`.
- `tutorials/Translational/freeDecay/`: cases 01-04.
- `tutorials/Translational/wind/`: cases 05-07.
- `tutorials/Rotational/freeDecay/`: cases 08-10.
- `tutorials/Rotational/wind/`: case 11.
- `tutorials/FieldsForMapping/`: steady fields mapped by the wind-excited wing cases to reduce the initial flow transient.
- `tutorials/FSI_applicationExample/`: case 12, the 3:2 rectangular-prism application.

## Compile the libraries

From a shell in which OpenFOAM ESI is loaded:

```sh
cd src
./Allwmake
```

The compiled libraries are written to `$FOAM_USER_LIBBIN`. To remove the installed libraries and generated `lnInclude` directories, keep OpenFOAM loaded and run:

```sh
./Allwclean
```

## Run and post-process a case

Each case contains an `Allrun` script that generates the mesh, prepares the initial fields, decomposes the domain, and runs the configured solver in parallel. Wind-excited wing cases also map the initial solution from `tutorials/FieldsForMapping`.

From the selected case directory:

```sh
./Allrun
```

After the simulation finishes, run the case-specific plotting script from the same directory, preferably in a separate shell with a standard Python environment:

```sh
python3 plotGeneration.py
```

The script reads the rigid-body state from `postProcessing/sixDoFRigidBodyState/0/sixDoFRigidBodyState.dat`. Wind-excited cases also read `postProcessing/forceCoeffsWing/0/coefficient.dat`. If OpenFOAM writes to a start-time directory other than `0`, update those paths in `plotGeneration.py`.

Every supplied plotting script writes:

- `image.pdf`
- `image.png`

For free-decay cases, these files contain the displacement or rotation time history. For wind-excited cases, they contain the motion and the relevant force coefficient (`Cl` for heave and `Cm` for pitch). The case-12 plot contains `y/H` and `Cl` histories.

## Case-to-manuscript mapping

The mapping below refers to the manuscript titled *Implementation and verification of nonlinear spring-damper restraints in OpenFOAM's six-DOF rigid-body motion framework*. A case reproduces the OpenFOAM-library contribution to the cited result. Comparative preCICE, standard-`linearSpring`, analytical, and central-finite-difference datasets are not included in this archive unless stated otherwise.

### Free-decay and verification cases

| Case | Manuscript result | Active DOF | Spring and damper setup | Interpolation | Default `deltaT` | Main supplied output |
| --- | --- | --- | --- | --- | --- | --- |
| `01_wing_linearSpring_linearDamping` | Section 4.1; Figure 4; Tables 2, 3, and 6 | Heave, global `y` translation (reported as `x` in the manuscript) | Linear stiffness represented by `stiffnessList`, `k = 4000 N/m`; linear damping, `c1 = 6 N s/m` | Linear; `EXTRAPOLATE` | `1.0e-5 s` | Translational free-decay history in `image.pdf` and `image.png` |
| `02_wing_listSpring_linearDamping` | Companion free-decay case for the tabulated nonlinear law used in Section 5.2 and Figure 13; no dedicated response figure or table | Heave, global `y` translation | Nonlinear `stiffnessList`; linear damping, `c1 = 6 N s/m` | Linear; `EXTRAPOLATE` | `1.0e-6 s` | Translational free-decay history; the stiffness list itself defines the Figure 13 law |
| `03_wing_DuffingSpring_noDamping` | Section 4.2; Figure 6; Tables 7, 8, and 11 | Heave, global `y` translation | Duffing stiffness, `A = 4000 N/m`, `B = 200000 N/m^3`; zero damping (`c1 = c2 = 0`) | Not applicable | `1.0e-6 s` | Translational undamped history in `image.pdf` and `image.png` |
| `04_wing_DuffingSpring_MannDamping` | Section 4.3; Figure 8; Tables 12 and 14 | Heave, global `y` translation | Duffing stiffness, `A = 4000 N/m`, `B = 200000 N/m^3`; Mann damping, `c1 = 6`, `c2 = 5000` | Not applicable | `1.0e-6 s` | Translational nonlinear damped history; the central-difference reference is not generated by the supplied script |
| `08_wing_linearSpring_linearDamping` | Section 4.1; Figure 5; Tables 4, 5, and 6 | Pitch, rotation about global `z` | Linear stiffness and damping represented by lists, `k = 8000 N m/rad`, `c = 20 N m s/rad` | Linear; `EXTRAPOLATE` | `1.0e-5 s` | Rotational free-decay history in `image.pdf` and `image.png` |
| `09_wing_Duffing_noDamping` | Section 4.2; Figure 7; Tables 9, 10, and 11 | Pitch, rotation about global `z` | Duffing stiffness, `A = 8000 N m/rad`, `B = 200000 N m/rad^3`; zero damping (`c1 = c2 = 0`) | Not applicable | `1.0e-6 s` | Rotational undamped history in `image.pdf` and `image.png` |
| `10_wing_listSpring_RayleighDamping` | Section 4.3; Figures 9 and 10; Tables 13 and 14 | Pitch, rotation about global `z` | Nonlinear `stiffnessList`; Rayleigh damping, `c1 = 20`, `c2 = 50` | Linear; `EXTRAPOLATE` | `1.0e-6 s` | Rotational nonlinear damped history; the central-difference reference is not generated by the supplied script |

### Wind-excited and application cases

| Case | Manuscript result | Active DOF | Spring and damper setup | Interpolation | Default `deltaT` | Main supplied output |
| --- | --- | --- | --- | --- | --- | --- |
| `05_wing_linearSpring_linearDamping` | Section 5.1; Figures 11 and 12 | Heave, global `y` translation (reported as `x` in the manuscript) | Linear stiffness and velocity damping represented by lists, `k = 4000 N/m`, `c = 2 N s/m` | Linear; `EXTRAPOLATE` | `1.0e-5 s` | Non-dimensional heave and lift-coefficient histories; preCICE and standard-`linearSpring` series are not included |
| `06_wing_listSpring_linearDamping` | Section 5.2; Figures 13 and 14 | Heave, global `y` translation | Nonlinear `stiffnessList`; linear damping, `c1 = 2 N s/m` | Linear; `EXTRAPOLATE` | `1.0e-6 s` | Non-dimensional heave and lift-coefficient histories; the preCICE series is not included |
| `07_wing_DuffingSpring_MannDamping` | Section 5.3; Figures 15 and 16; Table 15 | Heave, global `y` translation | Duffing stiffness, `A = 4000 N/m`, `B = 200000 N/m^3`; Mann damping, `c1 = 6`, `c2 = 5000` | Not applicable | `1.0e-5 s` | Non-dimensional heave and lift-coefficient histories |
| `11_wing_listSpring_RayleighDamping` | Section 5.4; Figure 17 | Pitch, rotation about global `z` | Nonlinear `stiffnessList`; Rayleigh damping, `c1 = 20`, `c2 = 50` | Linear; `EXTRAPOLATE` | `1.0e-6 s` | Pitch and moment-coefficient histories |
| `12_prism3_2_listSpring_linearDamping` | Section 6; setup in Figure 2 and Tables 16-19; results in Figures 18-23 and Tables 20-21 | Heave, global `y` translation | Nonlinear `stiffnessList`; linear damping, `c1 = 4.594375 N s/m` | Linear; `EXTRAPOLATE` | `8.0e-4 s` | `y/H` and lift-coefficient histories; mesh, spectra, phase averages, statistics, Q criterion, and preCICE comparison require additional post-processing/data |

### Time-step sensitivity runs

The archive provides one default `deltaT` per case. The manuscript figures and tables listed below combine multiple simulations. To reproduce them, copy the relevant case to a separate working directory for each time step, change `deltaT` in `system/controlDict`, and run `./Allclean` followed by `./Allrun`. Preserve each run's `postProcessing` directory and log files before starting another run.

| Manuscript result | Base cases | Required time steps |
| --- | --- | --- |
| Figures 4 and 5; Tables 3, 5, and 6 | 01 and 08 | `7.2e-5`, `1.0e-5`, and `1.0e-6 s` |
| Figures 6 and 7; Tables 8, 10, and 11 | 03 and 09 | `7.2e-5`, `1.0e-5`, and `1.0e-6 s` |
| Figure 11 | 05 | `7.2e-5`, `1.0e-5`, and `1.0e-6 s` |
| Figure 14 | 06 | `1.0e-5` and `1.0e-6 s` for the OpenFOAM `nonLinearSpring` curves |
| Figure 15 | 07 | `7.2e-5`, `1.0e-5`, and `1.0e-6 s` |

Cases 04 and 10 use `1.0e-6 s` for the Section 4.3 numerical-verification comparison. The application-case values shown in the mapping table are the defaults supplied in the archive.

## FSI application resource note

`tutorials/FSI_applicationExample/12_prism3_2_listSpring_linearDamping` is substantially more expensive than the wing tutorials. Its `Allrun` script uses the number of subdomains set in `system/decomposeParDict` (eight in the supplied archive), but a publication-quality production run may require more cores and several days. Adjust the decomposition to the available hardware and record the OpenFOAM version, core count, time step, end time, and any dictionary changes with the reproduced results.

## Disclaimer

The nonlinear spring-damper libraries are provided as is, without warranty. The developers are not liable for errors, failures, or damages arising from their use. Users assume full responsibility for applying the code in their own work.
