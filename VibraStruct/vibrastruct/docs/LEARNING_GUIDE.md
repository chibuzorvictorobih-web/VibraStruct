# Understanding and presenting VibraStruct

## What each term means
Mass stores kinetic energy; stiffness stores elastic energy; viscous damping removes mechanical energy. A natural frequency is a free-vibration frequency. A mode shape describes the relative displacements in that vibration pattern. Mode-shape sign and magnitude are arbitrary until a normalization is specified; this project uses mass normalization for calculations and peak normalization for plots.

Resonance occurs when forcing couples strongly to a mode near its natural frequency. The displacement peak of a damped SDOF system is slightly below the undamped natural frequency: omega_peak = omega_n sqrt(1 - 2 zeta²), valid for zeta < 1/sqrt(2). Damping reduces peaks and dissipates energy. It does not change the undamped eigenproblem used here.

## Read the figures
1. Resonance: compare peak heights at each damping ratio.
2. Modes: identify masses moving together or in opposition. Plot lines connect mass stations; they are not a continuous beam deformation.
3. Transient: follow the initial 2 mm disturbance at mass 3 as it travels through the coupled system.
4. Energy: undamped energy should remain constant within integration error; viscous damping reduces it.
5. Coupled response: understand why force location and response location affect peak amplitudes.
6. Sensitivity: all three curves follow square-root stiffness scaling because all springs are scaled together.

## Suggested hands-on work before using this in an application
- Reproduce the baseline and compare summary.json with the supplied version.
- Double all masses and predict the frequency change before running it.
- Change only spring 2; explain why the simple global scaling rule no longer applies.
- Change the force to mass 1 and compare modal participation.
- Reduce integration step size and check energy drift.
- Add an investigation you can explain and retain your own results and notes.

## Questions to be ready to answer
- How is each spring assembled into the stiffness matrix?
- Why solve a generalized eigenproblem instead of taking eigenvalues of K alone?
- Why is modal mass normalization useful?
- Why do the middle-mode damping ratio and the endpoint ratios differ under Rayleigh damping?
- What distinguishes numerical verification from experimental validation?
- What measurements would be needed for a real support structure?

## Equations to derive independently
For SDOF: omega_n = sqrt(k/m), c = 2 zeta m omega_n and X/F = 1/(k - m omega² + i c omega). For two equal masses and equal springs in this fixed-free topology, omega² = (k/m)(3 ± sqrt(5))/2. These provide independent reference checks.

## Reading
Standard background: Daniel J. Inman, Engineering Vibration; Anil K. Chopra, Dynamics of Structures. Consult the official NumPy, SciPy (scipy.linalg.eigh and scipy.integrate.solve_ivp) and Matplotlib documentation for API details. These are further reading, not evidence of experimental validation.
