"""Reproduce all figures, tables and the technical report."""
from pathlib import Path
import json
import platform
import numpy as np
import scipy
from scipy.linalg import expm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from vibrastruct import sdof,receptance,chain,modes,rayleigh,response,energy,harmonic_response


def main():
    out=Path(__file__).resolve().parent/'results'; out.mkdir(exist_ok=True)
    plt.rcParams.update({'figure.dpi':140,'axes.grid':True,'grid.alpha':0.25})
    def save(name):
        plt.tight_layout(); plt.savefig(out/name); plt.close()
    m,k,z=250.,1e6,0.03
    wn,c=sdof(m,k,z); f=np.linspace(0,30,2401); omega=2*np.pi*f
    plt.figure(figsize=(8,4.5))
    for damping in [0.01,0.03,0.10]:
        plt.plot(f,1000*100*np.abs(receptance(m,k,damping,omega)),label=f'Damping ratio {damping:.0%}')
    plt.axvline(wn/(2*np.pi),color='k',ls='--',label='Undamped natural frequency')
    plt.xlabel('Excitation frequency [Hz]'); plt.ylabel('Displacement amplitude [mm], 100 N force')
    plt.title('SDOF resonance and damping'); plt.legend(); save('01_resonance.png')
    np.savetxt(out/'sdof_frequency_response.csv',np.c_[f,np.abs(receptance(m,k,z,omega)),np.angle(receptance(m,k,z,omega))],delimiter=',',header='frequency_Hz,receptance_m_per_N,phase_rad',comments='')
    masses=[250,180,120]; springs=[1e6,0.8e6,0.6e6]
    M,K=chain(masses,springs); w,P=modes(M,K); C=rayleigh(M,K,0.03)
    t=np.linspace(0,3,3001); x0=np.array([0,0,0.002]); v0=np.zeros(3)
    x,v=response(M,K,C,t,x0,v0); xu,vu=response(M,K,np.zeros_like(M),t,x0,v0)
    E=energy(M,K,x,v); Eu=energy(M,K,xu,vu)
    fig,axs=plt.subplots(1,3,figsize=(10,4),sharey=True)
    for j,ax in enumerate(axs):
        shape=P[:,j]/np.max(np.abs(P[:,j])); ax.plot(np.r_[0,shape],range(4),'o-')
        ax.set_title(f'Mode {j+1}: {w[j]/(2*np.pi):.2f} Hz'); ax.set_xlabel('Relative displacement'); ax.set_yticks(range(4))
    axs[0].set_ylabel('Mass station (0 = ground)'); save('02_modes.png')
    plt.figure(figsize=(8,4.5))
    for j in range(3): plt.plot(t,x[:,j]*1000,label=f'Mass {j+1}')
    plt.xlabel('Time [s]'); plt.ylabel('Displacement [mm]'); plt.title('Damped free vibration'); plt.legend(); save('03_transient.png')
    plt.figure(figsize=(8,4.5)); plt.plot(t,E,label='Damped'); plt.plot(t,Eu,'--',label='Undamped')
    plt.xlabel('Time [s]'); plt.ylabel('Mechanical energy [J]'); plt.title('Energy balance'); plt.legend(); save('04_energy.png')
    freq=np.linspace(0,30,3001); F=np.array([0,0,100.]); H=harmonic_response(M,K,C,2*np.pi*freq,F)
    plt.figure(figsize=(8,4.5))
    for j in range(3): plt.plot(freq,np.abs(H[:,j])*1000,label=f'Mass {j+1}')
    for a in w/(2*np.pi): plt.axvline(a,color='grey',ls=':',alpha=0.5)
    plt.xlabel('Excitation frequency [Hz]'); plt.ylabel('Amplitude [mm]'); plt.title('Coupled system: 100 N harmonic force at mass 3'); plt.legend(); save('05_coupled_response.png')
    scales=np.linspace(0.5,1.5,41); sweep=np.array([modes(M,s*K)[0]/(2*np.pi) for s in scales])
    plt.figure(figsize=(8,4.5))
    for j in range(3): plt.plot(scales,sweep[:,j],label=f'Mode {j+1}')
    plt.xlabel('Global stiffness multiplier'); plt.ylabel('Natural frequency [Hz]'); plt.title('Stiffness sensitivity'); plt.legend(); save('06_sensitivity.png')
    np.savetxt(out/'modes.csv',np.c_[np.arange(1,4),w/(2*np.pi),P.T],delimiter=',',header='mode,frequency_Hz,phi_mass1,phi_mass2,phi_mass3',comments='')
    np.savetxt(out/'transient.csv',np.c_[t,x,v,E],delimiter=',',header='time_s,x1_m,x2_m,x3_m,v1_m_s,v2_m_s,v3_m_s,energy_J',comments='')
    np.savetxt(out/'coupled_response.csv',np.c_[freq,np.abs(H)],delimiter=',',header='frequency_Hz,amplitude1_m,amplitude2_m,amplitude3_m',comments='')
    np.savetxt(out/'stiffness_sensitivity.csv',np.c_[scales,sweep],delimiter=',',header='stiffness_multiplier,f1_Hz,f2_Hz,f3_Hz',comments='')
    # Independent exact linear-system benchmark and a refined integration.
    A=np.block([[np.zeros_like(M),np.eye(len(M))],[-np.linalg.solve(M,K),-np.linalg.solve(M,C)]])
    sample_idx=np.arange(0,len(t),100)
    exact=np.array([expm(A*(t[j]-t[0]))@np.r_[x0,v0] for j in sample_idx])
    xr,vr=response(M,K,C,t,x0,v0,rtol=1e-11,atol=1e-13,max_step=0.0005)
    baseline_error=float(np.max(np.abs(x[sample_idx]-exact[:,:len(M)])))
    refined_error=float(np.max(np.abs(xr[sample_idx]-exact[:,:len(M)])))
    modal_zeta=np.diag(P.T@C@P)/(2*w)
    modal_H=np.array([P@((P.T@F)/(w*w-a*a+2j*modal_zeta*w*a)) for a in 2*np.pi*freq])
    modal_error=float(np.max(np.abs(H-modal_H))/np.max(np.abs(H)))
    np.savetxt(out/'modal_damping.csv',np.c_[np.arange(1,4),modal_zeta],delimiter=',',header='mode,damping_ratio',comments='')
    np.savetxt(out/'coupled_phase.csv',np.c_[freq,np.angle(H)],delimiter=',',header='frequency_Hz,phase1_rad,phase2_rad,phase3_rad',comments='')
    residual=np.linalg.norm(K@P-M@P@np.diag(w*w))/np.linalg.norm(K@P)
    metrics={'sdof_frequency_Hz':wn/(2*np.pi),'sdof_damping_Ns_m':c,'chain_frequencies_Hz':(w/(2*np.pi)).tolist(),'eigen_residual_relative':residual,'mass_orthogonality_error':float(np.linalg.norm(P.T@M@P-np.eye(3))),'undamped_energy_relative_drift':float(np.max(np.abs(Eu-Eu[0]))/Eu[0]),'damped_final_energy_fraction':float(E[-1]/E[0]),'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__}
    metrics.update({'modal_damping_ratios':modal_zeta.tolist(),
                    'baseline_displacement_error_vs_expm_m':baseline_error,
                    'refined_displacement_error_vs_expm_m':refined_error,
                    'max_displacement_change_on_refinement_m':float(np.max(np.abs(xr-x))),
                    'harmonic_modal_vs_direct_relative_error':modal_error})
    checks={'eigen_residual':residual<1e-10,
            'orthogonality':metrics['mass_orthogonality_error']<1e-10,
            'undamped_energy_drift':metrics['undamped_energy_relative_drift']<1e-6,
            'transient_exact_error':baseline_error<1e-8,
            'refined_exact_error':refined_error<1e-9,
            'harmonic_modal_agreement':modal_error<1e-10,
            'damped_energy_nonincrease':bool(np.max(np.diff(E))<1e-10)}
    checks={name:bool(passed) for name,passed in checks.items()}
    (out/'verification.json').write_text(json.dumps(checks,indent=2)+'\n', encoding='utf-8')
    if not all(checks.values()): raise RuntimeError(f'Verification failed: {checks}')
    (out/'summary.json').write_text(json.dumps(metrics,indent=2)+'\n', encoding='utf-8')
    report=f'''# VibraStruct — Technical Report

## Purpose and scope
An educational numerical study of resonance, modal behaviour and damping in an idealised mechanical support system. Parameters are synthetic: no refinery measurements or proprietary data are used. This is a computational demonstration, not a calibrated model of industrial equipment.

## Models and assumptions
The SDOF baseline has mass {m:g} kg, stiffness {k:g} N/m and viscous damping ratio {z:g}. The three-degree-of-freedom chain has masses {masses} kg and springs {springs} N/m. Its topology is ground–spring 1–mass 1–spring 2–mass 2–spring 3–mass 3; the right end is free. Every mass translates along one axis. Small displacements, linear springs, constant lumped masses and viscous damping are assumed. Rotations, nonlinear contacts, foundation flexibility and fluid forces are excluded.

## Governing equations
M x'' + C x' + K x = f(t). Natural modes satisfy K phi = omega² M phi and are normalized so Phiᵀ M Phi = I. Frequency is omega/(2 pi). Rayleigh damping uses C = alpha M + beta K, with beta = 2 zeta/(omega_1+omega_n) and alpha = beta omega_1 omega_n. The computed modal damping ratios are {", ".join(f"{a:.4%}" for a in modal_zeta)}. The endpoint target is 3%; the middle mode differs. Harmonic response is X = [K - omega² M + i omega C]⁻¹ F. Mechanical energy is (vᵀ M v + xᵀ K x)/2.

## Numerical methods
SciPy's symmetric generalized eigensolver computes the modes. Adaptive RK45 directly integrates free response, with relative tolerance 1e-9, absolute tolerance 1e-11 and maximum step {np.min(np.diff(t)):.6f} s. A refined run uses relative tolerance 1e-11, absolute tolerance 1e-13 and maximum step 0.0005 s. The free response begins with displacements {x0.tolist()} m and velocities {v0.tolist()} m/s. A separate harmonic study applies the force amplitude vector {F.tolist()} N over {freq[0]:g}–{freq[-1]:g} Hz, with a cosine time convention. The SDOF sweep uses 0–30 Hz and three damping ratios. These are independent load cases.

## Results
The SDOF undamped natural frequency is **{wn/(2*np.pi):.4f} Hz**. The chain frequencies are **{', '.join(f'{a:.4f}' for a in w/(2*np.pi))} Hz**. Increasing damping suppresses resonance peaks. Coupling produces three distinct modes, with participation depending on where force and displacement are evaluated. Uniform stiffness scaling changes each frequency by the square root of the multiplier; doubling stiffness therefore increases frequency by sqrt(2). It does not double frequency.

![Resonance](../results/01_resonance.png)
![Modes](../results/02_modes.png)
![Transient](../results/03_transient.png)
![Energy](../results/04_energy.png)
![Coupled harmonic response](../results/05_coupled_response.png)
![Sensitivity](../results/06_sensitivity.png)

## Verification
Relative eigen-equation residual: {residual:.3e}. Mass orthogonality error: {metrics['mass_orthogonality_error']:.3e}. Maximum relative energy drift for the undamped run: {metrics['undamped_energy_relative_drift']:.3e}. The damped energy fraction after {t[-1]-t[0]:g} s is {E[-1]/E[0]:.6f}. The test suite additionally compares against the exact SDOF free-decay solution, a closed-form two-mass eigenvalue benchmark, energy conservation and the static frequency-response limit. These verify the implementation; they do not establish real-world model validity.

### Independent coupled-system checks
A matrix exponential of the constant first-order state matrix supplies a separate exact linear-system solution at 31 sampled times. Baseline maximum displacement error is {baseline_error:.3e} m; refined maximum error is {refined_error:.3e} m. The largest baseline/refined displacement difference over all output times is {metrics['max_displacement_change_on_refinement_m']:.3e} m. Modal superposition and direct dynamic-stiffness solutions agree to a maximum relative difference of {modal_error:.3e} over the harmonic sweep. Relative error here divides maximum absolute disagreement by maximum absolute response over all frequencies and masses. Acceptance thresholds and pass/fail outcomes are recorded in verification.json. These checks are useful because a self-consistent implementation can still contain modeling errors.

## Interpretation and limitations
Resonance depends on the relationship between excitation and system frequencies. Reducing vibration can involve changing stiffness, mass or damping, but this study does not select a safe operating speed for any actual machine. Harmonic peaks from a finite sweep are sampled estimates. Rayleigh damping is a mathematical approximation, not measured material damping. The lumped chain demonstrates structural dynamics and does not provide continuum stress or finite-element predictions. Real application would require measured geometry, supports, forcing and damping, followed by experimental validation.

## Reproducibility
Run `python run_analysis.py` from the project folder. CSV tables use SI units and summary.json records library versions. All results are generated from the supplied code using synthetic inputs. Numerical verification is documented above; experimental validation is outside the scope of this study.
'''
    (out.parent/'docs'/'TECHNICAL_REPORT.md').write_text(report, encoding='utf-8')
    print(json.dumps(metrics,indent=2))

if __name__=='__main__': main()
