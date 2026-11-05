# Review notes — version 1.1

The baseline equations and original reported natural frequencies were confirmed. This revision improves robustness and verification; it does not change the physical model.

## Changes
- Reject exactly singular undamped SDOF resonance with an explanatory error, instead of producing a non-finite complex result.
- Use scale-aware damping symmetry and nonnegative-dissipation checks, including very small negative damping matrices.
- Validate mechanical-energy history shapes and harmonic inputs.
- Expose integration tolerance and maximum-step controls; factor the constant mass matrix once per transient solve.
- Provide reusable coupled harmonic response, preserving complex force and response amplitudes.
- Make mode signs deterministic for reproducible plots.
- Verify coupled transients against the matrix exponential and coupled harmonic response against modal superposition.
- Add refined-integration results, actual modal damping ratios, harmonic phase data and machine-readable pass/fail checks.
- Generate baseline model values in the report from the variables used for calculation.
- Record exact scientific package versions used for the supplied results.

## Evidence and limits
Twenty unit tests pass. The complete analysis additionally checks modal residuals, orthogonality, undamped energy drift, damped energy decay, independent transient accuracy and agreement between two harmonic solution methods. Read results/summary.json and results/verification.json for measured errors and outcomes. A refinement study is numerical verification, not experimental validation. Neither this review nor passing tests establishes correctness for every possible input or makes the model suitable for industrial design decisions.
