"""Linear structural dynamics; SI units throughout."""
import numpy as np
from scipy.linalg import eigh, cho_factor, cho_solve
from scipy.integrate import solve_ivp


def positive(value, name):
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f'{name} must be positive and finite')


def sdof(mass, stiffness, damping_ratio):
    """Return natural angular frequency and viscous damping coefficient."""
    positive(mass, 'mass'); positive(stiffness, 'stiffness')
    if not np.isfinite(damping_ratio) or damping_ratio < 0:
        raise ValueError('damping_ratio must be finite and nonnegative')
    wn = np.sqrt(stiffness / mass)
    return wn, 2 * damping_ratio * mass * wn


def receptance(mass, stiffness, damping_ratio, omega):
    """Complex displacement/force, using exp(+i omega t), in m/N."""
    _, c = sdof(mass, stiffness, damping_ratio)
    omega = np.asarray(omega, dtype=float)
    if np.any(~np.isfinite(omega)) or np.any(omega < 0):
        raise ValueError('frequencies must be finite and nonnegative')
    denominator = stiffness - mass * omega**2 + 1j * c * omega
    if np.any(denominator == 0):
        raise ValueError('undamped resonance has no finite steady-state response')
    return 1 / denominator


def chain(masses, springs):
    """Fixed-left/free-right chain: n masses and n springs, first grounded."""
    masses = np.asarray(masses, dtype=float); springs = np.asarray(springs, dtype=float)
    if masses.ndim != 1 or masses.size == 0 or springs.shape != masses.shape:
        raise ValueError('masses and springs must be equal nonempty 1D arrays')
    if np.any(~np.isfinite(masses)) or np.any(masses <= 0) or np.any(~np.isfinite(springs)) or np.any(springs <= 0):
        raise ValueError('masses and springs must be finite and positive')
    n = len(masses); K = np.zeros((n,n)); K[0,0] = springs[0]
    for j in range(1,n):
        K[j-1,j-1] += springs[j]; K[j,j] += springs[j]
        K[j-1,j] -= springs[j]; K[j,j-1] -= springs[j]
    return np.diag(masses), K


def matrices(M, K):
    M = np.asarray(M, dtype=float); K = np.asarray(K, dtype=float)
    if M.ndim != 2 or M.shape[0] == 0 or M.shape[0] != M.shape[1] or K.shape != M.shape:
        raise ValueError('M and K must be equally sized square matrices')
    for A in (M,K):
        if not np.all(np.isfinite(A)) or not np.allclose(A,A.T,rtol=1e-12,atol=1e-12):
            raise ValueError('M and K must be finite and symmetric')
        try: np.linalg.cholesky(A)
        except np.linalg.LinAlgError as e: raise ValueError('M and K must be positive definite') from e
    return M,K


def modes(M,K):
    """Ascending angular frequencies and mass-normalized mode columns."""
    M,K = matrices(M,K)
    values, phi = eigh(K,M)
    if np.any(values <= 0):
        raise ValueError('nonpositive eigenvalues; check conditioning and constraints')
    # Deterministic signs simplify reproducible plots; eigenvector sign is arbitrary.
    for j in range(phi.shape[1]):
        if phi[np.argmax(np.abs(phi[:, j])), j] < 0:
            phi[:, j] *= -1
    return np.sqrt(values), phi


def rayleigh(M,K,zeta=0.02):
    """Fit equal damping ratio at lowest and highest natural frequencies."""
    M,K = matrices(M,K)
    if not np.isfinite(zeta) or zeta < 0: raise ValueError('zeta must be nonnegative')
    w,_ = modes(M,K)
    beta = 2*zeta/(w[0]+w[-1]); alpha = beta*w[0]*w[-1]
    return alpha*M + beta*K


def damping_matrix(C, shape):
    """Validate symmetry and nonnegative dissipation at floating-point precision."""
    C = np.asarray(C, dtype=float)
    if C.shape != shape or not np.all(np.isfinite(C)):
        raise ValueError('C must be finite and match M')
    scale = max(float(np.linalg.norm(C, ord=2)), np.finfo(float).tiny)
    if np.max(np.abs(C-C.T)) > 1e-12*scale:
        raise ValueError('C must be symmetric')
    C = (C+C.T)/2
    if np.min(np.linalg.eigvalsh(C)) < -100*np.finfo(float).eps*scale:
        raise ValueError('C must be positive semidefinite')
    return C


def harmonic_response(M, K, C, omega, force):
    """Complex displacement vector for each angular frequency (rad/s).

    Uses exp(+i omega t); force is a constant complex amplitude vector in N.
    Returns an array of shape (number of frequencies, number of masses).
    A singular undamped resonance raises ValueError.
    """
    M, K = matrices(M, K)
    C = damping_matrix(C, M.shape)
    omega = np.atleast_1d(np.asarray(omega, dtype=float))
    force = np.asarray(force, dtype=complex)
    if omega.ndim != 1 or not len(omega) or not np.all(np.isfinite(omega)) or np.any(omega < 0):
        raise ValueError('omega must be a nonempty finite nonnegative vector')
    if force.shape != (len(M),) or not np.all(np.isfinite(force)):
        raise ValueError('force must be a finite amplitude vector matching M')
    try:
        return np.array([np.linalg.solve(K-w*w*M+1j*w*C, force) for w in omega])
    except np.linalg.LinAlgError as exc:
        raise ValueError('singular dynamic stiffness: no unique finite harmonic response') from exc


def response(M,K,C,t,x0,v0,force=None,*,rtol=1e-9,atol=1e-11,max_step=None):
    """Integrate M xdd + C xd + K x = f(t); output x and v at t.

    Initial conditions apply at t[0]. Force is a callable returning an n-vector.
    The default maximum internal step is the smallest output interval.
    Override max_step (seconds), rtol and atol for convergence studies.
    """
    M,K = matrices(M,K); n = len(M); C = damping_matrix(C, M.shape)
    positive(rtol, 'rtol'); positive(atol, 'atol')
    t = np.asarray(t,dtype=float); x0=np.asarray(x0,dtype=float); v0=np.asarray(v0,dtype=float)
    if t.ndim != 1 or len(t)<2 or not np.all(np.isfinite(t)) or np.any(np.diff(t)<=0):
        raise ValueError('time must be finite and strictly increasing')
    if x0.shape != (n,) or v0.shape != (n,) or not np.all(np.isfinite(np.r_[x0,v0])):
        raise ValueError('initial conditions must be finite vectors matching M')
    if force is not None and not callable(force):
        raise ValueError('force must be callable or None')
    if max_step is None:
        max_step = float(np.min(np.diff(t)))
    positive(max_step, 'max_step')
    factor = cho_factor(M)
    def rhs(time,y):
        f = np.zeros(n) if force is None else np.asarray(force(time),dtype=float)
        if f.shape != (n,) or not np.all(np.isfinite(f)): raise ValueError('force must return a finite n-vector')
        return np.r_[y[n:],cho_solve(factor,f-C@y[n:]-K@y[:n])]
    sol=solve_ivp(rhs,(t[0],t[-1]),np.r_[x0,v0],t_eval=t,rtol=rtol,atol=atol,max_step=max_step)
    if not sol.success: raise RuntimeError(sol.message)
    return sol.y[:n].T,sol.y[n:].T


def energy(M,K,x,v):
    """Mechanical energy [J] for histories of shape (samples, masses)."""
    M,K = matrices(M,K)
    x=np.asarray(x,dtype=float); v=np.asarray(v,dtype=float)
    if x.ndim != 2 or x.shape != v.shape or x.shape[1] != len(M) or not np.all(np.isfinite(x)) or not np.all(np.isfinite(v)):
        raise ValueError('x and v must be finite matching histories with n columns')
    return 0.5*np.einsum('ti,ij,tj->t',v,M,v)+0.5*np.einsum('ti,ij,tj->t',x,K,x)
