"""Additional independent checks and edge cases from the project review."""
import unittest
import numpy as np
from scipy.linalg import expm
from vibrastruct import *

class ReviewValidation(unittest.TestCase):
    def test_undamped_resonance_is_explicit(self):
        with self.assertRaisesRegex(ValueError,'no finite'):
            receptance(1,4,0,2)
    def test_complex_harmonic_matches_sdof(self):
        w=np.array([0,4,20,40])
        result=harmonic_response([[2]],[[800]],[[8]],w,[3j])
        np.testing.assert_allclose(result[:,0],3j*receptance(2,800,0.1,w))
    def test_coupled_static_equilibrium(self):
        M,K=chain([2,3,4],[800,500,600]); F=np.array([0,0,5])
        H=harmonic_response(M,K,rayleigh(M,K),[0],F)
        np.testing.assert_allclose(K@H[0],F,atol=1e-13)
    def test_coupled_modal_harmonic(self):
        M,K=chain([2,3,4],[800,500,600]); w,P=modes(M,K); C=rayleigh(M,K)
        F=np.array([0,0,5]); om=np.linspace(0,35,101); z=np.diag(P.T@C@P)/(2*w)
        modal=np.array([P@((P.T@F)/(w*w-a*a+2j*z*w*a)) for a in om])
        np.testing.assert_allclose(harmonic_response(M,K,C,om,F),modal,rtol=1e-11,atol=1e-12)
    def test_coupled_transient_matrix_exponential(self):
        M,K=chain([2,3],[800,500]); C=rayleigh(M,K)
        A=np.block([[np.zeros((2,2)),np.eye(2)],[-np.linalg.solve(M,K),-np.linalg.solve(M,C)]])
        t=np.linspace(2,3,101); initial=np.array([0,0.01,0.02,0])
        x,v=response(M,K,C,t,initial[:2],initial[2:])
        exact=np.array([expm(A*(a-t[0]))@initial for a in t])
        np.testing.assert_allclose(np.c_[x,v],exact,atol=1e-9,rtol=1e-6)
    def test_damping_rejects_small_negative_matrix(self):
        with self.assertRaises(ValueError):
            response([[1]],[[1]],[[-1e-12]],[0,1],[0],[0])
    def test_damping_rejects_asymmetry(self):
        with self.assertRaises(ValueError):
            response(np.eye(2),np.eye(2),[[1,1e-7],[0,1]],[0,1],[0,0],[0,0])
    def test_controls_and_force_validation(self):
        for kwargs in [{'max_step':0},{'rtol':-1},{'atol':np.nan},{'force':[1]}]:
            with self.assertRaises(ValueError): response([[1]],[[1]],[[0]],[0,1],[0],[0],**kwargs)
        with self.assertRaises(ValueError): response([[1]],[[1]],[[0]],[0,1],[0],[0],lambda t: [1,2])
    def test_energy_invalid_shape(self):
        with self.assertRaises(ValueError): energy([[1]],[[1]],[1],[0])
    def test_harmonic_invalid_inputs(self):
        for om,F in [([-1],[1]),([np.nan],[1]),([],[1]),([1],[1,2])]:
            with self.assertRaises(ValueError): harmonic_response([[1]],[[1]],[[0]],om,F)
    def test_sdof_rayleigh(self):
        M,K=chain([2],[800]); C=rayleigh(M,K,0.03)
        self.assertAlmostEqual(C[0,0],sdof(2,800,0.03)[1])

if __name__=='__main__': unittest.main()
