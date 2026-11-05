import unittest
import numpy as np
from vibrastruct import *

class Validation(unittest.TestCase):
    def test_sdof_exact_decay(self):
        m,k,z=2.,800.,0.04; wn,c=sdof(m,k,z); t=np.linspace(0,2,2001)
        x,v=response([[m]],[[k]],[[c]],t,[0.01],[0.])
        wd=wn*np.sqrt(1-z*z)
        exact=0.01*np.exp(-z*wn*t)*(np.cos(wd*t)+z*wn/wd*np.sin(wd*t))
        np.testing.assert_allclose(x[:,0],exact,atol=1e-9,rtol=1e-6)
    def test_two_mass_closed_form(self):
        M,K=chain([2,2],[800,800]); w,P=modes(M,K)
        np.testing.assert_allclose(w*w,400*np.array([(3-np.sqrt(5))/2,(3+np.sqrt(5))/2]),rtol=1e-12)
        np.testing.assert_allclose(P.T@M@P,np.eye(2),atol=1e-12)
    def test_modal_residual(self):
        M,K=chain([250,180,120],[1e6,8e5,6e5]); w,P=modes(M,K)
        np.testing.assert_allclose(K@P,M@P@np.diag(w*w),atol=1e-9)
    def test_static_limit(self):
        self.assertAlmostEqual(receptance(2,800,0.1,0).real,1/800)
    def test_stiffness_scaling(self):
        M,K=chain([2,3],[800,500]); w,_=modes(M,K); w2,_=modes(M,2*K)
        np.testing.assert_allclose(w2,np.sqrt(2)*w)
    def test_energy(self):
        M,K=chain([2,3],[800,500]); t=np.linspace(0,2,2001)
        x,v=response(M,K,np.zeros_like(M),t,[0,0.01],[0,0]); E=energy(M,K,x,v)
        self.assertLess(np.max(np.abs(E/E[0]-1)),1e-7)
        x,v=response(M,K,rayleigh(M,K),t,[0,0.01],[0,0]); E=energy(M,K,x,v)
        self.assertLess(np.max(np.diff(E)),1e-10); self.assertLess(E[-1],E[0])
    def test_harmonic_time_domain(self):
        m,k,z=2,800,0.1; wn,c=sdof(m,k,z); om=12.; H=receptance(m,k,z,om)*5
        t=np.linspace(0,1,1001)
        x,v=response([[m]],[[k]],[[c]],t,[H.real],[-om*H.imag],lambda s: [5*np.cos(om*s)])
        np.testing.assert_allclose(x[:,0],np.real(H*np.exp(1j*om*t)),atol=1e-9)
    def test_rayleigh_targets(self):
        M,K=chain([2,3,4],[800,500,600]); w,P=modes(M,K); C=rayleigh(M,K,0.03)
        z=np.diag(P.T@C@P)/(2*w)
        np.testing.assert_allclose(z[[0,-1]],0.03)
    def test_invalid_inputs(self):
        for args in [(-1,800,0.1),(1,0,0.1),(1,800,-1),(np.nan,1,0)]:
            with self.assertRaises(ValueError): sdof(*args)
        with self.assertRaises(ValueError): chain([1],[1,2])
        with self.assertRaises(ValueError): modes([[1]], [[-1]])
        with self.assertRaises(ValueError): response([[1]],[[1]],[[0]],[1,0],[0],[0])

if __name__=='__main__': unittest.main()
