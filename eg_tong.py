"""Normalisation fixed once and for all: Tong's (1/2pi Polyakov; H = 1/g^2 + N/(2r); TN radius g, dual radius 1/g),
HLM's elliptic-genus conventions (eqs. 3.17, 3.39, 4.19).
Tong's A1 ALF GLSM = ALE GLSM with U(1)_F gauged + twisted hyper (Stueckelberg after dualising theta).
At the EG level (HLM 3.15: V_F,Phi_F and Psi,Gamma cancel up to the winding/momentum sum) this is
   E_Tong(g^2) = g^2 int_C dudubar/tau2  E_ALE(zR, z1 + u, z2 - u) e^{-g^2 pi |u|^2/tau2}
with E_ALE from HLM (3.39) N=2 and the shift of the U(1)_F fugacity mapped by (4.19).
Compare with HLM's own A1 Taub-NUT formula (3.17), computed independently in eg_hlm.py."""
import numpy as np, mpmath as mp
from eg_tcp1 import th1
tau = 0.31 + 0.9j; t2 = tau.imag
z, x1, x2 = 0.21 + 0.03j, 0.13 + 0.07j, 0.29 + 0.05j            # HLM's (z, xi1, xi2), same as eg_hlm.py
zR, z1, z2 = z - x2, x1 - x2, -(x1 + x2)                          # (4.19)

def E_ALE(zR, z1, z2):                                            # HLM (3.39), N = 2
    tot = 0
    for a in (0, 1):
        for b in (0, 1):
            c = (a + b*tau)/2
            tot += th1(tau, c - zR + z1)*th1(tau, c + zR - z2)/(th1(tau, c + z1)*th1(tau, c - z2))
    return complex(tot/2)

E_ALE_v = np.vectorize(lambda u: E_ALE(zR, z1 + u, z2 - u))
print("E_ALE (u = 0):", E_ALE(zR, z1, z2))
for g2 in [3, 12, 48]:
    w = np.sqrt(t2/(g2*np.pi))
    rs = (np.arange(160) + 0.5)/160*6*w; ths = 2*np.pi*(np.arange(96) + 0.5)/96
    R, T = np.meshgrid(rs, ths, indexing='ij')
    U = R*np.exp(1j*T)
    E = g2/t2*np.sum(E_ALE_v(U)*np.exp(-(g2*np.pi/t2)*R**2)*R)*(rs[1]-rs[0])*(ths[1]-ths[0])
    print(f"g^2 = {g2:>2}: E_Tong = {E:.8f}")
print("HLM (3.17) values from eg_hlm.py: g^2=3: 1.95856306+0.11207262j, g^2=12: 1.92071966+0.16335005j, g^2=48: 1.87601480+0.23572642j")
