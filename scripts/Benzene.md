---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.8
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Calculate Geometry

```{code-cell} ipython3
from pyscf import gto, dft
mol = gto.M(atom = """
C   1.209     0.698   0.000
C   0.000     1.396   0.000
C  -1.209     0.698   0.000
C  -1.209    -0.698   0.000
C  -0.000    -1.396   0.000
C   1.209    -0.698   0.000
H   2.159     1.247   0.000
H   0.000     2.493   0.000
H  -2.159     1.247   0.000
H  -2.159    -1.247   0.000
H  -0.000    -2.493   0.000
H   2.159    -1.247   0.000
    """,  # in Angstrom
    basis='6-31g',
    symmetry=True)
mf = dft.RKS(mol)
mf.xc ='b3lyp'

from pyscf.geomopt.geometric_solver import optimize

#Run Geomety optimization
mol_eq = optimize(mf, **conv_params)
print(mol)
```

# Mean Field Calculation

Now we setup a DFT ground state calculation using the `pyscf.dft` module. The steps are as follows.
    - we create our mean field object `mf` 
    - choose our exchange correlation functional (in this case we happened to choose `b3lyp`)
    - run the ground state calculation by calling `mf.kernel()`
    
```{code-cell} ipython3
mf = dft.RKS(mol_eq)
mf.xc ='b3lyp'
mf.kernel()
```

# Linear Response/Time Dependent Density Functional Theory (TDDFT)

Finally we use the `pyscf.tddft` module to calculate the excited states of benzene. The output will give us a list of excited state energies as well as their transition dipole moments. 

```{code-cell} ipython3
from pyscf import tddft
mytd = tddft.TDDFT(mf)
mytd.nstates = 15
mytd.kernel()
mytd.analyze(verbose=4)
mytd.transition_dipole()
```

# Relevant Transition  Dipole Moments

+++

# Excited State Definitions
- Singlet excitation energies and oscillator strengths 
- Excited State   1:  ???      1.03262 eV   1200.68 nm  f=nan
- Excited State   2:  B3u      5.57739 eV    222.30 nm  f=0.0000
-      20 -> 23        0.49957
-      21 -> 22        0.50021
- Excited State   3:  ???      6.42869 eV    192.86 nm  f=0.0000
-      20 -> 22       -0.49802
-      21 -> 23        0.49824
- Excited State   4:  ???      7.46446 eV    166.10 nm  f=0.5860
-      20 -> 23       -0.49834
-      21 -> 22        0.49771
- Excited State   5:  ???      7.46501 eV    166.09 nm  f=0.5865
-      20 -> 22        0.49813
-      21 -> 23        0.49791
- Excited State   6:   Au      7.93727 eV    156.21 nm  f=0.0000
-      18 -> 23       -0.48598
-      19 -> 22       -0.51071


## Transition electric dipole moments (AU) 
- state          X           Y           Z        Dip. S.      Osc.
-  2        -0.0010     -0.0000     -0.0000      0.0000      0.0000
-  3         0.0000      0.0001     -0.0000      0.0000      0.0000
-  4        -1.7900     -0.0003      0.0000      3.2043      0.5860
-  5         0.0014      1.7907     -0.0000      3.2066      0.5865
-  6        -0.0000      0.0000      0.0000      0.0000      0.0000
-  7         0.0000     -0.0000     -0.0144      0.0002      0.0000

+++

# Conversion to eV

+++

- Using Hartree Atomic units
- bohr radius = electron charge = 1
- Approx 1.89 Bohr Radius = 1 Angstrom
- Thus divide au electric dip. moments by 1.89 to get eAnstroms

+++

# Reproduce Results Wang 2021, https://doi.org/10.1063/5.0036283 

- __Excitation Energy:__ 6.93 eV
- __Electric Dipole Moment:__ 0.96 eAngstroms

## Our Results
- __Excitation Energy (state 4):__ 7.46446 eV
- __Electric Dipole Moment (state 4):__ 1.7900 eBohr = 0.94 eAngstroms

+++

# Plug DFT Results into QuTip

+++

# See here for conversions and constants https://physics.nist.gov/cuu/Constants/index.html
- 1 Hartree in eV = 27.211 386 245 988 eV 
- 1 au electric permitivity = $ e^2/(a_0 E_h) = 1.11 *10^{-10} F m-1$ 
- vacuum permitivity in SI = $8.85 * 10^{-12} F m-1$
- vacuum permitivity in au = $\frac{8.85 * 10^{-12} F m-1}{1.11 *10^{-10} F m-1} [e^2/(a_0 E_h)]$
- vacuum permitivity in au = $7.97* 10^{-2} [e^2/(a_0 E_h)]$
- Hartree Energy $E_h$ in eV = 27.211 386 245 988 eV
- Excitation energy = 7.46 eV = 7.46/(27.2) E_h
- Bohr radius a_0 in SI = 5.291 772 109 03 x 10-11 m ~ 5.291 * 10^(-2) nm
- 1 au volume in nm^3  = a_0^3 = 148.12 * 10^(-6) nm^3 ~ 1.48 * 10^(-4) nm^3
- 10^4 nm^3 in au = 10^(4) nm^3 * (a_0^3)/(1.48 * 10^(-4) nm^3) = 1/(1.48) a_0^3 ~ .68 a_0^3
- Wang 2021 uses a cavity strength parameter$\vec{\lambda}_{k}$
- $\vec{\lambda}_{k}\equiv\sqrt{\frac{2}{\hbar \omega_k}}\vec{E}_{k} e$
- Wang 2020 defines coupling stregth $g_{ik} = -\frac{e \vec{E}_k \cdot \langle g | \vec{R} | e_i\rangle}{\hbar} =  -\sqrt{\frac{\omega_k}{2\hbar}} \frac{ \vec{\lambda}_k}{e} \cdot \vec{d}_i$
- Wang 2020 uses $\vec{\lambda}_c = (.001,.002,.008) \frac{eV^{1/2}}{nm}$
- Convert cavity strength to au: $\frac{eV^{1/2}}{nm} \cdot \sqrt{\frac{1 E_h }{27.211 eV}}\cdot \frac{ 5.291 * 10^(-2) nm}{a_0}$
- 1 unit of $\frac{eV^{1/2}}{nm} ~ .0102 [E_h^{1/2} /nm]$

+++

## $g = -\frac{e\langle g|x| 4 \rangle \mathcal{E}_{\vec{k}_x} }{\hbar}$
#### $g = 1.79 [e a_0] * 1.6 [\frac{E_h}{\hbar e a_0}]$
#### $g ~ 2.8 \frac{E_h}{\hbar}$

+++

# Conversion Code (use Wang 2020 vals and try to be consistent)

```{code-cell} ipython3
#unit constants
import numpy as np
auDictionary = { 'unit_energy':{ 'eV':27.211 },
                 'unit_permitivity':{'SI':1.11265005545*10**(-10) },
                'unit_length' : { 'SI': 5.29177210903*10**(-11),
                                 'A':0.529177,
                                 'nm':.0529177210903},
                'unit_cavity_strength': { 'root_eV_nm':
                                         np.sqrt(1/27.211)/(5.29177210903*10**(-2))}
               }
constantsDictionary = { 'vacuum_permitivity' : {'SI' :  8.8541878128 * 10**(-12),
                                                'au' : 1.11265005545* 10**(-10)}
                        
                      }

wang2020Dictionary = { 'transition_dipole':
                      {'eA':0.96,
                       'enm' : 0.096,
                       'au':0.96/auDictionary['unit_length']['A']
                      },
                      'cavity_strengths': 
                      { 'root_eV_nm' : [.001,.002,.008],
                        'root_Eh_a0' : [ 
                            x*auDictionary['unit_cavity_strength']['root_eV_nm']**(-1)
                            for x in [.001,.002,.008] 
                        ]
                                          
                      },
                      'excited_state_energy':
                      { 'eV': 6.93,
                        'au': 6.93/auDictionary['unit_energy']['eV']
                      },
                      'loss_rates':
                      { 'eV': [.001,.004,.008],
                        'hartrees': [ x/auDictionary['unit_energy']['eV'] for x in [.001,.004,.008]]
                      }
                     }
                       
def hbar_g_wang2020(hbar_w_cav,lambda_cav,transition_dipole):
                       #w_cav = (1/hbar) hbar_w_cav
                       #hbar * g = -hbar*sqrt(w_cav/(2hbar)) lmbda/e *d
                       #         = -hbar*sqrt(hbar_w_cav/(2 hbar**2)) * lmbda/e *d
                       #        = -sqrt(hbar_w_cav/2) * lmda/e * d          
                       result = -np.sqrt(hbar_w_cav/2)* lambda_cav*transition_dipole #e=1
                       return result
```

```{code-cell} ipython3
hbar_g_wang2020(6.93,.001,0.096)/.001
```

```{code-cell} ipython3
from qutip import *
import numpy as np
N = 15
wa = wang2020Dictionary['excited_state_energy']['eV']
wc = wa
cavity_strength = wang2020Dictionary['cavity_strengths']['root_eV_nm'][2]
transition_dipole = wang2020Dictionary['transition_dipole']['enm']
g = hbar_g_wang2020(wa,cavity_strength,transition_dipole)
kappa = wang2020Dictionary['loss_rates']['eV'][2]



# Jaynes-Cummings Hamiltonian
a  = tensor(destroy(N), qeye(2))
sm = tensor(qeye(N), destroy(2))
H = wc * a.dag() * a + wa * sm.dag() * sm + g * (a.dag() * sm + a * sm.dag())

# collapse operators
n_th = 0.25
c_ops = [np.sqrt(kappa * (1 + n_th)) * a, np.sqrt(kappa * n_th) * a.dag()]


wlist = np.linspace(wc-1, wc+1, 10000)

spec = qutip.spectrum(H, wlist, c_ops, a.dag(), a)
```

```{code-cell} ipython3
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 1)


ax.plot(wlist -wc ,spec,label="qutip_eval")

plt.axvline(x=-g,color='r',ls='--',label="poles (-g,g)")
plt.axvline(x=g,color='r',ls='--')



plt.title('Wang2021 Parameters')
plt.xlabel('eV')
plt.ylabel('Intensity')
ax.set_xlim([-.02, .02])

ax.legend()

plt.show()
```

```{code-cell} ipython3
(8.8541878128 * 10**(-12) )*(4*np.pi)
```

# Application of DFT and TDDFT to Jaynes-Cummings Hamiltonian

```{code-cell} ipython3

```

+++ {"tags": []}

### $H_{Dipole} = -e \vec{r} \cdot \vec{E}= -e\langle g|x| 4 \rangle(\sigma_{g4} + \sigma_{4g}) \mathcal{E}_{kx}(a + a^\dagger)$
- $\mathcal{E}_{\vec{k}_x} = (\hbar \nu_k /( 2 \varepsilon_0 V))^{1/2}$ (Scully, ch. 6)
- $\mathcal{E}_{\vec{k}_x} = \left ( \frac{7.46/27.2 [E_h]}{2 (7.97* 10^{-2} V [e^2/(a_0 E_h)]}\right)^{1/2}$
- $\mathcal{E}_{\vec{k}_x} = \left ( \frac{.274 [E_h^2] }{1.594 * 10^{-1} V [e^2 /a_0]}\right)^{1/2}$
- Mode Volume in May 2020 seems to be of order 10^4nm^3 so let $V=10^4 nm^3 = .68 [a_0^3] $
- $\mathcal{E}_{\vec{k}_x} = \left ( \right \frac{.274 [E_h^2]}{.108 [e^2 a_0^2]} )^{1/2}$
- $\mathcal{E}_{\vec{k}_x} = \sqrt{.274/.108} [\frac{E_h}{e a_0}]$
- $\mathcal{E}_{\vec{k}_x} ~ 1.6[\frac{E_h}{e a_0}] $
can we get 70meV coupling?

```{code-cell} ipython3
def scully_hbar_g(hbar_w_cav,vacuum_permitivity,cavity_volume,transition_dipole):
    efield_strength = np.sqrt(hbar_w_cav/(2.0*vacuum_permitivity*cavity_volume))
    return -transition_dipole*efield_strength
                              
```

```{code-cell} ipython3
#using au
wa = 7.46446 
wc = wa
transition_dipole = .094
e_0 = (18.9*27.2)* constantsDictionary['vacuum_permitivity']['au'] 
cav_vol = 20**3/(auDictionary['unit_length']['nm']**3)
g = scully_hbar_g(wa,e_0,cav_vol,-transition_dipole)
kappa = .06
```

```{code-cell} ipython3

# Jaynes-Cummings Hamiltonian
a  = tensor(destroy(N), qeye(2))
sm = tensor(qeye(N), destroy(2))
H = wc * a.dag() * a + wa * sm.dag() * sm + g * (a.dag() * sm + a * sm.dag())

# collapse operators
n_th = 0.25
c_ops = [np.sqrt(kappa * (1 + n_th)) * a, np.sqrt(kappa * n_th) * a.dag()]


wlist = np.linspace(wc-3, wc +3, 1000)

spec = qutip.spectrum(H, wlist, c_ops, a.dag(), a)
```

```{code-cell} ipython3
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 1)


ax.plot(wlist -wc ,spec,label="qutip_eval")

plt.axvline(x=-g,color='r',ls='--',label="poles (-g,g)")
plt.axvline(x=g,color='r',ls='--')



plt.title('Scully + DFT')
plt.xlabel('eV')
plt.ylabel('Intensity')
ax.set_xlim([-.5, .5])

ax.legend()

plt.show()
```

```{code-cell} ipython3
g
```
