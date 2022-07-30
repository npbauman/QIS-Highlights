This repository is dedicated to highlighting developments within the Quantum Information Science (QIS) Initiative. It also provides some links, software, and data pertinent to the projects as well. 

**Table of Contents**
- [Development of the NWChem-QDK interface](#Development-of-the-NWChem-QDK-interface)
- [Downfolding techniques for dimension reduction of correlated electronic Hamiltonians](#Downfolding-techniques-for-dimension-reduction-of-correlated-electronic-Hamiltonians)
- [Quantum simulations of excited state with active-space downfolded Hamiltonians](#Quantum-simulations-of-excited-state-with-active-space-downfolded-Hamiltonians)
- [Resource-efficient VQE algorithms with downfolded Hamiltonians](#Resource-efficient-VQE-algorithms-with-downfolded-Hamiltonians)
- [Coupled Cluster Green’s function formulations based on the effective Hamiltonians](#coupled-cluster-greens-function-formulations-based-on-the-effective-hamiltonians)
- [A Filon-like integration strategy for calculating exact exchange in periodic boundary conditions](#a-filon-like-integration-strategy-for-calculating-exact-exchange-in-periodic-boundary-conditions)
- [Sub-system quantum dynamics using coupled cluster downfolding techniques](#Sub-system-quantum-dynamics-using-coupled-cluster-downfolding-techniques)
- [Exploiting chemistry and molecular systems for quantum information science](#Exploiting-chemistry-and-molecular-systems-for-quantum-information-science)
- [Quantum computing for high-energy states - QPE simulations of core-level states](#Quantum-computing-for-high-energy-states---QPE-simulations-of-core-level-states)
- [Quantum algorithms for simulating the lattice Schwinger Model](#Quantum-Algorithms-for-Simulating-the-Lattice-Schwinger-Model)
- [Quantum algorithms for connected moments expansion](#Quantum-algorithms-for-connected-moments-expansion)
- [Variational Quantum Eigensolver for Approximate Diagonalization of Downfolded Hamiltonians using Generalized Unitary Coupled Cluster Ansatz](#Variational-Quantum-Eigensolver-for-Approximate-Diagonalization-of-Downfolded-Hamiltonians-using-Generalized-Unitary-Coupled-Cluster-Ansatz)
- [Reduced-size plane-wave based representation of many-body methods for quantum computing](#Reduced-size-plane-wave-based-representation-of-many-body-methods-for-quantum-computing)
- [Benchmarking adaptive variational quantum eigensolvers](#Benchmarking-adaptive-variational-quantum-eigensolvers)
- [Even more efficient quantum computations of chemistry through tensor hypercontractions](#Even-more-efficient-quantum-computations-of-chemistry-through-tensor-hypercontractions)
- [Dynamical self-energy mapping for quantum computing](#dynamical-self-energy-mapping-for-quantum-computing)
- [Simulating quantum materials with digital quantum computers](#Simulating-quantum-materials-with-digital-quantum-computers)
- [Variational quantum solver employing the PDS energy functional](#Variational-quantum-solver-employing-the-PDS-energy-functional)
- [Dimensionality reduction of many-body problem using coupled-cluster sub-system flow equations](#dimensionality-reduction-of-many-body-problem-using-coupled-cluster-sub-system-flow-equations)
- [Improving the accuracy and efficiency of quantum connected moments expansions](#improving-the-accuracy-and-efficiency-of-quantum-connected-moments-expansions)
- [Numerical Simulations of Noisy Quantum Circuits for Computational Chemistry](#numerical-simulations-of-noisy-quantum-circuits-for-computational-chemistry)
- [Coupled cluster downfolding methods: The effect of double commutator terms on the accuracy of ground-state energies](#coupled-cluster-downfolding-methods-the-effect-of-double-commutator-terms-on-the-accuracy-of-ground-state-energies)




## Development of the NWChem-QDK interface
<p align="center">
  <img width="500" src="https://github.com/npbauman/BES-QIS-at-PNNL/blob/a14900427b39cae1ccb64c3d2272b8d2782e5465/Figures/QuantumWorkflow.png"> <br>
<sub><sup>Detailed workflow for simulating quantum chemistry on quantum computers.</sup></sub>
</p>
  
"Q# and NWChem: Tools for Scalable Quantum Chemistry on Quantum Computers," Guang Hao Low, Nicholas P. Bauman, Christopher E. Granade, Bo Peng, Nathan Wiebe, Eric J. Bylaska, Dave Wecker, Sriram Krishnamoorthy, Martin Roetteler, Karol Kowalski, Matthias Troyer and Nathan A. Baker, [*arXiv*:1904.01131](https://arxiv.org/abs/1904.01131).

**Challenge:** Enabling new quantum computing simulations in chemistry requires tight integration of computational chemistry infrastructure with novel quantum algorithms for ground- and excited-state simulations.

**Approach and Results:** Our team has developed NWChem-QDK (Quantum Development Kit (QDK) developed and maintained by Microsoft Research group) interface to perform Quantum Phase Estimator (QPE) simulations for chemical processes. NWChem provides a reach infrastructure for the characterization of second quantized forms of electronic Hamiltonians and the initial characterization of ground- and excited-state wavefunction obtained with high-accuracy coupled-cluster calculations on classical systems. QDK offers a variety of quantum simulation algorithms ranging from Trotter-Suzuki expansion to various qubitization approaches. Quantum simulations using NWChem and QDK can also be performed employing the web interface EMSL Arrows Quantum Editor.

**Significance and Impact:** Using the NWChem-QDK interface, we demonstrated the efficiency of the QPE approach in describing strongly correlated ground and excited states of molecular systems. The QPE algorithm, with proper initial estimates of the electronic wave functions, was able to describe potential energy surfaces for the ground state and several low-lying excited states, which also involved challenging doubly excited states. Using the H<sub>10</sub> benchmark system in the STO-3G basis set, we demonstrated the advantages of using QPE in achieving highly accurate energy estimates in the strongly correlated regime.

https://nwchemgit.github.io/ (NWChem Documentation)  
https://arrows.emsl.pnnl.gov/api/qsharp_chem (EMSL Arrows Quantum Editor)  
https://nwchemgit.github.io/EMSL_Arrows.html (EMSL Arrows Documentation)  
https://docs.microsoft.com/en-us/azure/quantum/ (Quantum Development Kit Documentation)  



## Downfolding techniques for dimension reduction of correlated electronic Hamiltonians
<p align="center">
  <img width="400" src="https://github.com/npbauman/BES-QIS-at-PNNL/blob/d71bfd80cfc7c194d80bb1fabb1315a6a853780d/Figures/DownfoldingAbstract.jpg"> <br>
<sub><sup>A schematic representation of the coupled cluster (CC) downfolding technique, where the CC Ansatz is used to generate effective/downfolded Hamiltonians in small-dimensionality active spaces.</sup></sub>
</p>
  
"Downfolding of many-body Hamiltonians using active-space models: Extension of the sub-system embedding sub-algebras approach to unitary coupled cluster formalisms," Nicholas P. Bauman, Eric J. Bylaska, Sriram Krishnamoorthy, Guang Hao Low, Nathan Wiebe, Christopher E. Granade, Martin Roetteler, Matthias Troyer and Karol Kowalski, [*J. Chem. Phys.* **151**, 014107 (2019)](https://doi.org/10.1063/1.5094643); [*arXiv*:1902.01553](https://arxiv.org/abs/1902.01553).

**Challenge:** Limited quantum resources preclude simulations of complex and realistic chemical processes. Therefore, quantum computing is in high demand for efficient techniques for re-representing quantum many-body problems in reduced dimensionality spaces before reaching maturity.

**Approach and Results:** We have extended the sub-system embedding sub-algebras coupled cluster (SES-CC) theory to the downfolding procedure based on the double unitary CC (DUCC) formalism to address this challenge. In contrast to the standard single-reference SES-CC formulations, the DUCC approach results in a Hermitian form of the effective Hamiltonian in an active-space, which provides a rigorous separation of external cluster amplitudes that describe dynamical correlation effects from those corresponding to the internal (within the active space) excitations that define the components of eigenvectors associated with the energy of the entire system. 

**Significance and Impact:** We have extended the sub-system embedding sub-algebras coupled cluster (SES-CC) theory to the downfolding procedure based on the double unitary CC (DUCC) formalism to address this challenge. In contrast to the standard single-reference SES-CC formulations, the DUCC approach results in a Hermitian form of the effective Hamiltonian in an active-space, which provides a rigorous separation of external cluster amplitudes that describe dynamical correlation effects from those corresponding to the internal (within the active space) excitations that define the components of eigenvectors associated with the energy of the entire system.

https://nwchemgit.github.io/ (NWChem Documentation)  
https://docs.microsoft.com/en-us/azure/quantum/ (Quantum Development Kit Documentation)





##  Quantum simulations of excited state with active-space downfolded Hamiltonians
<p align="center">
  <img width="400" src="https://github.com/npbauman/BES-QIS/blob/fa24a1e79352f4a8084548514c12154ae1e2207b/Figures/Excited_States_H2.png"> <br>
<sub><sup> A typical distribution of energies for a strongly correlated system obtained from several simulations with the QPE algorithm. This particular spread of energies corresponds to H<sub>2</sub>.</sup></sub>
</p>

"Quantum simulations of excited state with active-space downfolded Hamiltonians," Nicholas P. Bauman, Guang Hao Low and Karol Kowalski, [*J. Chem. Phys.* **151**, 234114 (2019)](https://doi.org/10.1063/1.5128103); [*arXiv*:1909.06404](https://arxiv.org/abs/1909.06404).

**Challenge:** Even though impressive progress has been achieved in developing wavefunction-based excited-state approaches, problems with the description of complicated states dominated by high-rank excitations still exist. Of particular interest is the application to strongly correlated molecular systems characterized by small energy gaps between occupied and unoccupied orbitals where multiple electronic states defined by complex collective phenomena lie.

**Approach and Results:** The stochastic nature of the quantum phase estimation (QPE) algorithm opens the opportunity to discover/identify novel and exotic excited states that classical methods struggle in describing. The probability of finding such states through repeated simulations is driven by an initial guess of the true configuration structure of a sought-after excited state, even if it is limited or postulated. We demonstrated that not only does the QPE capture several excited states of various complexity, but a distribution of eigenstate energies is collected from a single initial state through repeated QPE simulations. This contrasts with other quantum algorithms, such as VQE approaches, which only provide energy estimates for a single targeted state. The QPE algorithm was combined with the excited-state extension of the double unitary coupled cluster (DUCC) formalism to construct active-space representations of downfolded Hamiltonians that can be used to reproduce a large portion of excited-state correlation effect and are amenable for quantum computers.

**Significance and Impact:** It was demonstrated that the QPE algorithm is an efficient tool for testing various "excited-state" hypotheses for strongly correlated systems, which usually pose a significant challenge for existing many-body formalisms. Using QPE techniques, one can obtain a spectrum of all states that have non-negligible overlap with the hypothesis state, a fact is well-known in quantum computing but deserves broader exposure. The stochastic character of QPE may be instrumental in studies of excited-state processes, especially in strongly correlated or metallic-like systems. In addition, the promising results of the DUCC excited-state extension pave the way for future developments of the method.

https://nwchemgit.github.io/ (NWChem Documentation)  
https://docs.microsoft.com/en-us/azure/quantum/ (Quantum Development Kit Documentation)






## Resource-efficient VQE algorithms with downfolded Hamiltonians
<p align="center">
  <img width="400" src="https://github.com/npbauman/BES-QIS/blob/cfa60174e5876ea633879b4145b293318583cb85/Figures/VQE-DUCC_Li2.png"> <br>
<sub><sup> Ground state energy of Li<sub>2</sub> with the cc-pvtz basis. The errors, ε, of the VQE-DUCC and the active-space CCSDTQ calculations are with respect to the 60-orbital CCSDTQ energies.</sup></sub>
</p>

"Resource-Efficient Chemistry on Quantum Computers with the Variational Quantum Eigensolver and the Double Unitary Coupled-Cluster Approach," Mekena Metcalf, Nicholas P. Bauman, Karol Kowalski and Wibe A. de Jong, [*J. Chem. Theory Comput.* **16**, 6165 (2020)](https://doi.org/10.1021/acs.jctc.0c00421); [*arXiv*:2004.07721](https://arxiv.org/abs/2004.07721).

**Challenge:** In modeling many-body problems, the biggest challenge confronted is that the number of qubits scales linearly with the molecular basis's size. This poses a significant limitation on the basis sets' size and the number of correlated electrons included in quantum simulations of chemical processes.

**Approach and Results:** To address this issue and enable more realistic simulations on NISQ computers, we employed the double unitary coupled-cluster method to effectively downfold correlation effects into the reduced-size orbital space, commonly referred to as the active space. Using downfolding and VQE techniques, we demonstrated that effective Hamiltonians could capture the effect of the whole orbital space in small-size active spaces, especially when natural orbitals are employed. 

**Significance and Impact:** The DUCC/VQE framework was used to solve the ground-state energy of H<sub>2</sub>, Li<sub>2</sub>, and BeH<sub>2</sub> on the cc-pVTZ basis using the reduced-size active spaces. The VQE formalism has also been extended to the generalized unitary CC (GUCC) ansatz and applied to benchmark systems/processes (N<sub>2</sub>, H<sub>2</sub>O, and C<sub>2</sub>H<sub>4</sub>) described by the downfolded Hamiltonians. The preliminary data indicate that simple downfolding procedures based on single commutator expansion and small active spaces can recover around 90% of correlation energy calculated when all orbitals are correlated.

https://nwchemgit.github.io/ (NWChem Documentation)  
https://qiskit.org/ (Qiskit)





## Coupled Cluster Green’s function formulations based on the effective Hamiltonians
<p align="center">
  <img width="400" src="https://github.com/npbauman/BES-QIS/blob/f1a9e3e87afdff033662f8c1f3d48a02be7940ee/Figures/DUCC-GFCCSD_N2.png"> <br>
<sub><sup> Spectral functions of N<sub>2</sub> in the valence energy regimes directly computed by the closed-shell GFCCSD and DUCC-GFCCSD methods with cc-pVDZ basis set.</sup></sub>
</p>

"Coupled Cluster Green's function formulations based on the effective Hamiltonians," Nicholas P. Bauman, Bo Peng and  Karol Kowalski, [*Mol. Phys.* **118** (2020)](https://doi.org/10.1080/00268976.2020.1725669); [*arXiv*:1910.00394](https://arxiv.org/abs/1910.00394).
  
**Challenge:** The recently developed double unitary coupled cluster (DUCC) formalism for for integrating out high-energy wave-function components from low-energy ones in the effective (or downfolded) Hamiltonians was shown to be promising. While in its prime, the utility of the formalism still needed affirming by combining it with existing and powerful many-body methods.

**Approach and Results:** The Green’s function coupled cluster (GFCC) approach remains an active area of development and was chosen to investigate the performance of the DUCC downfolding proceedure. This combined approach (DUCC-GFCC) was applied to H<sub>2</sub>O, N<sub>2</sub>, CO, and trans-1,3-butadiene and shown to provide a significant reduction of numerical effort and good agreement with the corresponding all-orbital GFCC methods in energy windows that are consistent with the choice of active space. Spectral functions were also shown to systematically improve with larger active spaces while maintaining a significant reduction in dimensionality, unlike standard truncated Hamiltonians where errors were sporadic. 

**Significance and Impact:** We demonstrated that the utilization of the effective Hamiltonian stemming from the DUCC downfolding procedure can be used to reproduce the main features of the standard GFCCSD spectral function while significantly reducing the cost of the GFCC calculations. A growing interest in quantum computing algorithms for correlated Green's function makes reduced-dimension DUCC-GFCC formulations a possible target for early quantum computing applications.

https://nwchemgit.github.io/ (NWChem Documentation)  





## A Filon-like integration strategy for calculating exact exchange in periodic boundary conditions

"A Filon-like integration strategy for calculating exact exchange in periodic boundary conditions: a plane-wave DFT implementation," Eric J. Bylaska, Kevin Waters, Eric D. Hermes, Judit Zádor and Kevin M Rosso, [*Mater. Theory* **4**, 3 (2020)](https://doi.org/10.1186/s41313-020-00019-9).

**Abstract:** An efficient and accurate approach for calculating exact exchange and other two-electron integrals has been developed for periodic electronic structure methods. Traditional approaches used for integrating over the Brillouin zone in band structure calculations, e.g. trapezoidal or Monkhorst-Pack, are not accurate enough for two-electron integrals. This is because their integrands contain multiple singularities over the double integration of the Brillouin zone, which with simple integration methods lead to very inaccurate results. A common approach to this problem has been to replace the Coulomb interaction with a screened Coulomb interaction that removes singularities from the integrands in the two-electron integrals, albeit at the inelegance of having to introduce a screening factor which must precomputed or guessed. Instead of introducing screened Coulomb interactions in an ad hoc way, the method developed in this work derives an effective screened potential using a Filon-like integration approach that is based only on the lattice parameters. This approach overcomes the limitations of traditionally defined screened Coulomb interactions for calculating two-electron integrals, and makes chemistry many-body calculations tractable in periodic boundary conditions. This method has been applied to several systems for which conventional DFT methods do not work well, including the reaction pathways for the addition of H<sub>2</sub> to phenol and Au<sub>20</sub><sup>-</sup> nanoparticle, and the electron transfer of a charge trapped state in the Fe(II) containing mica, annite.


<!---**Challenge:**--->

<!---**Approach and Results:**--->

<!---**Significance and Impact:**--->






## Sub-system quantum dynamics using coupled cluster downfolding techniques

"Sub-system quantum dynamics using coupled cluster downfolding techniques," Karol Kowalski and Nicholas P. Bauman, [*J. Chem. Phys.* **152**, 244127 (2020)]( https://doi.org/10.1063/5.0008436); [*arXiv*:2003.09566](https://arxiv.org/abs/2003.09566).
  
**Challenge:** As quantum computing rapidly develops, there is a growing interest in simulations involving imaginary time evolution. However, like their standard counterparts, the scaling of the number of qubits with the system size imposes a significant limitation on the basis set size and number of correlated electrons. 

**Approach and Results:** We presented an extension of the sub-system embedding sub-algebra coupled cluster (SESCC) formalism and the double unitary coupled cluster (DUCC) Ansatz to the time domain. Prior to this work, development of the DUCC formalism was predicated on assuming the exactness of the double unitary CC ansatz. We were able to complement the time domain discussion by showing that the exact wave function can be represented by the double unitary exponential ansatz with general-type anti-Hermitian many-body cluster operators representing internal and external excitations. This result corresponds to the general property of the exact wave function proven at the level of SESCC formalism.

**Significance and Impact:** The DUCC formalism provides a rigorous many-body characterization of the time-dependent action functional to describe the dynamics of the entire system in time modes captured by the corresponding active space. This approach and corresponding approximations can not only reduce the cost of time-domain CC simulations for larger molecular applications. However, they can also be employed in the imaginary time evolution, which has recently been intensively studied in the context of quantum computing.





## Exploiting chemistry and molecular systems for quantum information science
<p align="center">
  <img width="600" src="https://github.com/npbauman/BES-QIS/blob/655788a273ab82d0cae05acc2564874715a15491/Figures/NatRev_Figure.png"> <br>
</p>

"Exploiting chemistry and molecular systems for quantum information science," Michael R. Wasielewski, Malcolm D. E. Forbes, Natia L. Frank, Karol Kowalski, Gregory D. Scholes, Joel Yuen-Zhou, Marc A. Baldo, Danna E. Freedman, Randall H. Goldsmith, Theodore Goodson III, Martin L. Kirk, James K. McCusker, Jennifer P. Ogilvie, David A. Shultz, Stefan Stoll and K. Birgitta Whaley , [*Nat. Rev. Chem.* **2**, 490 (2020)](https://doi.org/10.1038/s41570-020-0200-5).
  
**Abstract:** The power of chemistry to prepare new molecules and materials has driven the quest for new approaches to solve problems having global societal impact, such as in renewable energy, healthcare and information science. In the latter case, the intrinsic quantum nature of the electronic, nuclear and spin degrees of freedom in molecules offers intriguing new possibilities to advance the emerging field of quantum information science. In this Perspective, which resulted from discussions by the co-authors at a US Department of Energy workshop held in November 2018, we discuss how chemical systems and reactions can impact quantum computing, communication and sensing. Hierarchical molecular design and synthesis, from small molecules to supramolecular assemblies, combined with new spectroscopic probes of quantum coherence and theoretical modelling of complex systems, offer a broad range of possibilities to realize practical quantum information science applications.





## Quantum computing for high-energy states - QPE simulations of core-level states
<p align="center">
  <img width="600" src="https://github.com/npbauman/BES-QIS/blob/812cd57d73f44015f805fcf310805d62db7cba23/Figures/Core-Level.png"> <br>
<sub><sup> QPE simulations of doubly excited core-level states of H<sub>2</sub>O in cc-pVDZ basis set using initial guess representing singlet combination of Slater determinants corresponding to simultaneous excitation of electrons from core and valence levels to virtual valence level.</sup></sub>
</p>

"Toward Quantum Computing for High-Energy Excited States in Molecular Systems: Quantum Phase Estimations of Core-Level States," Nicholas P. Bauman, Hongbin Liu, Eric J. Bylaska, Sriram Krishnamoorthy, Guang Hao Low, Christopher E. Granade, Nathan Wiebe, Nathan Baker, Bo Peng, Martin Roetteler, Matthias Troyer and Karol Kowalski, [*J. Chem. Theory Comput.* **17**, 201 (2021)](https://doi.org/10.1021/acs.jctc.0c00909); [*arXiv*:2007.06185](https://arxiv.org/abs/2007.06185).

**Challenge:** The experimental effort at the BES-supported light source facilities is contingent on the availability of predictive modeling tools to interpret x-ray experimental spectra. However, the theoretical description and numerical identification of complicated core-level states still pose a significant challenge for approximate methods.

**Approach and Results:** The stochastic nature of the QPE algorithm can facilitate the discovery/identification of core-level states when the knowledge about the true configuration structure of a sought-after excited state is limited or postulated. For this purpose, we developed an algorithm where through repeated simulations, one can accumulate samples from the distribution of eigenstates energies. The desired error in each energy estimate is inversely proportional to the number of applications of the time evolution operator in the QPE algorithm. This is in contrast to VQE approaches, which only provide energy estimates for a single targeted state. 

**Significance and Impact:** The PNNL-Microsoft team has demonstrated that extension of the QPE algorithm to high-energy core-level states of various spin, spatial symmetries, and complexity is possible. This is the first demonstration of quantum algorithms' potential for identifying shake-up/satellite states and their future role in supporting various x-ray spectroscopies.

https://nwchemgit.github.io/ (NWChem Documentation)  
https://docs.microsoft.com/en-us/azure/quantum/ (Quantum Development Kit Documentation)  




## Quantum algorithms for simulating the lattice Schwinger Model

"Quantum Algorithms for Simulating the Lattice Schwinger Model," Alexander F. Shaw, Pavel Lougovski, Jesse R. Stryker and Nathan Wiebe, [*Quantum* **4**, 306 (2020)](	https://doi.org/10.22331/q-2020-08-10-306); [*arXiv*:2002.11146](https://arxiv.org/abs/2002.11146).

**Abstract:** The Schwinger model (quantum electrodynamics in 1+1 dimensions) is a testbed for the study of quantum gauge field theories. We give scalable, explicit digital quantum algorithms to simulate the lattice Schwinger model in both NISQ and fault-tolerant settings. In particular, we perform a tight analysis of low-order Trotter formula simulations of the Schwinger model, using recently derived commutator bounds, and give upper bounds on the resources needed for simulations in both scenarios. In lattice units, we find a Schwinger model on **N/2** physical sites with coupling constant **x<sup>−1/2</sup>** and electric field cutoff **x<sup>−1/2</sup>Λ** can be simulated on a quantum computer for time **2xT** using a number of **T**-gates or CNOTs in **O(N<sup>3/2</sup>T<sup>3/2</sup>Λ√x)** for fixed operator error. This scaling with the truncation **Λ** is better than that expected from algorithms such as qubitization or QDRIFT. Furthermore, we give scalable measurement schemes and algorithms to estimate observables which we cost in both the NISQ and fault-tolerant settings by assuming a simple target observable–the mean pair density. Finally, we bound the root-mean-square error in estimating this observable via simulation as a function of the diamond distance between the ideal and actual CNOT channels. This work provides a rigorous analysis of simulating the Schwinger model, while also providing benchmarks against which subsequent simulation algorithms can be tested.

<!---**Challenge:**--->

<!---**Approach and Results:**--->

<!---**Significance and Impact:**--->




## Quantum algorithms for connected moments expansion
<p align="center">
  <img width="400" src="https://github.com/npbauman/BES-QIS/blob/812cd57d73f44015f805fcf310805d62db7cba23/Figures/CMXCircuit.png"> <br>
<sub><sup> A schematic representation of the circuit depth reduction in the quantum CMX algorithms.</sup></sub>
</p>

"Quantum simulations employing connected moments expansion," Karol Kowalski and Bo Peng, [*J. Chem. Phys.* **153**, 201102 (2020)](https://doi.org/10.1063/5.0030688); [*arXiv*:2009.05709](https://arxiv.org/abs/2009.05709).

**Challenge:** Further advancement of quantum computing is contingent on enabling models that avoid deep circuits and the excessive use of CNOT gates and associated noise effects.

**Approach and Results:** We developed a quantum algorithm employing finite-order connected moment expansions (CMX) and affordable procedures for initial-state preparation to address this problem. We demonstrated the performance of our approach employing several quantum variants of CMX through the classical emulations of quantum circuits on the H<sub>2</sub> molecule potential energy surface and the Anderson model with a broad range of correlation strength. A good agreement with exact solutions can be maintained even at the dissociation and strong correlation limits. An essential aspect of this research was integrating a new class of correlated energy functionals with low-depth VQE expansions to provide a source for the trail functions. We demonstrated that the accuracy of this combined approach is superior to the standard VQE methods. Due to the natural resilience of the quantum CMX algorithms to the effect of noise, they are ideally suited to be a target for early quantum simulations on NISQ devices.

**Significance and Impact:** This approach accurately reconstructs a molecular system's total energy using fewer cycles of calculation and reduced numbers of qubits in inherently error-prone quantum circuits.






## Variational Quantum Eigensolver for Approximate Diagonalization of Downfolded Hamiltonians using Generalized Unitary Coupled Cluster Ansatz
<p align="center">
  <img width="400" src="https://github.com/npbauman/BES-QIS/blob/17b2162f6bd2f36d7402d0fad42238896c87a840/Figures/GUCC-VQE_N2.png"> <br>
<sub><sup> Total energies for the dissociation of N<sub>2</sub>.</sup></sub>
</p>

"Variational Quantum Eigensolver for Approximate Diagonalization of Downfolded Hamiltonians using Generalized Unitary Coupled Cluster Ansatz," Nicholas P. Bauman, Jaroslav Chládek, Libor Veis, Jiří Pittner Karol Kowalski, [*arXiv*:2011.01985](https://arxiv.org/abs/2011.01985).
  
**Challenge:** The applicability of quantum algorithms is either limited by the number of variational variables and the numbers of necessary measurements or by the circuit depths. Their accuracy is reliant on underlying approximations. Understanding how these approximations perform, especially when coupled with techniques for reducing the system's dimensionality, is crucial for applications on quantum computers.

**Approach and Results:** We investigated the utilization of variational quantum solver (VQE) and the recently introduced generalized unitary coupled cluster (GUCC) formalism for the diagonalization of downfolded/effective Hamiltonians in active spaces. We also considered various solvers to identify solutions of the GUCC equations using N<sub>2</sub>, H<sub>2</sub>O, and C<sub>2</sub>H<sub>4</sub> as benchmark systems to illustrate the combined framework's performance. Our numerical data suggest that UGCCSD can yield energies and wave functions close to the FCI solution, certainly superior to the UCCSD Ansatz. We also showed that the initial guesses in the form of matrix product states, which can be efficiently prepared on a quantum register, can improve the VQE-GUCC method's accuracy for multireference systems.

**Significance and Impact:** This work promotes the GUCC methodology as a reliable formalism for quantum computing applications. We observed that initial guesses are an influential component of quantum calculations and warrant further investigation. 

https://nwchemgit.github.io/ (NWChem Documentation)  





## Reduced-size plane-wave based representation of many-body methods for quantum computing

"Quantum Solvers for Plane-Wave Hamiltonians: Abridging Virtual Spaces Through the Optimization of Pairwise Correlations,"  Eric Bylaska, Duo Song, Nicholas P. Bauman, Karol Kowalski, Daniel Claudino and Travis S. Humble, [*Front. Chem.*](https://www.frontiersin.org/articles/10.3389/fchem.2021.603019/abstract); [*arXiv*:2009.00080](https://arxiv.org/abs/2009.00080).

**Challenge:** The proper choice of molecular basis set is an essential factor contributing to the efficiency of quantum algorithms in applications to quantum chemistry. This issue is critical in VQE unitary CC and QPE formulations, where a compact representation of virtual orbitals plays a crucial role in the accuracy and efficiency of quantum algorithms. 

**Approach and Results:** Our approach was to define orbital space so that virtual sub-space can capture a significant amount of electron-electron correlation in the system. Due to its efficiency in reaching the complete-basis-set-limit and potential application in quantum computing, the plane-wave basis set is an ideal choice. However, the virtual orbitals in a pseudopotential plane-wave Hartree-Fock calculation are often scattering states that interact very weakly with the filled orbitals because of Coulomb repulsion. As a result, very little correlation energy is captured from them. To overcome these limitations, we have been developing new algorithms to define virtual spaces by optimizing orbitals from small pairwise CI Hamiltonians. We term resulting orbitals as correlation optimized virtual orbitals (COVOs). With these procedures, we have been able to derive virtual spaces containing only a few orbitals that can capture a significant amount of correlation. 

**Significance and Impact:** Using these derived basis sets for quantum computing calculations targeting full CI (FCI) quality-results opens up the door to many-body calculations for pseudopotential plane-wave basis set methods. In conjunction with downfolding methods, the COVOs provide yet another mechanism for dimensionality reduction and maintaining the desired level of accuracy in quantum simulations for chemical processes.

https://nwchemgit.github.io/ (NWChem Documentation)  
https://arrows.emsl.pnnl.gov/api/qsharp_chem (EMSL Arrows Quantum Editor)  
https://nwchemgit.github.io/EMSL_Arrows.html (EMSL Arrows Documentation)  
https://docs.microsoft.com/en-us/azure/quantum/ (Quantum Development Kit Documentation)  





## Benchmarking adaptive variational quantum eigensolvers
<p align="center">
  <img width="400" src="https://github.com/npbauman/BES-QIS/blob/d56ca1b91c28e0555653ee522d547d8b60258094/Figures/VQE_NaH.png"> <br>
<sub><sup> Potential energy curves of NaH computed with the STO-3G basis set.</sup></sub>
</p>

"Benchmarking Adaptive Variational Quantum Eigensolvers,"  Daniel Claudino, Jerimiah Wright, Alexander J. McCaskey and Travis S. Humble, [*Front. Chem.* (2020)](https://doi.org/10.3389/fchem.2020.606863); [*arXiv*:2011.01279](https://arxiv.org/abs/2011.01279).

**Challenge:** The long-term success of quantum computing will depend on the efficiency with which the algorithms solve important problems, and this work addresses the question of what makes a good algorithm for electronic structure calculations from quantum chemistry. While the results are limited to the ground-states of diatomic molecules solved using noiseless simulation, they provide ground truth for how well such variational methods may perform. In addition, this work supports the broader goal of establishing a diverse library of quantum algorithms.

**Approach and Results:** The accuracy of quantum computing methods to calculate the electronic ground states and potential energy curves for selected diatomic molecules, namely H<sub>2</sub>, NaH, and KH has been investigated. Using numerical simulation, it was found that multiple methods provide good estimates of the energy and ground state, but only some methods were shown to be robust to the underlying optimization methods. An important finding from this work is that current, gradient-based optimizations is more economical and delivers superior performance than analogous simulations carried out with gradient-free optimizers. The results also identify small errors in the prepared state fidelity which show an increasing trend with molecular size.

**Significance and Impact:** The long-term success of quantum computing will depend on the efficiency with which the algorithms solve important problems, and this work addresses the question of what makes a good algorithm for electronic structure calculations from quantum chemistry. While the results are limited to the ground-states of diatomic molecules solved using noiseless simulation, they provide ground truth for how well such variational methods may perform. In addition, this work supports the broader goal of establishing a diverse library of quantum algorithms.






## Even more efficient quantum computations of chemistry through tensor hypercontractions

"Even more efficient quantum computations of chemistry through tensor hypercontractions," Joonho Lee, Dominic W. Berry, Craig Gidney, William J. Huggins, Jarrod R. McClean, Nathan Wiebe and Ryan Babbush, [*arXiv*:2011.03494](https://arxiv.org/abs/2011.03494).

**Abstract:** We describe quantum circuits with only **O(N)** Toffoli complexity that block encode the spectra of quantum chemistry Hamiltonians in a basis of **N** arbitrary (e.g., molecular) orbitals. With **O(λ/ϵ)** repetitions of these circuits one can use phase estimation to sample in the molecular eigenbasis, where **λ** is the 1-norm of Hamiltonian coefficients and **ϵ** is the target precision. This is the lowest complexity that has been shown for quantum computations of chemistry within an arbitrary basis. Furthermore, up to logarithmic factors, this matches the scaling of the most efficient prior block encodings that can only work with orthogonal basis functions diagonalizing the Coloumb operator (e.g., the plane wave dual basis). Our key insight is to factorize the Hamiltonian using a method known as tensor hypercontraction (THC) and then to transform the Coulomb operator into an isospectral diagonal form with a non-orthogonal basis defined by the THC factors. We then use qubitization to simulate the non-orthogonal THC Hamiltonian, in a fashion that avoids most complications of the non-orthogonal basis. We also reanalyze and reduce the cost of several of the best prior algorithms for these simulations in order to facilitate a clear comparison to the present work. In addition to having lower asymptotic scaling spacetime volume, compilation of our algorithm for challenging finite-sized molecules such as FeMoCo reveals that our method requires the least fault-tolerant resources of any known approach. By laying out and optimizing the surface code resources required of our approach we show that FeMoCo can be simulated using about four million physical qubits and under four days of runtime, assuming 1 μs cycle times and physical gate error rates no worse than 0.1%. 

<!---**Challenge:**--->

<!---**Approach and Results:**--->

<!---**Significance and Impact:**--->






## Dynamical self-energy mapping for quantum computing

"Dynamical Self-energy Mapping (DSEM) for quantum computing," Diksha Dhawan, Mekena Metcalf and Dominika Zgid, [*arXiv*:2010.05441](https://arxiv.org/abs/2010.05441).

**Abstract:** For noisy intermediate-scale quantum (NISQ) devices only a moderate number of qubits with a limited coherence is available thus enabling only shallow circuits and a few time evolution steps in the currently performed quantum computations. Here, we present how to bypass this challenge in practical molecular chemistry simulations on NISQ devices by employing a classical-quantum hybrid algorithm allowing us to produce a sparse Hamiltonian which contains only **O(n<sup>2</sup>)** terms in a Gaussian orbital basis when compared to the **O(n<sup>4</sup>)** terms of a standard Hamiltonian, where **n** is the number of orbitals in the system. Classical part of this hybrid entails parameterization of the sparse, fictitious Hamiltonian in such a way that it recovers the self-energy of the original molecular system. Quantum machine then uses this fictitious Hamiltonian to calculate the self-energy of the system. We show that the developed hybrid algorithm yields very good total energies for small molecular test cases while reducing the depth of the quantum circuit by at least an order of magnitude when compared with simulations involving a full Hamiltonian. 

<!---**Challenge:**--->

<!---**Approach and Results:**--->

<!---**Significance and Impact:**--->





## Simulating quantum materials with digital quantum computers

"Simulating Quantum Materials with Digital Quantum Computers," Lindsay Bassman, Miroslav Urbanek, Mekena Metcalf, Jonathan Carter, Alexander F. Kemper and Wibe de Jong, [*arXiv*:2101.08836](https://arxiv.org/abs/2101.08836).

**Abstract:** Quantum materials exhibit a wide array of exotic phenomena and practically useful properties. A better understanding of these materials can provide deeper insights into fundamental physics in the quantum realm as well as advance technology for entertainment, healthcare, and sustainability. The emergence of digital quantum computers (DQCs), which can efficiently perform quantum simulations that are otherwise intractable on classical computers, provides a promising path forward for testing and analyzing the remarkable, and often counter-intuitive, behavior of quantum materials. Equipped with these new tools, scientists from diverse domains are racing towards achieving physical quantum advantage (i.e., using a quantum computer to learn new physics with a computation that cannot feasibly be run on any classical computer). The aim of this review, therefore, is to provide a summary of progress made towards this goal that is accessible to scientists across the physical sciences. We will first review the available technology and algorithms, and detail the myriad ways to represent materials on quantum computers. Next, we will showcase the simulations that have been successfully performed on currently available DQCs, emphasizing the variety of properties, both static and dynamic, that can be studied with this nascent technology. Finally, we work through two examples of how to map a materials problem onto a DQC, with full code included in the Supplementary Material. It is our hope that this review can serve as an organized overview of progress in the field for domain experts and an accessible introduction to scientists in related fields interested in beginning to perform their own simulations of quantum materials on DQCs. 

<!---**Challenge:**--->

<!---**Approach and Results:**--->

<!---**Significance and Impact:**--->





## Variational quantum solver employing the PDS energy functional

"Variational quantum solver employing the PDS energy functional," Bo Peng and Karol Kowalski, [*arXiv*:2101.08526](https://arxiv.org/abs/2101.08526).
  
**Challenge:** To capture a subtle balance between static and dynamical correlations effect in quantum computing, the conventional variational quantum eigensolver (VQE) approach requires incorporating more and more parameters in the anasatz state preparation, which on the other hand makes the VQE approach vulnerable to complex quantum circuit and staggered optimization on classical machines.

**Approach and Results:** Here we find that a new class of objective cost function based on the Peeters-Devreese-Soldatov (PDS) formulation can be employed in a conventional VQE framework. In comparison with the conventional VQE and the static PDS approach, this new variational quantum solver offers an effective approach that helps navigate the optimization dynamics on the upper bound energy landscape of the target state, thus let the dynamics be free from getting trapped in the local minima that refer to different states. We demonstrate the performance of the proposed variational quantum solver for toy models, H<sub>2</sub> molecule, and strongly correlated planar H<sub>4</sub> system in some challenging situations. In all the case studies, the proposed variational quantum approach outperforms the usual VQE and static PDS calculations even at the lowest order.

**Significance and Impact:** The proposed variational quantum solver employing the PDS functionals helps achieve high accuracy at finding the ground state and its energy through the rotation of the trial wave function of modest quality, thus improves the accuracy and efficiency of the quantum simulation. 





## Dimensionality reduction of many-body problem using coupled-cluster sub-system flow equations
"Dimensionality reduction of many-body problem using coupled-cluster sub-system flow equations: classical and quantum computing perspective," Karol Kowalski, [*arXiv*:2102.05783](https://arxiv.org/abs/2102.05783).
  
**Challenge:** Efficient probing of the large sub-spaces by quantum algorithms emerges as a critical challenge for accurate quantum simulations of chemical systems.

**Approach and Results:** To address this challenge, we have extended the active-space downfolding formalisms into quantum flow algorithm  (QFA), where the high-dimensionality quantum problem is partitioned into coupled low-dimensionality eigenvalue problems amenable to quantum computing. This formalism precisely defines the communication protocol  between sub-problems involved in the flow and can be easily extended to the time domain. 

**Significance and Impact:** The QFA formalism enables to probe large sub-spaces of Hilbert space in a size-consistent manner without significant increase in quantum resources. It is also ideally suited to develop quantum algorithms for local formulations of quantum chemical methods. 





## Improving the accuracy and efficiency of quantum connected moments expansions
<p align="center">
  <img width="400" src="https://github.com/npbauman/BES-QIS/blob/a347c6fa28dfb40cfc0d7c91cf938a1f03a4d0dd/Figures/CMX-H2.png"> <br>
<sub><sup> Potential energy curves of H<sub>2</sub> in the 6-31G basis.</sup></sub>
</p>

"Improving the accuracy and efficiency of quantum connected moments expansions," Daniel Claudino, Bo Peng, Nicholas P. Bauman, Karol Kowalski and Travis S. Humble, [*arXiv*:2103.09124](https://arxiv.org/abs/2103.09124).

**Abstract:** The still-maturing noisy intermediate-scale quantum (NISQ) technology faces strict limitations on the algorithms that can be implemented efficiently. In quantum chemistry, the variational quantum eigensolver (VQE) algorithm has become ubiquitous, using the functional form of the ansatz as a degree of freedom, whose parameters are found variationally in a feedback loop between the quantum processor and its conventional counterpart. Alternatively, a promising new avenue has been unraveled by the quantum variants of techniques grounded on expansions of the moments of the Hamiltonian, among which two stand out: the connected moments expansion (CMX) [Phys. Rev. Lett. 58, 53 (1987)] and the Peeters-Devreese-Soldatov (PDS) functional [J. Phys. A 17, 625 (1984); Int. J. Mod. Phys. B 9, 2899], the latter based on the standard moments **<H<sup>k</sup>>**. Contrasting with VQE-based methods and provided the quantum circuit prepares a state with non-vanishing overlap with the true ground state, CMX often converges to the ground state energy, while PDS is guaranteed to converge by virtue of being variational. However, for a finite CMX/PDS order, the circuit may significantly impact the energy accuracy. Here we use the ADAPT-VQE algorithm to test shallow circuit construction strategies that are not expected to impede their implementation in the present quantum hardware while granting sizable accuracy improvement in the computed ground state energies. We also show that we can take advantage of the fact that the terms in the connected moments are highly recurring in different powers, incurring a sizable reduction in the number of necessary measurements. By coupling this measurement caching with a threshold that determines whether a given term is to be measured based on its associated scalar coefficient, we observe a further reduction in the number of circuit implementations while allowing for tunable accuracy. 

**Challenge:** The accuracy of the ground-state energy based on expansion of the moments of the Hamiltonian for a given expansion order is dependent on the degree of overlap between the prepared state and the true ground state. Increasing this overlap needs to be done in order to strike a balance between accuracy and feasibility of the circuits to be implemented. The accuracy in the energy estimates increases with the expansion order, at the expense of a potential exponential increase in the number of terms to be measured. 

**Approach and Results:** We take advantage of the intrinsic iterative circuit construction from the ADAPT-VQE algorithm to prepare a state upon a single ADAPT iteration, i.e., one many-body rotation away from the Hartree-Fock state. We find that this is enough to provide substantial improvement in the energy estimates, as it is shown that out approach is enough to retain a large overlap with the ground state, particularly in the regime where HF is not a good approximation to the ground state. For a given Hamiltonian, there is a limited number of unique terms that are found in the low orders that are recurrent for all moments, so caching the measured terms drastically decreases the number of necessary measurements. This can be further alleviated by setting a numerical threshold in the scalar multiplying the terms to be measured to signal terms whose contribution may be deemed numerically unimportant.

**Significance and Impact:** Our findings strengthen the case for quantum moments expansions as an efficient and accurate alternative to strategies based on the variational quantum eigensolver.

https://nwchemgit.github.io/ (NWChem Documentation)
https://xacc.readthedocs.io/en/latest/ (XACC Documentation)






## Numerical Simulations of Noisy Quantum Circuits for Computational Chemistry
<p align="center">
  <img width="400" src="https://github.com/danclaudino/QIS-Highlights/blob/master/Figures/noisy-vqe.png"> <br>
<sub><sup> Simulation results of the energy (top row) and state fidelity (bottom row) for NaH using variational quantum circuits in the presence of varying noise levels of hardware noise. Left column is the bare circuit while the right column uses randomized compiling.</sup></sub>
</p>

"Numerical Simulations of Noisy Quantum Circuits for Computational Chemistry," Meenambika Gowrishankar, Jerimiah Wright, Daniel Claudino, Thien Nguyen, Alexander McCaskey, and Travis S. Humble, [2021 IEEE International Conference on Quantum Computing and Engineering (QCE)](https://ieeexplore.ieee.org/document/9605344).

**Abstract:** This is a case study of the variational quantum eigensolver (VQE) method using numerical simulations to test the influence of noise on the accuracy of the underlying circuit ansatz. We investigate a computational chemistry application of VQE to calculate the electronic ground state and its energy for Sodium Hydride (NaH), a prototypical two-electron problem. Using a one-parameter ansatz derived from unitary coupled cluster (UCC) theory, we simulate the effects of noise on the energy expectation value and variance with respect to the ansatz parameter. These numerical simulations provide insights into the accuracy of the prepared quantum state and the efficiency of the classical optimizer that iteratively refines the ansatz. We conduct a comparative study between analytical results derived for the UCC ansatz in the absence of noise and the noisy numerical simulation results obtained using an isotropic depolarizing noise model for each gate. We also compare the relative increase in noise on logically equivalent UCC ansatz circuits generated by randomized compiling. Notably, we observe that the intrinsic variance in the energy due to the simplicity of the ansatz itself compares with the noise induced by the bare circuit. 

**Challenge:** Performance bounds for variational algorithms are known largely in the absence of noise. As the quantum processors evolve, it is imperative to develop a better understanding of the effects they introduce for such algorithms.

**Approach and Results:** We analyze the performance of noisy quantum computers to calculate the electronic energy of small molecules using numerical simulations. Our results quantify the influence of noise on energy accuracy and state fidelity.

**Significance and Impact:** Estimating the intrinsic accuracy of quantum computational methods is important for setting performance expectations of algorithms and performance requirements for hardware.




## Coupled cluster downfolding methods: The effect of double commutator terms on the accuracy of ground-state energies

"Coupled cluster downfolding methods: The effect of double commutator terms on the accuracy of ground-state energies," Nicholas P. Bauman and Karol Kowalski, [*J. Chem. Phys.* (2022)](https://doi.org/10.1063/5.0076260); [*arXiv*:2011.12077](https://doi.org/10.48550/arXiv.2110.12077).

**Abstract:** Downfolding coupled cluster techniques have recently been introduced into quantum chemistry as a tool for the dimensionality reduction of the many-body quantum problem. As opposed to earlier formulations in physics and chemistry based on the concept of effective Hamiltonians, the appearance of the downfolded Hamiltonians is a natural consequence of the single-reference exponential parameterization of the wave function. In this paper, we discuss the impact of higher-order terms originating in double commutators. In analogy to previous studies, we consider the case when only one- and two-body interactions are included in the downfolded Hamiltonians. We demonstrate the efficiency of the many-body expansions involving single and double commutators for the unitary extension of the downfolded Hamiltonians on the example of the beryllium atom, and bond-breaking processes in the Li<sub>2</sub> and H<sub>2</sub>O molecules. For the H<sub>2</sub>O system, we also analyze energies obtained with downfolding procedures as functions of the active space size.


**Significance and Impact:** For the H<sub>2</sub>O molecule, we demonstrated that expansion based on the inclusion of single, double, and selected triple commutators reachs the CCSDTQ level of accuracy, when all orbitals are correlated. A newly developed symbolic manipulation equation generator for the double unitary coupled cluster formalism (SyMan-DUCC) has been integrated with the exa-scale Tensor Algebra for Many-body Methods library (TAMM), which paves the way for downfolding systems composed of 1,000+ orbitals to arbitrary size active space. 

