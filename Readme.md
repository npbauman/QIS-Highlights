This repository is dedicated to highlighting developments within the Quantum Information Science (QIS) Initiative at the Pacific Northwest National Laboratory (PNNL). It also provides some links, software, and data pertinent to the projects as well. 

**Table of Contents**
- [Development of the NWChem-QDK interface](#Development-of-the-NWChem-QDK-interface)
- [Downfolding techniques for dimension reduction of correlated electronic Hamiltonians](#Downfolding-techniques-for-dimension-reduction-of-correlated-electronic-Hamiltonians)
- [Resource-efficient VQE algorithms with downfolded Hamiltonians](#Resource-efficient-VQE-algorithms-with-downfolded-Hamiltonians)
- [Quantum computing for high-energy states - QPE simulations of core-level states](#Quantum-computing-for-high-energy-states---QPE-simulations-of-core-level-states)
- [Quantum algorithms for connected moments expansion](#Quantum-algorithms-for-connected-moments-expansion)
- [Reduced-size plane-wave based representation of many-body methods for quantum computing](#Reduced-size-plane-wave-based-representation-of-many-body-methods-for-quantum-computing)
- [Benchmarking adaptive variational quantum eigensolvers](#Benchmarking-adaptive-variational-quantum-eigensolvers)






## Development of the NWChem-QDK interface
<p align="center">
  <img width="400" src="https://github.com/npbauman/BES-QIS-at-PNNL/blob/a14900427b39cae1ccb64c3d2272b8d2782e5465/Figures/QuantumWorkflow.png"> <br>
<sub><sup>Detailed workflow for simulating quantum chemistry on quantum computers.</sup></sub>
</p>
  
"Q# and NWChem: Tools for Scalable Quantum Chemistry on Quantum Computers," Guang Hao Low, Nicholas P. Bauman, Christopher E. Granade, Bo Peng, Nathan Wiebe, Eric J. Bylaska, Dave Wecker, Sriram Krishnamoorthy, Martin Roetteler, Karol Kowalski, Matthias Troyer, Nathan A. Baker, [*arXiv*:1904.01131](https://arxiv.org/abs/1904.01131)    

**Challenge:** Enabling new quantum computing simulations in chemistry requires tight integration of computational chemistry infrastructure with novel quantum algorithms for ground- and excited-state simulations.

**Approach and Results:** Our team has developed NWChem-QDK (Quantum Development Kit (QDK) developed and maintained by Microsoft Research group) interface to perform Quantum Phase Estimator (QPE) simulations for chemical processes. NWChem provides a reach infrastructure for the characterization of second quantized forms of electronic Hamiltonians and the initial characterization of ground- and excited-state wavefunction obtained with high-accuracy coupled-cluster calculations on classical systems. QDK offers a variety of quantum simulation algorithms ranging from Trotter-Suzuki expansion to various qubitization approaches. Quantum simulations using NWChem and QDK can also be performed employing the web interface EMSL Arrows Quantum Editor.

**Significance and Impact:** Using the NWChem-QDK interface, we demonstrated the efficiency of the QPE approach in describing strongly correlated ground and excited states of molecular systems. The QPE algorithm, with proper initial estimates of the electronic wave functions, was able to describe potential energy surfaces for the ground state and several low-lying excited states, which also involved challenging doubly excited states. Using the H10 benchmark system in the STO-3G basis set, we demonstrated the advantages of using QPE in achieving highly accurate energy estimates in the strongly correlated regime.

https://nwchemgit.github.io/ (NWChem Documentation)  
https://arrows.emsl.pnnl.gov/api/qsharp_chem (EMSL Arrows Quantum Editor)  
https://nwchemgit.github.io/EMSL_Arrows.html (EMSL Arrows Documentation)  
https://docs.microsoft.com/en-us/azure/quantum/ (Quantum Development Kit Documentation)  



## Downfolding techniques for dimension reduction of correlated electronic Hamiltonians
<p align="center">
  <img width="400" src="https://github.com/npbauman/BES-QIS-at-PNNL/blob/d71bfd80cfc7c194d80bb1fabb1315a6a853780d/Figures/DownfoldingAbstract.jpg"> <br>
<sub><sup>A schematic representation of the coupled cluster (CC) downfolding technique, where the CC Ansatz is used to generate effective/downfolded Hamiltonians in small-dimensionality active spaces.</sup></sub>
</p>
  
"Downfolding of many-body Hamiltonians using active-space models: Extension of the sub-system embedding sub-algebras approach to unitary coupled cluster formalisms," Nicholas P. Bauman, Eric J. Bylaska, Sriram Krishnamoorthy, Guang Hao Low, Nathan Wiebe, Christopher E. Granade, Martin Roetteler, Matthias Troyer, and Karol Kowalski, [*J. Chem. Phys.* **151**, 014107 (2019)](https://doi.org/10.1063/1.5094643); [*arXiv*:1902.01553](https://arxiv.org/abs/1902.01553)    

**Challenge:** Limited quantum resources preclude simulations of complex and realistic chemical processes. Therefore, quantum computing is in high demand for efficient techniques for re-representing quantum many-body problems in reduced dimensionality spaces before reaching maturity.

**Approach and Results:** We have extended the sub-system embedding sub-algebras coupled cluster (SES-CC) theory to the downfolding procedure based on the double unitary CC (DUCC) formalism to address this challenge. In contrast to the standard single-reference SES-CC formulations, the DUCC approach results in a Hermitian form of the effective Hamiltonian in an active-space, which provides a rigorous separation of external cluster amplitudes that describe dynamical correlation effects from those corresponding to the internal (within the active space) excitations that define the components of eigenvectors associated with the energy of the entire system. 

**Significance and Impact:** We have extended the sub-system embedding sub-algebras coupled cluster (SES-CC) theory to the downfolding procedure based on the double unitary CC (DUCC) formalism to address this challenge. In contrast to the standard single-reference SES-CC formulations, the DUCC approach results in a Hermitian form of the effective Hamiltonian in an active-space, which provides a rigorous separation of external cluster amplitudes that describe dynamical correlation effects from those corresponding to the internal (within the active space) excitations that define the components of eigenvectors associated with the energy of the entire system.

https://nwchemgit.github.io/ (NWChem Documentation)  
https://docs.microsoft.com/en-us/azure/quantum/ (Quantum Development Kit Documentation)



## Resource-efficient VQE algorithms with downfolded Hamiltonians
**Challenge:** In modeling many-body problems, the biggest challenge confronted is that the number of qubits scales linearly with the molecular basis's size. This poses a significant limitation on the basis sets' size and the number of correlated electrons included in quantum simulations of chemical processes.

**Approach and Results:** To address this issue and enable more realistic simulations on NISQ computers, we employed the double unitary coupled-cluster method to effectively downfold correlation effects into the reduced-size orbital space, commonly referred to as the active space. Using downfolding and VQE techniques, we demonstrated that effective Hamiltonians could capture the effect of the whole orbital space in small-size active spaces, especially when natural orbitals are employed. 

**Significance and Impact:** The DUCC/VQE framework was used to solve the ground-state energy of H2, Li2, and BeH2 on the cc-pVTZ basis using the reduced-size active spaces. The VQE formalism has also been extended to the generalized unitary CC (GUCC) ansatz and applied to benchmark systems/processes (N2, H2O, and C2H4) described by the downfolded Hamiltonians. The preliminary data indicate that simple downfolding procedures based on single commutator expansion and small active spaces can recover around 90% of correlation energy calculated when all orbitals are correlated.



## Quantum computing for high-energy states - QPE simulations of core-level states
**Challenge:** The experimental effort at the BES-supported light source facilities is contingent on the availability of predictive modeling tools to interpret x-ray experimental spectra. However, the theoretical description and numerical identification of complicated core-level states still pose a significant challenge for approximate methods.

**Approach and Results:** The stochastic nature of the QPE algorithm can facilitate the discovery/identification of core-level states when the knowledge about the true configuration structure of a sought-after excited state is limited or postulated. For this purpose, we developed an algorithm where through repeated simulations, one can accumulate samples from the distribution of eigenstates energies. The desired error in each energy estimate is inversely proportional to the number of applications of the time evolution operator in the QPE algorithm. This is in contrast to VQE approaches, which only provide energy estimates for a single targeted state. 

**Significance and Impact:** The PNNL-Microsoft team has demonstrated that extension of the QPE algorithm to high-energy core-level states of various spin, spatial symmetries, and complexity is possible. This is the first demonstration of quantum algorithms' potential for identifying shake-up/satellite states and their future role in supporting various x-ray spectroscopies.



## Quantum algorithms for connected moments expansion
**Challenge:** Further advancement of quantum computing is contingent on enabling models that avoid deep circuits (Fig.7) and the excessive use of CNOT gates and associated noise effects.

**Approach and Results:** We developed a quantum algorithm employing finite-order connected moment expansions (CMX) and affordable procedures for initial-state preparation to address this problem. We demonstrated the performance of our approach employing several quantum variants of CMX through the classical emulations of quantum circuits on the H2 molecule potential energy surface and the Anderson model with a broad range of correlation strength. A good agreement with exact solutions can be maintained even at the dissociation and strong correlation limits. An essential aspect of this research was integrating a new class of correlated energy functionals with low-depth VQE expansions to provide a source for the trail functions. We demonstrated that the accuracy of this combined approach is superior to the standard VQE methods. Due to the natural resilience of the quantum CMX algorithms to the effect of noise, they are ideally suited to be a target for early quantum simulations on NISQ devices.

**Significance and Impact:** This approach accurately reconstructs a molecular system's total energy using fewer cycles of calculation and reduced numbers of qubits in inherently error-prone quantum circuits.



## Reduced-size plane-wave based representation of many-body methods for quantum computing
**Challenge:** The proper choice of molecular basis set is an essential factor contributing to the efficiency of quantum algorithms in applications to quantum chemistry. This issue is critical in VQE unitary CC and QPE formulations, where a compact representation of virtual orbitals plays a crucial role in the accuracy and efficiency of quantum algorithms. 

**Approach and Results:** Our approach was to define orbital space so that virtual sub-space can capture a significant amount of electron-electron correlation in the system. Due to its efficiency in reaching the complete-basis-set-limit and potential application in quantum computing, the plane-wave basis set is an ideal choice. However, the virtual orbitals in a pseudopotential plane-wave Hartree-Fock calculation are often scattering states that interact very weakly with the filled orbitals because of Coulomb repulsion. As a result, very little correlation energy is captured from them. To overcome these limitations, we have been developing new algorithms to define virtual spaces by optimizing orbitals from small pairwise CI Hamiltonians. We term resulting orbitals as correlation optimized virtual orbitals (COVOs). With these procedures, we have been able to derive virtual spaces containing only a few orbitals that can capture a significant amount of correlation. 

**Significance and Impact:** Using these derived basis sets for quantum computing calculations targeting full CI (FCI) quality-results opens up the door to many-body calculations for pseudopotential plane-wave basis set methods. In conjunction with downfolding methods, the COVOs provide yet another mechanism for dimensionality reduction and maintaining the desired level of accuracy in quantum simulations for chemical processes.



## Benchmarking adaptive variational quantum eigensolvers
**Challenge:** The long-term success of quantum computing will depend on the efficiency with which the algorithms solve important problems, and this work addresses the question of what makes a good algorithm for electronic structure calculations from quantum chemistry. While the results are limited to the ground-states of diatomic molecules solved using noiseless simulation, they provide ground truth for how well such variational methods may perform. In addition, this work supports the broader goal of establishing a diverse library of quantum algorithms.

**Approach and Results:** The accuracy of quantum computing methods to calculate the electronic ground states and potential energy curves for selected diatomic molecules, namely H2, NaH, and KH has been investigated (Fig.8). Using numerical simulation, it was found that multiple methods provide good estimates of the energy and ground state, but only some methods were shown to be robust to the underlying optimization methods. An important finding from this work is that current, gradient-based optimizations is more economical and delivers superior performance than analogous simulations carried out with gradient-free optimizers. The results also identify small errors in the prepared state fidelity which show an increasing trend with molecular size.

**Significance and Impact:** The long-term success of quantum computing will depend on the efficiency with which the algorithms solve important problems, and this work addresses the question of what makes a good algorithm for electronic structure calculations from quantum chemistry. While the results are limited to the ground-states of diatomic molecules solved using noiseless simulation, they provide ground truth for how well such variational methods may perform. In addition, this work supports the broader goal of establishing a diverse library of quantum algorithms.

