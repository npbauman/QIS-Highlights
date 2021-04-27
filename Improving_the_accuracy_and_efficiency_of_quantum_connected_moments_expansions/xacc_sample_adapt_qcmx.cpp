#include "xacc.hpp"
#include "xacc_observable.hpp"
#include "xacc_service.hpp"
#include <iomanip>

int main(int argc, char **argv) {
  xacc::Initialize(argc, argv);
 
  // read hamiltonian from file
  std::ifstream hfile("/path/to/hamiltonian/file");
  std::stringstream ss;
  ss << hfile.rdbuf();
  // representation is string, either pauli or fermion
  auto H = xacc::quantum::getObservable(representation, ss.str());

  // accelerator
  auto acc = xacc::getAccelerator("qpp");

  // optimizer
  auto opt = xacc::getOptimizer("nlopt", {{"algorithm", "cobyla"}, {"ftol-abs", 1e-7}, {"xtol", 1e-5}});

  // set up adapt
  // nElectrons is the number of electrons
  // pool is the pauli pool
  auto adapt = xacc::getAlgorithm("adapt", {{"accelerator", accelerator}, 
					{"observable", H},
					{"n-electrons", nElectrons},
					{"sub-algorithm", "vqe"},
					{"pool", "qubit-pool"},
					{"optimizer", optmizer}});
  // allocate buffer and execute
  auto q = xacc::qalloc(H->nBits());
  adapt->execute(q);

  // get circuit string and energy
  auto xasmStr = q->getInformation("opt-circuit").as<std::string>();
  auto ansatz = xasm->compile(xasmStr)->getComposites()[0];
  auto adaptEnergy = q->getInformation("opt-val").as<double>();

  // cmx order
  auto cmx_order = 2;

  // get reference to qcmx algorithm and initialize
  auto qcmx = xacc::getService<xacc::Algorithm>("qcmx");
  qcmx->initialize({{"accelerator", accelerator},
		    {"observable", H},
		    {"ansatz", ansatz},
                    {"cmx-order", cmx_order}});
  
  // allocate buffer, execute
  auto buffer = xacc::qalloc(nOrbitals);
  qcmx->execute(buffer);

  // retrieve map with energies for all expansions up to cmx_order
  auto map = buffer->getInformation("energies").as<std::map<std::string, double>>();

  for (auto & item : map) {
    std::cout << std::setprecision(12) << item.first << " = " << item.second << "\n";
  } 

  xacc::Finalize();
  return 0;
}


  xacc::Finalize();
  return 0;
}

