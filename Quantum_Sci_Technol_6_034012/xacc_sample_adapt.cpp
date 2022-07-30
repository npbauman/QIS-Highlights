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
  auto initialEnergy = q->getInformation("opt-val").as<double>();

  // print energy, circuit depth, and circuit
  std::cout << std::setprecision(12) << "Energy = " << initialEnergy << "\n"; 
  std::cout << "Circuit depth = " << ansatz->depth() << "\n";
  std::cout << "Circuit\n" << ansatz->toString() << "\n";

  xacc::Finalize();
  return 0;
}

