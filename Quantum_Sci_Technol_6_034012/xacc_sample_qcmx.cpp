#include "xacc.hpp"
#include "xacc_observable.hpp"
#include "xacc_service.hpp"
#include "ObservableTransform.hpp"
#include <iomanip>

int main(int argc, char **argv) {
  xacc::Initialize(argc, argv);

  // read hamiltonian from file
  std::ifstream hfile("/path/to/hamiltonian/file");
  std::stringstream ss;
  ss << hfile.rdbuf();
  // representation is string, either pauli or fermion
  auto H = xacc::quantum::getObservable(representation, ss.str());

  // create circuit to prepare state (HF)
  auto provider = xacc::getService<xacc::IRProvider>("quantum");
  auto ansatz = provider->createComposite("initial-state");
  auto nOrbitals = H->nBits();

  for (std::size_t i = 0; i < ne / 2; i++) {
    ansatz->addInstruction(provider->createInstruction("X", {i}));
    ansatz->addInstruction(provider->createInstruction("X", {i + nOrbitals / 2}));
  }

  // accelerator
  auto accelerator = xacc::getAccelerator("qpp");

  // cmx order
  auto cmx_order = 2;

  // set threshold (optional)
  auto threshold = 0.0;

  // get reference to qcmx algorithm and initialize
  auto qcmx = xacc::getService<xacc::Algorithm>("qcmx");
  qcmx->initialize({{"accelerator", accelerator},
		    {"observable", H},
		    {"ansatz", ansatz},
                    {"cmx-order", cmx_order},
		    {"threshold", threshold}});
  
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

