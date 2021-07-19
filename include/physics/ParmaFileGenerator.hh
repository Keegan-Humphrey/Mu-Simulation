#ifndef MU__PHYSICS__PARMA_FILE_GENERATOR_HH
#define MU__PHYSICS__PARMA_FILE_GENERATOR_HH

#include "physics/Generator.hh"

#include "physics/Particle.hh"
#include "ui.hh"

#include <string>
#include <iostream>
#include <cstddef>
#include <vector>

namespace MATHUSLA { namespace MU { namespace Physics {

class ParmaFileGenerator : public Generator {
public:
  ParmaFileGenerator(const std::string &name, const std::string &description);

  virtual ~ParmaFileGenerator() = default;

  virtual void GeneratePrimaryVertex(G4Event *event);
  virtual void SetNewValue(G4UIcommand *command, G4String value);
  virtual std::ostream &Print(std::ostream &os = std::cout) const;
  virtual const Analysis::SimSettingList GetSpecification() const;

protected:
  virtual void GenerateCommands();

  std::size_t _event_counter = 0;
  std::vector<Particle> _particle_parameters;

  Command::StringArg *_ui_pathname;
};

} } } // namespace MATHUSLA::MU::Physics

#endif // MU__PHYSICS__PARMA_FILE_GENERATOR_HH