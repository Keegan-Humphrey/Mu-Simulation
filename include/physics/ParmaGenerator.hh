#ifndef MU__PHYSICS__PARMA_GENERATOR_HH
#define MU__PHYSICS__PARMA_GENERATOR_HH

#include "physics/Generator.hh"

#include "physics/Particle.hh"
#include "physics/Parma.hh"
#include "ui.hh"

#include <string>
#include <iostream>
#include <cstddef>
#include <vector>

namespace MATHUSLA { namespace MU { namespace Physics {

class ParmaGenerator : public Generator {
public:
  ParmaGenerator(Parma4::Parma* parma);
  ParmaGenerator(std::vector<std::string>& settings);

  ~ParmaGenerator() = default;

  void GeneratePrimaryVertex(G4Event *event);
  virtual GenParticleVector GetLastEvent() const;
  // void SetNewValue(G4UIcommand *command, G4String value);
  void SetParma(Parma4::Parma* parma);
  void SetParma(std::vector<std::string>& settings);

  virtual const Analysis::SimSettingList GetSpecification() const;

protected:
  static G4ThreadLocal Parma4::Parma* _parma;
  static G4ThreadLocal std::vector<std::string>* _parma_settings;
  static G4ThreadLocal bool _settings_on;

  GenParticleVector _last_event;
  std::uint_fast64_t _counter;
  
  double _x0_min, _x0_max, _y0_min, _y0_max, _z0;

  Command::DoubleUnitArg*     _ui_x0_min;
  Command::DoubleUnitArg*     _ui_x0_max;
  Command::DoubleUnitArg*     _ui_y0_min;
  Command::DoubleUnitArg*     _ui_y0_max;
  Command::DoubleUnitArg*     _ui_z0;
};

} } } // namespace MATHUSLA::MU::Physics

#endif // MU__PHYSICS__PARMA_GENERATOR_HH