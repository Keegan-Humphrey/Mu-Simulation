/*
 * src/physics/ParmaGenerator.cc
 *
 * Copyright 2018 Brandon Gomes
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <sstream>

#include "physics/ParmaGenerator.hh"
#include "physics/Parma.hh"

#include "geometry/Earth.hh"
#include "geometry/Cavern.hh"
#include "physics/Units.hh"
#include "util/string.hh"

namespace MATHUSLA { namespace MU {

namespace Physics { ////////////////////////////////////////////////////////////////////////////

//__Parma Generator Constructor_________________________________________________________________
ParmaGenerator::ParmaGenerator(Parma4::Parma* parma) 
  : Generator("parma", "Parma4 Cosmic Spectra Generator.") {
  _parma_settings = new std::vector<std::string>();
  SetParma(parma);
  GenerateCommands();
}
//----------------------------------------------------------------------------------------------

//__Parma Generator Constructor_________________________________________________________________
ParmaGenerator::ParmaGenerator(std::vector<std::string>& settings)
  : ParmaGenerator() {
  SetParma(settings);
  GenerateCommands();
}
//----------------------------------------------------------------------------------------------

ParmaGenerator::ParmaGenerator() : ParmaGenerator() {
  SetParma();
  GenerateCommands();
}

void ParmaGenerator::GenerateCommands() {
  _ui_x0_min = CreateCommand<Command::DoubleUnitArg>("x0_min", "Set minimum x0.");
  _ui_x0_min->SetParameterName("x0_min", false, false);
  _ui_x0_min->SetDefaultUnit("mm");
  _ui_x0_min->SetUnitCandidates("mm cm m");
  _ui_x0_min->AvailableForStates(G4State_PreInit, G4State_Idle);

  _ui_x0_max = CreateCommand<Command::DoubleUnitArg>("x0_max", "Set maximum x0.");
  _ui_x0_max->SetParameterName("x0_max", false, false);
  _ui_x0_max->SetDefaultUnit("mm");
  _ui_x0_max->SetUnitCandidates("mm cm m");
  _ui_x0_max->AvailableForStates(G4State_PreInit, G4State_Idle);

  _ui_y0_min = CreateCommand<Command::DoubleUnitArg>("y0_min", "Set minimum y0.");
  _ui_y0_min->SetParameterName("y0_min", false, false);
  _ui_y0_min->SetDefaultUnit("mm");
  _ui_y0_min->SetUnitCandidates("mm cm m");
  _ui_y0_min->AvailableForStates(G4State_PreInit, G4State_Idle);

  _ui_y0_max = CreateCommand<Command::DoubleUnitArg>("y0_max", "Set maximum y0.");
  _ui_y0_max->SetParameterName("y0_max", false, false);
  _ui_y0_max->SetDefaultUnit("mm");
  _ui_y0_max->SetUnitCandidates("mm cm m");
  _ui_y0_max->AvailableForStates(G4State_PreInit, G4State_Idle);

  _ui_z0 = CreateCommand<Command::DoubleUnitArg>("z0", "Set z0.");
  _ui_z0->SetParameterName("z0", false, false);
  _ui_z0->SetDefaultUnit("mm");
  _ui_z0->SetUnitCandidates("mm cm m");
  _ui_z0->AvailableForStates(G4State_PreInit, G4State_Idle);
}

namespace { ////////////////////////////////////////////////////////////////////////////////////

// // Setup Parma Randomness (TODO)
// Parma4::Parma* _setup_random(Parma4::Parma* parma) {
//   int threadno=G4Threading::G4GetThreadId();
//   pythia->readString("Random:setSeed = on");
//   std::ostringstream oss;
//   oss << "Random:seed = " << threadno;
//   pythia->readString(oss.str());
//   pythia->readString("Next:showScaleAndVertex = on");
//   return pythia;
// }

//__Reconstruct Parma Object from Old Object___________________________________________________
Parma4::Parma* _reconstruct_parma(Parma4::Parma* parma) {
  if (!parma) {
    return new Parma4::Parma();
  } else {
    auto out = new Parma4::Parma(parma.getIP(), parma.getSeed());
    return out;
  }
}
//----------------------------------------------------------------------------------------------
//Create Default Parma
Parma4::Parma* _create_parma() {
  auto parma = new Parma4::Parma();
  parma->init();
  return parma;
}
//__Create Parma from Settings_________________________________________________________________
Parma4::Parma* _create_parma(std::vector<std::string>* settings,
                                bool& settings_on) {
  auto parma = new Parma4::Parma();
  // for (const auto& setting : *settings)
  //   parma->readString(setting);
  // _setup_random(parma);
  parma->init();
  // settings_on = true;
  return parma;
}
//----------------------------------------------------------------------------------------------

//__Convert Parma output array to Geant particle_________________________________________________________
Particle _convert_particle(const double[3] particle) {  //comes in as {energy, zenith, azimuthal}.
  int pdgID;  //ID is 13 for mu-, -13 for mu+

  int parma_id = _parma.getIP();
  if (parma_id == 29) {
    pdgID = -13;
    const auto mass = 105.6583745 * MeV
  }
  else if (parma_id == 30) {
    pdgID = 13;
    const auto mass = 105.6583745 * MeV
  }
  else {
    std::cout << "Error: ParmaGenerator only supports muon+ or muon- at the moment. This run is invalid.";
  }

  const auto energy = particle[0] * MeV;
  const auto zenith_angle = particle[1] * rad;
  const auto azimuthal_angle = particle[2] * rad;

  const auto x0 = G4RandFlat::shoot(_x0_min, _x0_max);
  const auto y0 = G4RandFlat::shoot(_y0_min, _y0_max);
  const auto z0 = _z0;

  auto p_magnitude = std::sqrt(std::pow(energy, 2) - std::pow(mass, 2));

  auto sin_theta = std::sin(zenith_angle);

  auto p_x = p_magnitude * sin_theta * np.cos(azimuthal_angle);
  auto p_y = p_magnitude * sin_theta * np.sin(azimuthal_angle);
  auto p_z = p_magnitude * np.cos(zenith_angle);

  Particle out = new Particle(pdgID, x0, y0, z0, p_x, p_y, p_z);
  return out;
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////
//__Generate Initial Particles__________________________________________________________________
void ParmaGenerator::GeneratePrimaryVertex(G4Event* event) {
  ++_counter;
  double newParticle[3] = _parma.getParticle();
  Particle p = _convert_particle(newParticle);
  AddParticle(p, *event);
}
//----------------------------------------------------------------------------------------------

//__Get Last Event Data_________________________________________________________________________
GenParticleVector ParmaGenerator::GetLastEvent() const {
  return GenParticleVector{}; //empty. TODO: get actual last event data returned to fill into branches
}
//----------------------------------------------------------------------------------------------

void ParmaGenerator::SetParma(Parma4::Parma* parma) {
  if (!parma)
    return;
  _counter = 0ULL;
  _parma_settings->clear();
  _settings_on = false;
  _parma = _reconstruct_parma(parma);
  _parma->init();
}

void ParmaGenerator::SetParma(std::vector<std::string>& settings) {
  *_parma_settings = settings;
  _counter = 0ULL;
  _parma = _create_parma(_parma_settings, _settings_on);
}

void ParmaGenerator::SetParma() {
  _parma = _create_parma();
}

//__ParmaGenerator Specifications_______________________________________________________________
const Analysis::SimSettingList ParmaGenerator::GetSpecification() const {
  //THIS IS A PLACEHOLDER TAKEN FROM PythiaGenerator.cc!
  Analysis::SimSettingList config;
  if (_path.empty() && !_pythia_settings->empty()) {
    config = Analysis::IndexedSettings(SimSettingPrefix, "_SETTING_", *_pythia_settings);
  } else if (!_path.empty()) {
    config.emplace_back(SimSettingPrefix, "_CONFIG", _path);
  }

  Analysis::SimSettingList out;
  out.reserve(2UL + config.size());
  out.emplace_back(SimSettingPrefix, "", _name);
  out.emplace_back(SimSettingPrefix, "_FILTER",_filter->GetName());
  Analysis::SimSettingList filterspec=_filter->GetSpecification();
  out.insert(out.cend(), filterspec.begin(), filterspec.end());
  out.insert(out.cend(),
        std::make_move_iterator(config.begin()),
        std::make_move_iterator(config.end()));
  
  
  out.emplace_back(SimSettingPrefix, "_EVENTS", std::to_string(_counter));
  
  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace Physics */ //////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::MU */
