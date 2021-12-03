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
#include <Randomize.hh>
#include <tls.hh>

#include "physics/ParmaGenerator.hh"

#include "geometry/Earth.hh"
#include "geometry/Cavern.hh"
#include "physics/Units.hh"
#include "util/string.hh"

#include <sys/time.h>
#include <unistd.h>

namespace MATHUSLA { namespace MU {

namespace Physics { ////////////////////////////////////////////////////////////////////////////

//__Parma Generator Constructor_________________________________________________________________
ParmaGenerator::ParmaGenerator()
  : Generator("parma", "Parma4 Cosmic Spectra Generator.") {
  GenerateCommands();
  SetParma();
}
//__Parma Generator Constructor_________________________________________________________________
// ParmaGenerator::ParmaGenerator(Parma4::Parma* parma) : ParmaGenerator() {
//   _parma_settings = new std::vector<std::string>();
//   SetParma(parma);
//   GenerateCommands();
// }
//----------------------------------------------------------------------------------------------

//__Parma Generator Constructor_________________________________________________________________
// ParmaGenerator::ParmaGenerator(std::vector<std::string>& settings)
//   : ParmaGenerator() {
//   SetParma(settings);
//   GenerateCommands();
// }
//----------------------------------------------------------------------------------------------


void ParmaGenerator::GenerateCommands() {
  // _ui_id = CreateCommand<Command::IntegerArg>("id", "Set Particle Id.");
  // _ui_id->SetParameterName("id", false);
  // _ui_id->AvailableForStates(G4State_PreInit, G4State_Idle);

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
    auto out = new Parma4::Parma(parma->getIP(), parma->getSeed());
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
Parma4::Parma* _create_parma(int ip, int seed) {
  auto parma = new Parma4::Parma();
  parma->init(ip, seed);
  return parma;
}
//----------------------------------------------------------------------------------------------

//__Convert Parma output array to Geant particle with specified x0,y0,z0 bounds________________
Particle _convert_particle(const double *particle, const double *bounds) {  
  //particle comes in as {ip, energy, zenith, azimuthal}
  //bounds comes in as {x0 min, x0 max, y0 min, y0 max, z0}
  int pdgID;  //ID is 13 for mu-, -13 for mu+

  int parma_id = particle[0];
  auto mass = 0 * MeV;
  if (parma_id == 29) {
    pdgID = -13;
    mass = 105.6583745 * MeV;
  }
  else if (parma_id == 30) {
    pdgID = 13;
    mass = 105.6583745 * MeV;
  }
  else {
    std::cout << "Error: ParmaGenerator only supports muon+ or muon- at the moment. This run is invalid.";
    mass = 0;
  }

  const auto energy = particle[1] * MeV;
  const auto zenith_angle = particle[2] * rad;
  const auto azimuthal_angle = particle[3] * rad;

  const auto x0 = G4RandFlat::shoot(bounds[0], bounds[1]);
  const auto y0 = G4RandFlat::shoot(bounds[2], bounds[3]);
  const auto z0 = bounds[4];

  auto p_magnitude = std::sqrt(std::pow(energy, 2) - std::pow(mass, 2));

  auto sin_theta = std::sin(zenith_angle);

  auto p_x = p_magnitude * sin_theta * std::cos(azimuthal_angle);
  auto p_y = p_magnitude * sin_theta * std::sin(azimuthal_angle);
  auto p_z = p_magnitude * std::cos(zenith_angle);
  
  Particle out(pdgID, x0, y0, z0, p_x, p_y, p_z);
  return out;
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////
//__Generate Initial Particles__________________________________________________________________
void ParmaGenerator::GeneratePrimaryVertex(G4Event* event) {
  ++_counter; //loop over until find interesting particle
  std::cout << "\nSVE1 Making Parma Vertex\n";
  Particle p;
  while (true) {
    const double *newParticle = _parma->getParticle();
    p = _convert_particle(newParticle, _bounds);

    G4ThreeVector p_direction = p.p_unit();
    auto p_x = p_direction.getX();
    auto p_y = p_direction.getY();
    auto p_z = p_direction.getZ();
    auto x0 = p.x;
    auto y0 = p.y;
    auto z0 = p.z;
    auto s_floor = (_z_floor - z0) / p_z;
    auto xint = x0 + s_floor * p_x;
    auto yint = y0 + s_floor * p_y; // Calculating x and y intercepts with floor.
    if(-_box_side_length <= xint && xint <= _box_side_length) { // Checking if floor int in box
      if(-_box_side_length <= yint && yint <= _box_side_length) { // NOTE: assuming symm. abt 0!!
        break;
      }
    }
  }
  //_last_event gets stored in the ntuple
  _last_event.clear();
  _last_event.push_back(GenParticle(p));
  std::cout << "last event size: " << _last_event.size() << "\n";
  AddParticle(p, *event);
}
//----------------------------------------------------------------------------------------------

//__Get Last Event Data_________________________________________________________________________
GenParticleVector ParmaGenerator::GetLastEvent() const {
  return _last_event; //empty. TODO: get actual last event data returned to fill into branches
}
//----------------------------------------------------------------------------------------------

void ParmaGenerator::SetNewValue(G4UIcommand* command,
                                 G4String value) {
  if (command == _ui_id) {
    int uiid = _ui_id->GetNewIntValue(value);
    if (uiid == 13) {
      _parma->setIP(29);
    }
    if (uiid == -13) {
      _parma->setIP(30);
    }
    else {
      std::cout << "ParmaGenerator: The particle ID set was not a muon.\n";
      std::cout << "ParmaGenerator: Using muon+ instead.\n";
      _parma->setIP(29);
    }
  }
  else if (command == _ui_x0_min) {
    _bounds[0] = _ui_x0_min->GetNewDoubleValue(value);
  } else if (command == _ui_x0_max) {
    _bounds[1] = _ui_x0_max->GetNewDoubleValue(value);
  } else if (command == _ui_y0_min) {
    _bounds[2] = _ui_y0_min->GetNewDoubleValue(value);
  } else if (command == _ui_y0_max) {
    _bounds[3] = _ui_y0_max->GetNewDoubleValue(value);
  } else if (command == _ui_z0) {
    _bounds[4] = _ui_z0->GetNewDoubleValue(value);
  } else {
    Generator::SetNewValue(command, value);
  }
}

// void ParmaGenerator::SetParma(Parma4::Parma* parma) {
//   if (!parma)
//     return;
//   _counter = 0ULL;
//   _parma_settings->clear();
//   _settings_on = false;
//   _parma = _reconstruct_parma(parma);
//   _parma->init();
// }

// void ParmaGenerator::SetParma(std::vector<std::string>& settings) {
//   *_parma_settings = settings;
//   _counter = 0ULL;
//   _parma = _create_parma(_parma_settings, _settings_on);
// }

void ParmaGenerator::SetParma() {
  int parma_ip = 29; // default value of parma_ip for muon+
  struct timeval curTime;
  gettimeofday(&curTime, NULL);
  long int micro_sec = curTime.tv_usec;
  std::cout << "ParmaGenerator: seed = " << micro_sec << "\n";
  _parma = _create_parma(parma_ip, micro_sec);
}

//__ParmaGenerator Specifications_______________________________________________________________
const Analysis::SimSettingList ParmaGenerator::GetSpecification() const {
  Analysis::SimSettingList out;
  
  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace Physics */ //////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::MU */
