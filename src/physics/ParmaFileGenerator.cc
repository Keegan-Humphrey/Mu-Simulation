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

#include "physics/Generator.hh"

#include <ostream>

#include <Randomize.hh>
#include <tls.hh>

#include "physics/Units.hh"

#include <iostream>

namespace MATHUSLA { namespace MU {

namespace Physics { ////////////////////////////////////////////////////////////////////////////

//__Parma Generator Constructor_________________________________________________________________
ParmaGenerator::ParmaGenerator(const std::string& name,
                               const std::string& description,
                               const Particle& particle)
    : ParmaGenerator(name, description, particle, particle) {}
//----------------------------------------------------------------------------------------------

//__Parma Generator Constructor_________________________________________________________________
ParmaGenerator::ParmaGenerator(const std::string& name,
                               const std::string& description,
                               const Particle& min,
                               const Particle& max)
    : Generator(name, description, min), _min(min), _max(max) {
  GenerateCommands();
}
//----------------------------------------------------------------------------------------------

//__Generate UI Commands________________________________________________________________________
void ParmaGenerator::GenerateCommands() {
  _ui_pT_min = CreateCommand<Command::DoubleUnitArg>("pT_min", "Set Minimum Transverse Momentum.");
  _ui_pT_min->SetParameterName("pT_min", false, false);
  _ui_pT_min->SetRange("pT_min >= 0");
  _ui_pT_min->SetDefaultUnit("GeV/c");
  _ui_pT_min->SetUnitCandidates("eV/c keV/c MeV/c GeV/c");
  _ui_pT_min->AvailableForStates(G4State_PreInit, G4State_Idle);

  _ui_pT_max = CreateCommand<Command::DoubleUnitArg>("pT_max", "Set Maximum Transverse Momentum.");
  _ui_pT_max->SetParameterName("pT_max", false, false);
  _ui_pT_max->SetRange("pT_max >= 0");
  _ui_pT_max->SetDefaultUnit("GeV/c");
  _ui_pT_max->SetUnitCandidates("eV/c keV/c MeV/c GeV/c");
  _ui_pT_max->AvailableForStates(G4State_PreInit, G4State_Idle);

  _ui_eta_min = CreateCommand<Command::DoubleArg>("eta_min", "Set Minimum Pseudorapidity.");
  _ui_eta_min->SetParameterName("eta_min", false);
  _ui_eta_min->AvailableForStates(G4State_PreInit, G4State_Idle);

  _ui_eta_max = CreateCommand<Command::DoubleArg>("eta_max", "Set Maximum Pseudorapidity.");
  _ui_eta_max->SetParameterName("eta_max", false);
  _ui_eta_max->AvailableForStates(G4State_PreInit, G4State_Idle);

  _ui_phi_min = CreateCommand<Command::DoubleUnitArg>("phi_min", "Set Minimum Semi-Opening Angle.");
  _ui_phi_min->SetParameterName("phi_min", false, false);
  _ui_phi_min->SetDefaultUnit("deg");
  _ui_phi_min->SetUnitCandidates("degree deg radian rad milliradian mrad");
  _ui_phi_min->AvailableForStates(G4State_PreInit, G4State_Idle);

  _ui_phi_max = CreateCommand<Command::DoubleUnitArg>("phi_max", "Set Maximum Semi-Opening Angle.");
  _ui_phi_max->SetParameterName("phi_max", false, false);
  _ui_phi_max->SetDefaultUnit("deg");
  _ui_phi_max->SetUnitCandidates("degree deg radian rad milliradian mrad");
  _ui_phi_max->AvailableForStates(G4State_PreInit, G4State_Idle);

  _ui_ke_min = CreateCommand<Command::DoubleUnitArg>("ke_min", "Set Minimum Kinetic Energy.");
  _ui_ke_min->SetParameterName("ke_min", false, false);
  _ui_ke_min->SetRange("ke_min >= 0");
  _ui_ke_min->SetDefaultUnit("GeV");
  _ui_ke_min->SetUnitCandidates("eV keV MeV GeV");
  _ui_ke_min->AvailableForStates(G4State_PreInit, G4State_Idle);

  _ui_ke_max = CreateCommand<Command::DoubleUnitArg>("ke_max", "Set Maximum Kinetic Energy.");
  _ui_ke_max->SetParameterName("ke_max", false, false);
  _ui_ke_max->SetRange("ke_max >= 0");
  _ui_ke_max->SetDefaultUnit("GeV");
  _ui_ke_max->SetUnitCandidates("eV keV MeV GeV");
  _ui_ke_max->AvailableForStates(G4State_PreInit, G4State_Idle);
}
//----------------------------------------------------------------------------------------------

//__Generate Initial Particles__________________________________________________________________
void ParmaGenerator::GeneratePrimaryVertex(G4Event* event) {
  std::size_t particle_parameters_index = _particle_parameters.size();
  particle_parameters_index = _event_counter;

  ifstream pfile ("MuonP.txt");
  if (pfile.is_open()) {
    
  }

  ++_event_counter;
  AddParticle(_particle_parameters.at(particle_parameters_index), *event);
}
//----------------------------------------------------------------------------------------------

//__Get Last Event Data_________________________________________________________________________
GenParticleVector ParmaGenerator::GetLastEvent() const {
  return GenParticleVector{};
}
//----------------------------------------------------------------------------------------------

//__Parma Generator Messenger Set Value_________________________________________________________
// void ParmaGenerator::SetNewValue(G4UIcommand* command, G4String value) {
//   if (command == _ui_pathname) {
//     std::ifstream input_stream(value);
//     while (input_stream) {
//       const auto next_char = input_stream.peek();
//       if (next_char == std::ifstream::traits_type::eof()) {
//               break;
//       }
//       if (next_char == ' ' || next_char == '\t' || next_char == '\r' || next_char == '\n') {
//               input_stream.ignore();
//               continue;
//       }
//       if (next_char == '#') {
//               input_stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
//               continue;
//       }
//       _particle_parameters.emplace_back();
//       auto &new_parameters = _particle_parameters.back();
//       if ( ! (input_stream >> new_parameters.id
//                            >> new_parameters.x
//                            >> new_parameters.y
//                            >> new_parameters.z
//                            >> new_parameters.px
//                            >> new_parameters.py
//                            >> new_parameters.pz)) {
//         throw std::runtime_error("Unable to parse particle parameters file");
//       }
//     }
//     if ( ! input_stream) {
//       throw std::runtime_error("Unable to read particle parameters file");
//     }
//   } else {
//     Generator::SetNewValue(command, value);
//   }
// }
//----------------------------------------------------------------------------------------------

//__Parma Generator Information String__________________________________________________________
std::ostream& ParmaGenerator::Print(std::ostream& os) const {
  os << "Generator Info:\n  "
     << "Name:        " << _name        << "\n  "
     << "Description: " << _description << "\n  "
     << "Particle ID: " << _particle.id << "\n  ";

  if (_using_range_ke) {
    os << "avg ke:      " << G4BestUnit(0.5 * (_min.ke() + _max.ke()), "Energy") << "\n    "
       << "ke min:      " << G4BestUnit(_min.ke(), "Energy")                     << "\n    "
       << "ke max:      " << G4BestUnit(_max.ke(), "Energy")                     << "\n  ";
  } else {
    os << "avg pT:      " << G4BestUnit(0.5 * (_min.pT() + _max.pT()), "Momentum") << "\n    "
       << "pT min:      " << G4BestUnit(_min.pT(), "Momentum")                     << "\n    "
       << "pT max:      " << G4BestUnit(_max.pT(), "Momentum")                     << "\n  ";
  }

  os << "avg eta:     "  << 0.5 * (_min.eta() + _max.eta())                      << "\n    "
     << "eta min:     "  << _min.eta()                                           << "\n    "
     << "eta max:     "  << _max.eta()                                           << "\n  "
     << "avg phi:     "  << G4BestUnit(0.5 * (_min.phi() + _max.phi()), "Angle") << "\n    "
     << "phi min:     "  << G4BestUnit(_min.phi(), "Angle")                      << "\n    "
     << "phi max:     "  << G4BestUnit(_max.phi(), "Angle")                      << "\n  "
     << "vertex:      (" << G4BestUnit(_particle.t, "Time")                      << ", "
                         << G4BestUnit(_particle.x, "Length")                    << ", "
                         << G4BestUnit(_particle.y, "Length")                    << ", "
                         << G4BestUnit(_particle.z, "Length")                    << ")\n";
  return os;
}
//----------------------------------------------------------------------------------------------

//__ParmaGenerator Specifications_______________________________________________________________
const Analysis::SimSettingList ParmaGenerator::GetSpecification() const {
  return Analysis::Settings(SimSettingPrefix,
    "",         _name,
    "_PDG_ID",  std::to_string(_particle.id),
    (_using_range_ke ? "_KE_MIN" : "_PT_MIN"),
    (_using_range_ke ? std::to_string(_min.ke() / Units::Energy)   + " " + Units::EnergyString
                     : std::to_string(_min.pT() / Units::Momentum) + " " + Units::MomentumString),
    (_using_range_ke ? "_KE_MAX" : "_PT_MAX"),
    (_using_range_ke ? std::to_string(_max.ke() / Units::Energy)   + " " + Units::EnergyString
                     : std::to_string(_max.pT() / Units::Momentum) + " " + Units::MomentumString),
    "_P_MAG_MIN", std::to_string(_min.p_mag() / Units::Momentum) + " " + Units::MomentumString,
    "_P_MAG_MAX", std::to_string(_max.p_mag() / Units::Momentum) + " " + Units::MomentumString,
    "_ETA_MIN", std::to_string(_min.eta()),
    "_ETA_MAX", std::to_string(_max.eta()),
    "_PHI_MIN", std::to_string(_min.phi() / Units::Angle) + " " + Units::AngleString,
    "_PHI_MAX", std::to_string(_max.phi() / Units::Angle) + " " + Units::AngleString,
    "_VERTEX", "(" + std::to_string(_particle.t / Units::Time)   + ", "
                   + std::to_string(_particle.x / Units::Length) + ", "
                   + std::to_string(_particle.y / Units::Length) + ", "
                   + std::to_string(_particle.z / Units::Length) + ")");
}
//----------------------------------------------------------------------------------------------

} /* namespace Physics */ //////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::MU */
