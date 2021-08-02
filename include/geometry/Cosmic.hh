

#ifndef MU__GEOMETRY_COSMIC_HH
#define MU__GEOMETRY_COSMIC_HH
#pragma once
#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "analysis.hh"
#include "geometry/Construction.hh"

namespace MATHUSLA { namespace MU {

namespace Cosmic { ////////////////////////////////////////////////////////////////////////////////

class Scintillator {
public:
  Scintillator(const std::string& name,
               const double length,
               const double height,
               const double width,
               const double thickness);

  struct Material {
    static G4Material* Casing;
    static G4Material* Scintillator;
    static void Define();
  private:
    Material();
  };

  double GetLength() const { return _length; }
  double GetHeight() const { return _height; }
  double GetWidth() const { return _width; }
  double GetCasingThickness() const { return _thickness; }
  G4LogicalVolume* GetVolume() const { return _lvolume; }
  G4VPhysicalVolume* GetSensitiveVolume() const { return _sensitive; }

  void Register(G4VSensitiveDetector* detector);

  G4VPhysicalVolume* PlaceIn(G4LogicalVolume* parent,
                             const G4Transform3D& transform);

private:
  G4LogicalVolume* _lvolume;
  G4VPhysicalVolume* _sensitive;
  std::string _name;
  double _length, _height, _width, _thickness;
};

class Detector : public G4VSensitiveDetector {
public:
  Detector();

  void Initialize(G4HCofThisEvent* event);
  G4bool ProcessHits(G4Step* step, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

  static const bool DataPerEvent = true;
  static const std::string& DataName;
  static const Analysis::ROOT::DataKeyList DataKeys;
  static const Analysis::ROOT::DataKeyTypeList DataKeyTypes;
  static TTree* pre_data;

  static G4VPhysicalVolume* Construct(G4LogicalVolume* world);
  static G4VPhysicalVolume* ConstructEarth(G4LogicalVolume* world);
  static G4VPhysicalVolume* ConstructModule(G4LogicalVolume* detector, int tag_number, double detector_x, double detector_y, double detector_z);
  static G4VPhysicalVolume* ConstructScintillatorLayer(G4LogicalVolume* Module_volume, int module_number, int layer_number, double module_x_displacement, double module_y_displacement, double layer_z_displacement);
  static bool SaveAll;
  static bool SaveCut;

  static void WritePreData();
};

} /* namespace Cosmic */ //////////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::MU */

#endif /* MU__GEOMETRY_COSMIC_HH */
