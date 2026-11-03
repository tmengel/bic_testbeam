#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH
#include "G4VUserDetectorConstruction.hh"
#include "BicGeometry.hh"
#include <vector>
#include "G4ThreeVector.hh"
class G4GenericMessenger; class G4LogicalVolume; class G4Material;
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction(); ~DetectorConstruction() override;
  G4VPhysicalVolume *Construct() override;
  const BicGeometry &geom() const { return geo; }
  // z of the front face of each SFIL / the Baby BCAL / the AstroPix silicon, filled by Construct()
  std::vector<double> sfilZ, bbcalLayerZ; double astroZ = 0, bbcalZ = 0;
  // trigger paddles as placed (centres and half sizes), for the cosmic target mode of the generator
  std::vector<G4ThreeVector> paddleCentre; G4ThreeVector paddleHalf;
private:
  void buildPbScFi(G4LogicalVolume *world, const char *name, int det, int nlayers, int rows, double fiberLength, int ncells, double cellHeight, double zfront, double layerSpacing, bool edges, double yoffset = 0., double xc = 0., int layer0 = 0);
  BicGeometry geo; G4GenericMessenger *msg = nullptr;
  G4Material *lead = nullptr, *edgemat = nullptr, *carbon = nullptr, *al = nullptr, *core = nullptr, *clad = nullptr, *air = nullptr, *scint = nullptr, *si = nullptr, *cu = nullptr, *kapton = nullptr, *epoxy = nullptr;
};
#endif
