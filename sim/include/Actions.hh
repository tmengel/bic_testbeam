#ifndef ACTIONS_HH
#define ACTIONS_HH
#include "G4VUserActionInitialization.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include <vector>
#include "EventRecord.hh"
#include <string>
class G4ParticleGun; class G4GenericMessenger; class TFile; class TTree; class DetectorConstruction;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction(const DetectorConstruction *d = nullptr); ~PrimaryGeneratorAction() override;
  void GeneratePrimaries(G4Event *) override;
  double energy = 4000., sigmaE = 0., x = 0., y = 0., sigmaXY = 5., z = -500., sigmaTheta = 0.; std::string particle = "e+";
  // cosmic mode (2026-10-30): particles start on a horizontal plane at y = cosmicY, uniform in x [cosmicX0, cosmicX1] and z [cosmicZ0, cosmicZ1],
  // zenith angle with the sea-level cos^2 distribution weighted by the flux through a horizontal plane (cos theta = u^(1/4)), azimuth uniform,
  // kinetic energy log-flat between cosmicEmin and cosmicEmax (the MIP response depends only weakly on it); particle from /bic/gun/particle (mu-)
  int cosmic = 0; double cosmicY = 250., cosmicX0 = -400., cosmicX1 = 400., cosmicZ0 = -350., cosmicZ1 = 750., cosmicEmin = 500., cosmicEmax = 20000.;
  // cosmic target mode (cosmicTarget 1, with cosmic 1): the same zenith / azimuth / energy distributions, but every track passes through one of
  // the trigger paddles as placed by the detector (flux-correct: a paddle is chosen with probability proportional to its projected area for the
  // drawn direction, the entry face by its projected area, the entry point uniformly on that face); the track starts cosmicBack mm before the
  // entry point (or at the plane y = cosmicY if that is closer). This is the population the paddle-triggered cosmic runs select, which a
  // start plane of practical size misses for the near-horizontal tracks that cross a paddle AND the stack.
  int cosmicTarget = 0; double cosmicBack = 800.;
private:
  G4ParticleGun *gun; G4GenericMessenger *msg; const DetectorConstruction *det;
};

class RunAction : public G4UserRunAction
{
public:
  RunAction(const DetectorConstruction *det); ~RunAction() override;
  void BeginOfRunAction(const G4Run *) override; void EndOfRunAction(const G4Run *) override;
  EventRecord rec; TFile *file = nullptr; TTree *tree = nullptr; std::string outname = "bic_sim.root";
  // event display: the steps of the first stepsEvents events go to a second tree 'steps'
  int stepsEvents = 0; TTree *steps = nullptr; int st_event = 0;
  std::vector<float> st_x, st_y, st_z, st_e, st_ke; std::vector<int> st_pdg, st_trk, st_par;
  void stepsClear() { st_x.clear(); st_y.clear(); st_z.clear(); st_e.clear(); st_ke.clear(); st_pdg.clear(); st_trk.clear(); st_par.clear(); }
  double pixelThreshold() const;
private:
  const DetectorConstruction *det; G4GenericMessenger *msg;
};

class EventAction : public G4UserEventAction
{
public:
  EventAction(RunAction *r) : run(r) {}
  void BeginOfEventAction(const G4Event *) override; void EndOfEventAction(const G4Event *) override;
  RunAction *run;
};

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction(RunAction *r, const DetectorConstruction *d);
  void UserSteppingAction(const G4Step *) override;
private:
  RunAction *run; const DetectorConstruction *det;
};

class ActionInitialization : public G4VUserActionInitialization
{
public:
  ActionInitialization(const DetectorConstruction *d) : det(d) {}
  void Build() const override;
private:
  const DetectorConstruction *det;
};
#endif
