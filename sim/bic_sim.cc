// bic_sim: Geant4 model of the Hall D BIC test-beam setup (SFIL1-5, AstroPix module, Baby BCAL, trigger paddles)
//   bic_sim run.mac [seed]          batch
//   bic_sim                         interactive (vis.mac)
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4PhysListFactory.hh"
#include "Randomize.hh"
#include "DetectorConstruction.hh"
#include "Actions.hh"
#include <cstdlib>
int main(int argc, char **argv)
{
  if (argc > 2) CLHEP::HepRandom::setTheSeed(atol(argv[2]));
  G4RunManager *rm = new G4RunManager;
  DetectorConstruction *det = new DetectorConstruction; rm->SetUserInitialization(det);
  G4PhysListFactory factory; rm->SetUserInitialization(factory.GetReferencePhysList("FTFP_BERT_EMZ"));
  rm->SetUserInitialization(new ActionInitialization(det));
  G4VisManager *vis = new G4VisExecutive("quiet"); vis->Initialize();
  G4UImanager *ui = G4UImanager::GetUIpointer();
  if (argc > 1) ui->ApplyCommand(G4String("/control/execute ") + argv[1]);
  else { G4UIExecutive *ux = new G4UIExecutive(argc, argv); ui->ApplyCommand("/control/execute vis.mac"); ux->SessionStart(); delete ux; }
  delete vis; delete rm; return 0;
}
