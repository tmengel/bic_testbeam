#include "Actions.hh"
#include "DetectorConstruction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4GenericMessenger.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHandle.hh"
#include "G4NavigationHistory.hh"
#include "Randomize.hh"
#include "TFile.h"
#include "TTree.h"
#include "TNamed.h"
#include <cmath>

PrimaryGeneratorAction::PrimaryGeneratorAction(const DetectorConstruction *d) : det(d)
{
  gun = new G4ParticleGun(1);
  msg = new G4GenericMessenger(this, "/bic/gun/", "beam");
  msg->DeclareProperty("energy", energy, "kinetic energy in MeV"); msg->DeclareProperty("sigmaE", sigmaE, "gaussian energy spread in MeV");
  msg->DeclareProperty("x", x, "beam centre x in mm"); msg->DeclareProperty("y", y, "beam centre y in mm (cell 3 centre = -(2.5 - 2) x cellHeight... see README)"); msg->DeclareProperty("sigmaXY", sigmaXY, "gaussian spot size in mm");
  msg->DeclareProperty("z", z, "start z in mm"); msg->DeclareProperty("sigmaTheta", sigmaTheta, "gaussian divergence in mrad"); msg->DeclareProperty("particle", particle, "particle name");
  msg->DeclareProperty("cosmic", cosmic, "1 = cosmic mode: from a horizontal plane above the setup, cos^2 zenith, log-flat energy (see Actions.hh)");
  msg->DeclareProperty("cosmicY", cosmicY, "cosmic mode: height of the start plane in mm"); msg->DeclareProperty("cosmicX0", cosmicX0, "cosmic mode: x range of the start plane [mm]"); msg->DeclareProperty("cosmicX1", cosmicX1, "");
  msg->DeclareProperty("cosmicZ0", cosmicZ0, "cosmic mode: z range of the start plane [mm]"); msg->DeclareProperty("cosmicZ1", cosmicZ1, "");
  msg->DeclareProperty("cosmicEmin", cosmicEmin, "cosmic mode: kinetic energy range [MeV], log-flat"); msg->DeclareProperty("cosmicEmax", cosmicEmax, "");
  msg->DeclareProperty("cosmicTarget", cosmicTarget, "cosmic mode: 1 = every track passes through a trigger paddle (flux-correct face sampling, see Actions.hh)");
  msg->DeclareProperty("cosmicBack", cosmicBack, "cosmic target mode: distance [mm] before the paddle entry point at which the track starts");
}
PrimaryGeneratorAction::~PrimaryGeneratorAction() { delete gun; delete msg; }
void PrimaryGeneratorAction::GeneratePrimaries(G4Event *ev)
{
  gun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(particle));
  if (cosmic)
    {
      const double e = cosmicEmin * std::pow(cosmicEmax / cosmicEmin, G4UniformRand());
      const double ct = std::pow(G4UniformRand(), 0.25), st = std::sqrt(1 - ct * ct), ph = 2 * M_PI * G4UniformRand();
      gun->SetParticleEnergy(e * MeV);
      const G4ThreeVector dir(st * std::cos(ph), -ct, st * std::sin(ph));
      if (cosmicTarget && det && !det->paddleCentre.empty())
        {
          // entry face weights |d_i| x face area (the same for every paddle: equal sizes), so the paddle is drawn uniformly
          const G4ThreeVector h = det->paddleHalf; const int ip = std::min<int>(det->paddleCentre.size() - 1, int(G4UniformRand() * det->paddleCentre.size()));
          const double w[3] = {std::fabs(dir.x()) * 4 * h.y() * h.z(), std::fabs(dir.y()) * 4 * h.x() * h.z(), std::fabs(dir.z()) * 4 * h.x() * h.y()};
          double u = G4UniformRand() * (w[0] + w[1] + w[2]); int ax = u < w[0] ? 0 : (u < w[0] + w[1] ? 1 : 2);
          double q[3] = {(2 * G4UniformRand() - 1) * h.x(), (2 * G4UniformRand() - 1) * h.y(), (2 * G4UniformRand() - 1) * h.z()};
          const double hh[3] = {h.x(), h.y(), h.z()}, dd[3] = {dir.x(), dir.y(), dir.z()};
          q[ax] = dd[ax] > 0 ? -hh[ax] : hh[ax];                       // entry face: the side the track comes from
          G4ThreeVector entry = det->paddleCentre[ip] + G4ThreeVector(q[0], q[1], q[2]);
          double back = cosmicBack * mm; if (dir.y() < -1e-6) back = std::min(back, (entry.y() - cosmicY * mm) / dir.y());   // both negative: distance back to the plane y = cosmicY
          G4ThreeVector start = entry - back * dir;
          while ((std::fabs(start.x()) > 1.4 * m || std::fabs(start.y()) > 1.4 * m || std::fabs(start.z()) > 2.4 * m) && back > 1 * mm) { back /= 2; start = entry - back * dir; }
          gun->SetParticlePosition(start);
        }
      else gun->SetParticlePosition(G4ThreeVector((cosmicX0 + (cosmicX1 - cosmicX0) * G4UniformRand()) * mm, cosmicY * mm, (cosmicZ0 + (cosmicZ1 - cosmicZ0) * G4UniformRand()) * mm));
      gun->SetParticleMomentumDirection(dir);
      gun->GeneratePrimaryVertex(ev); return;
    }
  const double e = sigmaE > 0 ? G4RandGauss::shoot(energy, sigmaE) : energy;
  gun->SetParticleEnergy(e * MeV);
  const double px = x + (sigmaXY > 0 ? G4RandGauss::shoot(0, sigmaXY) : 0), py = y + (sigmaXY > 0 ? G4RandGauss::shoot(0, sigmaXY) : 0);
  gun->SetParticlePosition(G4ThreeVector(px * mm, py * mm, z * mm));
  const double tx = sigmaTheta > 0 ? G4RandGauss::shoot(0, sigmaTheta * 1e-3) : 0, ty = sigmaTheta > 0 ? G4RandGauss::shoot(0, sigmaTheta * 1e-3) : 0;
  gun->SetParticleMomentumDirection(G4ThreeVector(tx, ty, 1).unit());
  gun->GeneratePrimaryVertex(ev);
}

double RunAction::pixelThreshold() const { return det->geom().pixelThreshold; }
RunAction::RunAction(const DetectorConstruction *d) : det(d)
{
  msg = new G4GenericMessenger(this, "/bic/", "output");
  msg->DeclareProperty("output", outname, "output ROOT file");
  msg->DeclareProperty("stepsEvents", stepsEvents, "write every step of the first N events to the tree 'steps' (event display); 0 = off");
}
RunAction::~RunAction() { delete msg; }
void RunAction::BeginOfRunAction(const G4Run *)
{
  file = new TFile(outname.c_str(), "RECREATE");
  tree = new TTree("sim", "BIC test-beam simulation: one entry per beam particle; energies in MeV, times in ns, lengths in mm");
  EventRecord &r = rec;
  if (stepsEvents > 0)
    {
      steps = new TTree("steps", "steps of the first events (event display)");
      steps->Branch("event", &st_event, "event/I"); steps->Branch("x", &st_x); steps->Branch("y", &st_y); steps->Branch("z", &st_z); steps->Branch("edep", &st_e); steps->Branch("ke", &st_ke); steps->Branch("pdg", &st_pdg); steps->Branch("trk", &st_trk); steps->Branch("parent", &st_par);
    }
  tree->Branch("ebeam", &r.ebeam, "ebeam/D"); tree->Branch("x0", &r.x0, "x0/D"); tree->Branch("y0", &r.y0, "y0/D"); tree->Branch("z0", &r.z0, "z0/D"); tree->Branch("ux", &r.ux, "ux/D"); tree->Branch("uy", &r.uy, "uy/D"); tree->Branch("uz", &r.uz, "uz/D"); tree->Branch("pdg", &r.pdg, "pdg/I");
  tree->Branch("sfil_efib", r.sfil_efib, "sfil_efib[8][8]/D"); tree->Branch("sfil_evis", r.sfil_evis, "sfil_evis[8][8]/D"); tree->Branch("sfil_epb", r.sfil_epb, "sfil_epb[8][8]/D");
  tree->Branch("sfil_ln", r.sfil_ln, "sfil_ln[8][8]/D"); tree->Branch("sfil_ls", r.sfil_ls, "sfil_ls[8][8]/D"); tree->Branch("sfil_x", r.sfil_x, "sfil_x[8][8]/D"); tree->Branch("sfil_tn", r.sfil_tn, "sfil_tn[8][8]/D"); tree->Branch("sfil_ts", r.sfil_ts, "sfil_ts[8][8]/D");
  tree->Branch("sfil_e", r.sfil_e, "sfil_e[8]/D"); tree->Branch("sfil_etot", &r.sfil_etot, "sfil_etot/D");
  tree->Branch("bb_efib", r.bb_efib, "bb_efib[16][8]/D"); tree->Branch("bb_evis", r.bb_evis, "bb_evis[16][8]/D"); tree->Branch("bb_epb", r.bb_epb, "bb_epb[16][8]/D"); tree->Branch("bb_ln", r.bb_ln, "bb_ln[16][8]/D"); tree->Branch("bb_ls", r.bb_ls, "bb_ls[16][8]/D"); tree->Branch("bb_e", r.bb_e, "bb_e[16]/D"); tree->Branch("bb_etot", &r.bb_etot, "bb_etot/D");
  tree->Branch("paddle_e", r.paddle_e, "paddle_e[4]/D"); tree->Branch("e_escape", &r.e_escape, "e_escape/D"); tree->Branch("e_other", &r.e_other, "e_other/D"); tree->Branch("e_dep", &r.e_dep, "e_dep/D"); tree->Branch("astro_etot", &r.astro_etot, "astro_etot/D");
  tree->Branch("astro_n", &r.astro_n, "astro_n/I"); tree->Branch("astro_chip", &r.astro_chip); tree->Branch("astro_row", &r.astro_row); tree->Branch("astro_col", &r.astro_col); tree->Branch("astro_e", &r.astro_e); tree->Branch("astro_t", &r.astro_t);
}
void RunAction::EndOfRunAction(const G4Run *run)
{
  const BicGeometry &g = det->geom();
  TString info = Form("bic_sim: %d events. Geometry: %d SFILs, layer %.2f mm (%d rows, row pitch %.2f mm, fibre pitch %.2f mm, fibre r %.2f mm), %d cells x %.1f mm, fibre length %.0f mm, SFIL gap %.1f mm; AstroPix silicon at z = %.2f mm (slot %d), %d chips x %d pixels of %.2f mm; Baby BCAL at z = %.1f mm, %d segments x %d rows x %d cells x %.1f mm, fibres %.0f mm, Al %.1f/%.1f mm; radiator SciFiPb_PbGlue %.3f g/cm3 with %.2f mm edge skins; attenuation %.0f mm, v_eff %.1f mm/ns, Birks %.3f mm/MeV. SFIL1 front face at z = 0; x from North (-) to South (+); cell index 0 = top.",
                      run->GetNumberOfEvent(), g.nSfil, g.layerThickness(), g.rowsPerLayer, g.rowPitchZ, g.fiberPitchY, g.fiberRadius, g.nCells, g.cellHeight, g.fiberLength, g.sfilGap, det->astroZ, g.astroSlot, g.astroChips, g.astroPixels, g.pixelPitch, det->bbcalZ, g.bbcalLayers, g.bbcalRowsPerLayer, g.bbcalCells, g.bbcalCellHeight, g.bbcalFiberLength, g.bbcalFrontAl, g.bbcalBackAl, g.radiatorDensity / (CLHEP::g / CLHEP::cm3), g.edgeThickness, g.attenuationLength, g.vEff, g.birksConstant);
  info += g.layout == 1 ? Form(" Layout 1 = August survey (MeasureMap2_Complete.pdf): per-SFIL fibre lengths %.1f/%.1f/%.1f/%.1f/%.1f mm and x centres %.1f/%.1f/%.1f/%.1f/%.1f mm, SFIL z from the survey (the fibre length and SFIL gap above are not used), x = 0 at the AstroPix centre, Baby BCAL x centre %.1f mm.",
                                   g.sfilLength(0), g.sfilLength(1), g.sfilLength(2), g.sfilLength(3), g.sfilLength(4), g.sfilXc(0), g.sfilXc(1), g.sfilXc(2), g.sfilXc(3), g.sfilXc(4), g.bbcalXc())
                          : " Layout 0 = v3 (equal pitch, all centred at x = 0).";
  file->cd(); TNamed("info", info.Data()).Write(); if (steps) steps->Write(); tree->Write(); file->Close(); delete file; file = nullptr; tree = nullptr;
  G4cout << info << G4endl;
}

void EventAction::BeginOfEventAction(const G4Event *ev) { run->rec.clear(); run->stepsClear(); run->st_event = ev->GetEventID(); }
void EventAction::EndOfEventAction(const G4Event *ev)
{
  EventRecord &r = run->rec; const double thr = run->pixelThreshold();
  const G4PrimaryVertex *v = ev->GetPrimaryVertex(); if (v) { r.x0 = v->GetX0(); r.y0 = v->GetY0(); r.z0 = v->GetZ0(); const G4PrimaryParticle *p = v->GetPrimary(); r.ebeam = p->GetKineticEnergy(); r.pdg = p->GetPDGcode(); const G4ThreeVector u = p->GetMomentumDirection(); r.ux = u.x(); r.uy = u.y(); r.uz = u.z(); }
  for (int i = 0; i < EventRecord::NSFIL; i++) for (int j = 0; j < EventRecord::NCELL; j++) { if (r.sfil_evis[i][j] > 0) r.sfil_x[i][j] /= r.sfil_evis[i][j]; if (r.sfil_tn[i][j] > 1e8) r.sfil_tn[i][j] = -1; if (r.sfil_ts[i][j] > 1e8) r.sfil_ts[i][j] = -1; r.sfil_e[i] += r.sfil_efib[i][j] + r.sfil_epb[i][j]; }
  for (int i = 0; i < EventRecord::NSFIL; i++) r.sfil_etot += r.sfil_e[i];
  for (int i = 0; i < EventRecord::NBBL; i++) { for (int j = 0; j < EventRecord::NBBC; j++) r.bb_e[i] += r.bb_efib[i][j] + r.bb_epb[i][j]; r.bb_etot += r.bb_e[i]; }
  for (auto &kv : r.pixels) { if (kv.second.first < thr) continue; r.astro_chip.push_back(int(kv.first / 1000000)); r.astro_row.push_back(int((kv.first / 1000) % 1000)); r.astro_col.push_back(int(kv.first % 1000)); r.astro_e.push_back(kv.second.first); r.astro_t.push_back(kv.second.second); }
  r.astro_n = r.astro_e.size();
  if (run->steps && ev->GetEventID() < run->stepsEvents) run->steps->Fill();
  run->tree->Fill();
}

// Birks quenching on the step: E_vis = E / (1 + kB dE/dx), kB from the geometry parameters
static double birks(const G4Step *step, double kB)
{
  const double e = step->GetTotalEnergyDeposit(), l = step->GetStepLength();
  if (l <= 0 || step->GetTrack()->GetDefinition()->GetPDGCharge() == 0) return e;
  return e / (1. + kB * e / l);
}
SteppingAction::SteppingAction(RunAction *r, const DetectorConstruction *d) : run(r), det(d) {}
void SteppingAction::UserSteppingAction(const G4Step *step)
{
  EventRecord &r = run->rec; const BicGeometry &g = det->geom();
  const G4StepPoint *post = step->GetPostStepPoint();
  if (run->steps && run->st_event < run->stepsEvents)
    {
      const G4ThreeVector p = post->GetPosition(); const G4Track *tr = step->GetTrack();
      run->st_x.push_back(p.x() / mm); run->st_y.push_back(p.y() / mm); run->st_z.push_back(p.z() / mm); run->st_e.push_back(step->GetTotalEnergyDeposit() / MeV); run->st_ke.push_back(post->GetKineticEnergy() / MeV);
      run->st_pdg.push_back(tr->GetDefinition()->GetPDGEncoding()); run->st_trk.push_back(tr->GetTrackID()); run->st_par.push_back(tr->GetParentID());
    }
  if (post->GetStepStatus() == fWorldBoundary) { r.e_escape += post->GetKineticEnergy(); return; }
  const double edep = step->GetTotalEnergyDeposit(); if (edep <= 0) return;
  r.e_dep += edep;
  const G4TouchableHandle th = step->GetPreStepPoint()->GetTouchableHandle(); const G4String &name = th->GetVolume()->GetName();
  const G4ThreeVector pos = 0.5 * (step->GetPreStepPoint()->GetPosition() + post->GetPosition()); const double t = step->GetPreStepPoint()->GetGlobalTime();
  if (name == "sfil_core" || name == "bbcal_core")
    {
      const int copy = th->GetCopyNumber(1);   // the cladding placement carries det/layer/cell
      const int cell = (copy / 10000) % 10, dl = copy / 100000, det0 = dl / 100, layer = dl % 100;
      const double evis = birks(step, g.birksConstant);
      const double L = det0 == 0 ? g.sfilLength(layer) : g.bbcalFiberLength, xc = det0 == 0 ? g.sfilXc(layer) : g.bbcalXc();
      const double dN = pos.x() - xc + L / 2, dS = L / 2 - (pos.x() - xc);
      auto att = [&](double d) { return g.attenuationAlpha * std::exp(-d / g.attenuationLength1) + (1. - g.attenuationAlpha) * std::exp(-d / g.attenuationLength); };
      const double ln = evis * att(dN), ls = evis * att(dS);
      if (det0 == 0 && layer < EventRecord::NSFIL && cell < EventRecord::NCELL)
        {
          r.sfil_efib[layer][cell] += edep; r.sfil_evis[layer][cell] += evis; r.sfil_ln[layer][cell] += ln; r.sfil_ls[layer][cell] += ls; r.sfil_x[layer][cell] += evis * pos.x();
          const double tn = t + dN / g.vEff, ts = t + dS / g.vEff; if (tn < r.sfil_tn[layer][cell]) r.sfil_tn[layer][cell] = tn; if (ts < r.sfil_ts[layer][cell]) r.sfil_ts[layer][cell] = ts;
        }
      else if (det0 == 1 && layer < EventRecord::NBBL && cell < EventRecord::NBBC) { r.bb_efib[layer][cell] += edep; r.bb_evis[layer][cell] += evis; r.bb_ln[layer][cell] += ln; r.bb_ls[layer][cell] += ls; }
    }
  else if (name == "sfil_clad" || name == "bbcal_clad" || name == "sfil_pb" || name == "bbcal_pb")
    {
      int copy = th->GetCopyNumber(0); if (name == "sfil_clad" || name == "bbcal_clad") copy = copy / 100000; else copy = copy;   // pb copy = det*100 + layer
      const int det0 = copy / 100, layer = copy % 100;
      const double height = det0 == 0 ? g.sfilHeight() : g.bbcalHeight(), ch = det0 == 0 ? g.cellHeight : g.bbcalCellHeight;
      int cell = int((height / 2 - (pos.y() - (det0 == 1 ? g.bbcalYOffset : 0.))) / ch); if (cell < 0) cell = 0;
      if (det0 == 0 && layer < EventRecord::NSFIL) r.sfil_epb[layer][std::min(cell, EventRecord::NCELL - 1)] += edep;
      else if (det0 == 1 && layer < EventRecord::NBBL) r.bb_epb[layer][std::min(cell, EventRecord::NBBC - 1)] += edep;
    }
  else if (name == "astro_si")
    {
      const int chip = th->GetCopyNumber(0); const G4ThreeVector local = th->GetHistory()->GetTopTransform().TransformPoint(pos);
      const double active = g.astroPixels * g.pixelPitch;
      int col = int((local.x() + active / 2) / g.pixelPitch), row = int((local.y() + active / 2) / g.pixelPitch);
      col = std::max(0, std::min(g.astroPixels - 1, col)); row = std::max(0, std::min(g.astroPixels - 1, row));
      auto &px = r.pixels[chip * 1000000L + row * 1000L + col]; if (px.first == 0) px.second = t; px.first += edep; r.astro_etot += edep;
    }
  else if (name == "paddle") { r.paddle_e[th->GetCopyNumber(0)] += birks(step, g.birksConstant); }
  else r.e_other += edep;
}

void ActionInitialization::Build() const
{
  SetUserAction(new PrimaryGeneratorAction(det));
  RunAction *run = new RunAction(det); SetUserAction(run);
  SetUserAction(new EventAction(run));
  SetUserAction(new SteppingAction(run, det));
}
