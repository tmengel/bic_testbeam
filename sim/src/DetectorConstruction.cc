#include "DetectorConstruction.hh"
#include "G4GenericMessenger.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RotationMatrix.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include <cmath>

DetectorConstruction::DetectorConstruction()
{
  msg = new G4GenericMessenger(this, "/bic/geom/", "BIC test-beam geometry (set before /run/initialize)");
  msg->DeclareProperty("layout", geo.layout, "1 = August survey positions and lengths (default), 0 = v3 (equal pitch, centred, sfilGap / bbcalGap / astroSlot)");
  msg->DeclarePropertyWithUnit("fiberLength", "mm", geo.fiberLength, "SFIL fibre length (North-South), layout 0 only");
  msg->DeclarePropertyWithUnit("cellHeight", "mm", geo.cellHeight, "SFIL readout cell height");
  msg->DeclareProperty("nCells", geo.nCells, "SFIL cells per layer");
  msg->DeclareProperty("nSfil", geo.nSfil, "number of SFIL layers");
  msg->DeclareProperty("rowsPerLayer", geo.rowsPerLayer, "fibre rows per Pb/SciFi layer");
  msg->DeclarePropertyWithUnit("sfilGap", "mm", geo.sfilGap, "gap between SFILs");
  msg->DeclareProperty("astroSlot", geo.astroSlot, "AstroPix placed after SFIL n (0 = in front of SFIL1)");
  msg->DeclarePropertyWithUnit("astroX", "mm", geo.astroX, "AstroPix module centre x offset");
  msg->DeclarePropertyWithUnit("astroY", "mm", geo.astroY, "AstroPix module centre y offset");
  msg->DeclareProperty("bbcalLayers", geo.bbcalLayers, "Baby BCAL depth segments");
  msg->DeclareProperty("bbcalRowsPerLayer", geo.bbcalRowsPerLayer, "fibre rows per Baby BCAL depth segment");
  msg->DeclarePropertyWithUnit("edgeThickness", "mm", geo.edgeThickness, "Pb/glue edge skin per SFIL face");
  msg->DeclareProperty("bbcalCells", geo.bbcalCells, "Baby BCAL readout cells");
  msg->DeclarePropertyWithUnit("bbcalCellHeight", "mm", geo.bbcalCellHeight, "Baby BCAL cell height");
  msg->DeclarePropertyWithUnit("bbcalFiberLength", "mm", geo.bbcalFiberLength, "Baby BCAL fibre length");
  msg->DeclarePropertyWithUnit("attenuationLength1", "mm", geo.attenuationLength1, "short attenuation component");
  msg->DeclareProperty("attenuationAlpha", geo.attenuationAlpha, "weight of the short attenuation component (0 = single exponential)");
  msg->DeclarePropertyWithUnit("bbcalGap", "mm", geo.bbcalGap, "gap between SFIL5 and the Baby BCAL");
  msg->DeclarePropertyWithUnit("bbcalYOffset", "mm", geo.bbcalYOffset, "vertical offset of the Baby BCAL block");
  msg->DeclarePropertyWithUnit("paddleZ", "mm", geo.paddleZpos, "z of the trigger paddles");
  msg->DeclarePropertyWithUnit("attenuationLength", "mm", geo.attenuationLength, "fibre light attenuation length");
  msg->DeclareProperty("vEff", geo.vEff, "effective light speed in the fibre, mm/ns");
  msg->DeclarePropertyWithUnit("pixelThreshold", "keV", geo.pixelThreshold, "AstroPix pixel threshold on deposited energy");
}
DetectorConstruction::~DetectorConstruction() { delete msg; }

// one Pb/SciFi block: nlayers fibre layers along z, each rowsPerLayer rows; fibres along x; cells stacked in y (cell 0 at the top)
void DetectorConstruction::buildPbScFi(G4LogicalVolume *world, const char *name, int det, int nlayers, int rows, double fiberLength, int ncells, double cellHeight, double zfront, double layerSpacing, bool edges, double yoffset, double xc, int layer0)
{
  const double height = ncells * cellHeight, tl = rows * geo.rowPitchZ, te = edges ? geo.edgeThickness : 0;
  const int nfib = int(std::floor(height / geo.fiberPitchY));
  for (int l = 0; l < nlayers; l++)
    {
      const int L = layer0 + l;
      const double zc = zfront + te + l * (tl + 2 * te + layerSpacing) + tl / 2;
      if (det == 0) sfilZ.push_back(zc - tl / 2 - te); else bbcalLayerZ.push_back(zc - tl / 2 - te);
      if (edges)   // Pb/glue edge skins on both faces of the layer
        {
          G4Box *eb = new G4Box(G4String(name) + "_edge", fiberLength / 2, height / 2, te / 2); G4LogicalVolume *elv = new G4LogicalVolume(eb, edgemat, G4String(name) + "_edge_lv"); elv->SetVisAttributes(new G4VisAttributes(G4Colour(0.3, 0.3, 0.3, 0.3)));
          new G4PVPlacement(nullptr, G4ThreeVector(xc, yoffset, zc - tl / 2 - te / 2), elv, G4String(name) + "_edge", world, false, det * 100 + L);
          new G4PVPlacement(nullptr, G4ThreeVector(xc, yoffset, zc + tl / 2 + te / 2), elv, G4String(name) + "_edge", world, false, det * 100 + L);
        }
      G4Box *box = new G4Box(G4String(name) + "_pb", fiberLength / 2, height / 2, tl / 2);
      G4LogicalVolume *lv = new G4LogicalVolume(box, lead, G4String(name) + "_pb_lv");
      lv->SetVisAttributes(new G4VisAttributes(G4Colour(0.4, 0.4, 0.45, 0.3)));
      new G4PVPlacement(nullptr, G4ThreeVector(xc, yoffset, zc), lv, G4String(name) + "_pb", world, false, det * 100 + L);
      G4Tubs *ct = new G4Tubs(G4String(name) + "_clad", 0, geo.fiberRadius, fiberLength / 2, 0, 360 * deg);
      G4LogicalVolume *cladlv = new G4LogicalVolume(ct, clad, G4String(name) + "_clad_lv");
      G4Tubs *cr = new G4Tubs(G4String(name) + "_core", 0, geo.fiberRadius - geo.claddingThickness, fiberLength / 2, 0, 360 * deg);
      G4LogicalVolume *corelv = new G4LogicalVolume(cr, core, G4String(name) + "_core_lv");
      corelv->SetVisAttributes(new G4VisAttributes(G4Colour(0.2, 0.8, 0.9)));
      cladlv->SetVisAttributes(G4VisAttributes::GetInvisible());
      new G4PVPlacement(nullptr, G4ThreeVector(), corelv, G4String(name) + "_core", cladlv, false, 0);
      G4RotationMatrix *rot = new G4RotationMatrix(); rot->rotateY(90 * deg);   // tube axis along x
      for (int r = 0; r < rows; r++)
        {
          const double z = -tl / 2 + (r + 0.5) * geo.rowPitchZ; const double yoff = (r % 2) ? geo.fiberPitchY / 2 : 0;
          for (int k = 0; k < nfib; k++)
            {
              const double y = height / 2 - (k + 0.5) * geo.fiberPitchY - yoff;
              if (y - geo.fiberRadius < -height / 2) continue;
              const int cell = std::min(ncells - 1, int((height / 2 - y) / cellHeight));
              // copy number: det (0 SFIL, 1 Baby BCAL), layer, cell, then row x 100 + fibre index (< 10000: rows < 100, fibres per row < 100).
              // Until 2026-10-23 this was r * 1000 + k, which for rows 10-16 overflowed into the cell digit and booked 41 % of the fibres one cell too low.
              new G4PVPlacement(rot, G4ThreeVector(0, y, z), cladlv, G4String(name) + "_fib", lv, false, ((det * 100 + L) * 10 + cell) * 10000 + r * 100 + k);
            }
        }
    }
}

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  G4NistManager *nist = G4NistManager::Instance();
  air = nist->FindOrBuildMaterial("G4_AIR"); si = nist->FindOrBuildMaterial("G4_Si"); cu = nist->FindOrBuildMaterial("G4_Cu"); kapton = nist->FindOrBuildMaterial("G4_KAPTON"); al = nist->FindOrBuildMaterial("G4_Al");
  core = nist->BuildMaterialWithNewDensity("ScFiCore", "G4_POLYSTYRENE", 1.049 * g / cm3); core->GetIonisation()->SetBirksConstant(geo.birksConstant);
  clad = nist->FindOrBuildMaterial("G4_PLEXIGLASS"); scint = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"); scint->GetIonisation()->SetBirksConstant(geo.birksConstant);
  epoxy = new G4Material("Epoxy", 1.3 * g / cm3, 3); epoxy->AddElement(nist->FindOrBuildElement("C"), 15); epoxy->AddElement(nist->FindOrBuildElement("H"), 44); epoxy->AddElement(nist->FindOrBuildElement("O"), 7);
  // ePIC materials (compact/materials.xml): the radiator between the fibres is lead sheets + glue, the edge skins are glue-rich
  G4Material *glue = new G4Material("SciFiPb_Glue", 1.18 * g / cm3, 4); glue->AddElement(nist->FindOrBuildElement("H"), 0.091); glue->AddElement(nist->FindOrBuildElement("C"), 0.822); glue->AddElement(nist->FindOrBuildElement("N"), 0.032); glue->AddElement(nist->FindOrBuildElement("O"), 0.055);
  lead = new G4Material("SciFiPb_PbGlue", geo.radiatorDensity, 2); lead->AddMaterial(nist->FindOrBuildMaterial("G4_Pb"), geo.radiatorPbFraction); lead->AddMaterial(glue, 1. - geo.radiatorPbFraction);
  edgemat = new G4Material("SciFiPb_PbGlue_Edge", geo.edgeDensity, 2); edgemat->AddMaterial(nist->FindOrBuildMaterial("G4_Pb"), geo.edgePbFraction); edgemat->AddMaterial(glue, 1. - geo.edgePbFraction);
  carbon = new G4Material("CarbonFiber", 1.5 * g / cm3, 2); carbon->AddElement(nist->FindOrBuildElement("C"), 0.65); carbon->AddMaterial(epoxy, 0.35);

  G4Box *wb = new G4Box("world", 1.5 * m, 1.5 * m, 2.5 * m);
  G4LogicalVolume *world = new G4LogicalVolume(wb, air, "world_lv"); world->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4VPhysicalVolume *pw = new G4PVPlacement(nullptr, G4ThreeVector(), world, "world", nullptr, false, 0);

  const bool survey = geo.layout == 1;
  if (survey) { geo.paddleX = geo.svPadX * cm; geo.paddleY = geo.svPadY * cm; geo.paddleZ = geo.svPadZ * cm; }
  G4Box *pb = new G4Box("paddle", geo.paddleX / 2, geo.paddleY / 2, geo.paddleZ / 2);
  G4LogicalVolume *plv = new G4LogicalVolume(pb, scint, "paddle_lv"); plv->SetVisAttributes(new G4VisAttributes(G4Colour(1, 0.6, 0.1)));
  paddleCentre.clear(); paddleHalf = G4ThreeVector(geo.paddleX / 2, geo.paddleY / 2, geo.paddleZ / 2);
  if (survey)   // paddles where the survey puts them (P4 and P6 ~60 mm, P3 and P5 ~170 mm upstream of SFIL1), centred on the beam line in y
    for (int p = 0; p < 4; p++) { paddleCentre.push_back(G4ThreeVector(geo.svX(geo.svPadD[p] + geo.svPadX / 2), geo.astroY, geo.svZ(geo.svPadU[p] + geo.svPadZ / 2))); new G4PVPlacement(nullptr, paddleCentre.back(), plv, "paddle", world, false, p); }
  else
    {
      // v3: along x at the AstroPix chip groups they cover in the run-3062 data (P4 chip 0, P3 chips 2-3, P6 chips 5-6, P5 chip 8);
      // z: two planes as in the epic BeamTest file (P5 and P3 about 180 mm, P6 and P4 about 60 mm upstream of SFIL1)
      const double chipx[4] = {2.5, 0., 8., 5.5};   // P3 P4 P5 P6 in chip units
      const double pz[4] = {geo.paddleZpos, geo.paddleZpos + 120 * mm, geo.paddleZpos, geo.paddleZpos + 120 * mm};   // P3 P4 P5 P6
      for (int p = 0; p < 4; p++) { paddleCentre.push_back(G4ThreeVector(geo.astroX + (chipx[p] - 4) * geo.chipPitch, geo.astroY, pz[p] + (p % 2) * (geo.paddleZ + 1 * mm))); new G4PVPlacement(nullptr, paddleCentre.back(), plv, "paddle", world, false, p); }
    }

  // SFIL stack: SFIL1 front face at z = 0; survey: each SFIL with its own length, x centre and z
  sfilZ.clear(); bbcalLayerZ.clear();
  if (survey)
    for (int L = 0; L < geo.nSfil; L++) buildPbScFi(world, "sfil", 0, 1, geo.rowsPerLayer, geo.sfilLength(L), geo.nCells, geo.cellHeight, geo.sfilZfront(L), 0., true, 0., geo.sfilXc(L), L);
  else
    buildPbScFi(world, "sfil", 0, geo.nSfil, geo.rowsPerLayer, geo.fiberLength, geo.nCells, geo.cellHeight, 0., geo.sfilGap, true, 0.);
  // v3 only: carbon frame plates on both sides of every slot between SFILs (ePIC imaging-layer frame), and in the slot in front of SFIL1;
  // the survey and the epic BeamTest file have no frames between the SFILs
  if (!survey)
  {
    const double pitch = geo.layerThickness() + 2 * geo.edgeThickness + geo.sfilGap;
    G4Box *fb = new G4Box("cframe", geo.fiberLength / 2, geo.sfilHeight() / 2, geo.carbonFrameThickness / 2); G4LogicalVolume *flv = new G4LogicalVolume(fb, carbon, "cframe_lv"); flv->SetVisAttributes(new G4VisAttributes(G4Colour(0.2, 0.2, 0.2, 0.4)));
    for (int L = 0; L < geo.nSfil; L++)
      {
        const double zslotStart = L * pitch - geo.sfilGap;
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zslotStart + geo.carbonFrameThickness / 2), flv, "cframe", world, false, 2 * L);
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zslotStart + geo.sfilGap - geo.carbonFrameThickness / 2), flv, "cframe", world, false, 2 * L + 1);
      }
  }
  // AstroPix module in the gap after SFIL<astroSlot>
  const double tl = geo.layerThickness();
  const double stack = geo.siliconThickness + geo.electronicsThickness + geo.copperThickness + geo.kaptonThickness + geo.epoxyThickness + geo.carbonStaveThickness;
  const double tfull = tl + 2 * geo.edgeThickness;
  // survey: the module's downstream face (the stave) at the measured u, 6.3 mm upstream of SFIL1
  astroZ = survey ? geo.svZ(geo.svAstroU) - stack / 2 : (geo.astroSlot <= 0) ? -geo.sfilGap / 2 : geo.astroSlot * tfull + (geo.astroSlot - 1) * geo.sfilGap + geo.sfilGap / 2;
  {
    const double modLen = geo.astroChips * geo.chipPitch, active = geo.astroPixels * geo.pixelPitch;
    G4Box *sb = new G4Box("astro_si", active / 2, active / 2, geo.siliconThickness / 2);
    G4LogicalVolume *slv = new G4LogicalVolume(sb, si, "astro_si_lv"); slv->SetVisAttributes(new G4VisAttributes(G4Colour(0.9, 0.2, 0.2)));
    G4Box *eb = new G4Box("astro_el", geo.chipPitch / 2, geo.chipPitch / 2, geo.electronicsThickness / 2); G4LogicalVolume *elv = new G4LogicalVolume(eb, si, "astro_el_lv");
    G4Box *cb = new G4Box("astro_cu", modLen / 2, geo.chipPitch / 2, geo.copperThickness / 2); G4LogicalVolume *clv = new G4LogicalVolume(cb, cu, "astro_cu_lv");
    G4Box *kb = new G4Box("astro_kapton", modLen / 2, geo.chipPitch / 2, geo.kaptonThickness / 2); G4LogicalVolume *klv = new G4LogicalVolume(kb, kapton, "astro_kapton_lv");
    G4Box *xb = new G4Box("astro_epoxy", modLen / 2, geo.chipPitch / 2, geo.epoxyThickness / 2); G4LogicalVolume *xlv = new G4LogicalVolume(xb, epoxy, "astro_epoxy_lv");
    double z = astroZ - stack / 2;
    for (int c = 0; c < geo.astroChips; c++)   // chip 0 at the North end (-x)
      {
        const double x = geo.astroX + (c - (geo.astroChips - 1) / 2.) * geo.chipPitch;
        new G4PVPlacement(nullptr, G4ThreeVector(x, geo.astroY, z + geo.siliconThickness / 2), slv, "astro_si", world, false, c);
        new G4PVPlacement(nullptr, G4ThreeVector(x, geo.astroY, z + geo.siliconThickness + geo.electronicsThickness / 2), elv, "astro_el", world, false, c);
      }
    z += geo.siliconThickness + geo.electronicsThickness;
    new G4PVPlacement(nullptr, G4ThreeVector(geo.astroX, geo.astroY, z + geo.copperThickness / 2), clv, "astro_cu", world, false, 0); z += geo.copperThickness;
    new G4PVPlacement(nullptr, G4ThreeVector(geo.astroX, geo.astroY, z + geo.kaptonThickness / 2), klv, "astro_kapton", world, false, 0); z += geo.kaptonThickness;
    new G4PVPlacement(nullptr, G4ThreeVector(geo.astroX, geo.astroY, z + geo.epoxyThickness / 2), xlv, "astro_epoxy", world, false, 0); z += geo.epoxyThickness;
    G4Box *stb = new G4Box("astro_stave", modLen / 2, geo.chipPitch / 2, geo.carbonStaveThickness / 2); G4LogicalVolume *stlv = new G4LogicalVolume(stb, carbon, "astro_stave_lv");
    new G4PVPlacement(nullptr, G4ThreeVector(geo.astroX, geo.astroY, z + geo.carbonStaveThickness / 2), stlv, "astro_stave", world, false, 0);
  }
  // Baby BCAL downstream; survey: its east face (back of the downstream Al support, the survey benchmark) at u = 0
  const double bbDepth = geo.bbcalFrontAl + geo.bbcalLayers * geo.bbcalLayerThickness() + geo.bbcalBackAl;
  bbcalZ = survey ? geo.svZ(0.) - bbDepth : geo.nSfil * tfull + (geo.nSfil - 1) * geo.sfilGap + geo.bbcalGap;
  {
    const double xb = geo.bbcalXc();
    G4Box *ab = new G4Box("bbcal_alfront", geo.bbcalFiberLength / 2, geo.bbcalHeight() / 2, geo.bbcalFrontAl / 2); G4LogicalVolume *alv = new G4LogicalVolume(ab, al, "bbcal_alfront_lv");
    new G4PVPlacement(nullptr, G4ThreeVector(xb, geo.bbcalYOffset, bbcalZ + geo.bbcalFrontAl / 2), alv, "bbcal_al", world, false, 0);
    const double zpb = bbcalZ + geo.bbcalFrontAl;
    buildPbScFi(world, "bbcal", 1, geo.bbcalLayers, geo.bbcalRowsPerLayer, geo.bbcalFiberLength, geo.bbcalCells, geo.bbcalCellHeight, zpb, 0., false, geo.bbcalYOffset, xb);
    const double zback = zpb + geo.bbcalLayers * geo.bbcalLayerThickness();
    G4Box *bb = new G4Box("bbcal_alback", geo.bbcalFiberLength / 2, geo.bbcalHeight() / 2, geo.bbcalBackAl / 2); G4LogicalVolume *blv = new G4LogicalVolume(bb, al, "bbcal_alback_lv");
    new G4PVPlacement(nullptr, G4ThreeVector(xb, geo.bbcalYOffset, zback + geo.bbcalBackAl / 2), blv, "bbcal_al", world, false, 1);
  }
  if (survey)
    {
      G4cout << "survey layout (MeasureMap2_Complete.pdf), x from North (-) to South (+), x = 0 at the AstroPix centre:" << G4endl;
      for (int L = 0; L < geo.nSfil; L++) G4cout << "  SFIL" << L + 1 << ": x " << geo.sfilXc(L) - geo.sfilLength(L) / 2 << " .. " << geo.sfilXc(L) + geo.sfilLength(L) / 2 << " mm, z " << sfilZ[L] << " .. " << sfilZ[L] + tfull << " mm" << G4endl;
      G4cout << "  AstroPix: x " << geo.astroX - geo.astroModuleLength() / 2 << " .. " << geo.astroX + geo.astroModuleLength() / 2 << " mm, z " << astroZ - stack / 2 << " .. " << astroZ + stack / 2 << " mm" << G4endl;
      const char *pn[4] = {"P3", "P4", "P5", "P6"};
      for (int p = 0; p < 4; p++) G4cout << "  " << pn[p] << ": x centre " << geo.svX(geo.svPadD[p] + geo.svPadX / 2) << " mm, z centre " << geo.svZ(geo.svPadU[p] + geo.svPadZ / 2) << " mm" << G4endl;
      G4cout << "  Baby BCAL: x " << geo.bbcalXc() - geo.bbcalFiberLength / 2 << " .. " << geo.bbcalXc() + geo.bbcalFiberLength / 2 << " mm, z " << bbcalZ << " .. " << bbcalZ + bbDepth << " mm (gap to SFIL5 " << bbcalZ - sfilZ[geo.nSfil - 1] - tfull << " mm)" << G4endl;
    }

  G4cout << "BIC test-beam geometry: " << geo.nSfil << " SFILs of " << tl / mm << " mm (" << geo.rowsPerLayer << " rows) + 2 x " << geo.edgeThickness / mm << " mm edges, radiator " << lead->GetDensity() / (g / cm3) << " g/cm3, cells " << geo.nCells << " x " << geo.cellHeight / mm << " mm, fibres " << geo.fiberLength / mm << " mm; AstroPix silicon at z = " << astroZ / mm << " mm; Baby BCAL front at z = " << bbcalZ / mm << " mm, " << geo.bbcalLayers << " segments x " << geo.bbcalRowsPerLayer << " rows x " << geo.bbcalCells << " cells, Al " << geo.bbcalFrontAl / mm << " / " << geo.bbcalBackAl / mm << " mm" << G4endl;
  return pw;
}
