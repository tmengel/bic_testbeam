// All dimensions of the Hall D BIC test-beam setup (August 2026, runs 3060-3064), settable from a
// macro before /run/initialize (see macros/geometry_default.mac). Beam along +z. x runs along the
// fibres from the North end (-x) to the South end (+x); y is vertical with cell 1 at the top.
#ifndef BIC_GEOMETRY_HH
#define BIC_GEOMETRY_HH
#include "G4SystemOfUnits.hh"
struct BicGeometry
{
  // Pb/SciFi matrix (ePIC BIC design values, GlueX-style: X0 ~ 1.45 cm)
  double fiberRadius = 0.5 * mm, claddingThickness = 0.06 * mm, fiberPitchY = 1.35 * mm, rowPitchZ = 1.22 * mm;   // epic BeamTest branch: cladding 0.06 mm, fibre spacing 1.35 x 1.22 mm
  int    rowsPerLayer = 17;                          // 17 rows -> 20.74 mm ~ 1.4 X0 per SFIL
  double edgeThickness = 0.61 * mm;                  // Pb/glue edge skin on each face of a Pb/ScFi layer (ePIC bic.xml, RadiatorEdgeThickness)
  // radiator: lead sheets + glue, ePIC material SciFiPb_PbGlue (96.2 % Pb + 3.8 % glue by mass, 8.558 g/cm3); edges SciFiPb_PbGlue_Edge (4.231 g/cm3)
  double radiatorDensity = 8.558 * g / cm3, edgeDensity = 4.231 * g / cm3, radiatorPbFraction = 0.9622, edgePbFraction = 0.805;
  // SFIL test articles (1-4 passive summing boards, 5 single SiPMs). ASSUMED: fibre length, cell height, spacing
  int    nSfil = 5, nCells = 4;
  double fiberLength = 510. * mm;                    // epic BeamTest Summer2026: SFIL lengths 519 / 505 / 518 / 505 / 505 mm (SFIL1-5)
  double cellHeight = 19.05 * mm;                    // epic BeamTest Summer2026: slab height 76.2 mm (3 in) read in 4 grids -> 4 cells of 19.05 mm
  double sfilGap = 43.1 * mm;                        // epic BeamTest Summer2026: SFIL centres 65.1 mm apart (SFIL4-5 61.9 mm) -> air gap 65.1 - 20.74 - 2 x 0.61 mm; was 17 mm (v2 assumption
  double carbonFrameThickness = 2. * mm, carbonStaveThickness = 0.5 * mm;   // ePIC: carbon frame plates of the imaging slot, carbon stave under the AstroPix module
  // AstroPix: one 9-chip module, 35 x 35 pixels of 0.5 mm per chip, chip pitch 20 mm (ePIC: 2 cm chips + margins)
  int    astroChips = 9, astroPixels = 35; double pixelPitch = 0.5 * mm, chipPitch = 20. * mm;
  int    astroSlot = 0;                              // 0 = in front of SFIL1 (confirmed by the analyst and by the 1-2 pixel clusters in the data); n = in the gap after SFIL n
  double astroX = 0., astroY = 0.;                   // module centre offsets; chip 0 is the North end
  double siliconThickness = 500. * um, electronicsThickness = 150. * um, copperThickness = 100. * um, kaptonThickness = 200. * um, epoxyThickness = 100. * um;
  // Baby BCAL (GlueX-type module, arXiv:2504.03079): 22.2 cm of Pb/ScFi (184 lead sheets), 58 cm long fibres, trapezoid 8.4-11.8 cm wide
  // (modelled as a 10.1 cm box), aluminium plates 0.8 cm upstream and 3.175 cm downstream, 40 light guides per end read as
  // 10 depth rows x 4 transverse columns of 4x4 SiPM arrays. Modelled as 10 depth segments of 18 fibre rows (21.96 cm) x 4 cells.
  int    bbcalLayers = 10, bbcalRowsPerLayer = 18, bbcalCells = 4; double bbcalCellHeight = 25.25 * mm, bbcalFiberLength = 580. * mm, bbcalGap = 61. * mm;   // epic BeamTest Summer2026: SFIL5 back face to the Baby BCAL front Al plate
  double bbcalFrontAl = 8. * mm, bbcalBackAl = 31.75 * mm;
  double bbcalYOffset = 0.;   // vertical offset of the whole Baby BCAL block relative to the SFIL stack centre line (unknown in the setup; scanned)
  // trigger paddles (R580 PMT counters) upstream; ASSUMED sizes
  double paddleX = 23.8 * mm, paddleY = 137. * mm, paddleZ = 14. * mm, paddleZpos = -180. * mm;   // epic BeamTest Summer2026 paddle slabs (P5/P3 ~180 mm, P6/P4 ~60 mm upstream of SFIL1)
  // layout: 1 = the August survey (context/MeasureMap2_Complete.pdf, the same numbers as the SFILn_X/Z/Length/Thickness, BBCal_offsetX and
  // paddle constants of the epic BeamTest file; default since 2026-10-29), 0 = v3 (equal SFIL pitch from sfilGap, everything centred at
  // x = 0, Baby BCAL after bbcalGap, AstroPix from astroSlot, paddles from the data-driven chip map, carbon frames in every slot)
  int    layout = 1;
  // survey, in cm: d = distance from the beam-pipe centre to the element's edge nearest the pipe (its South end; North is away from the
  // pipe), u = distance upstream from the Baby BCAL east face (back face of its downstream Al support) to the element's downstream face;
  // W = envelope thickness along z, L = length along x. Model frame: x = 0 at the AstroPix module centre (chip 4), North negative;
  // z = 0 at the SFIL1 upstream face, each SFIL centred in its measured envelope. The measured SFIL heights (H1-H5 = 11.43 / 11.27 /
  // 11.27 / 11.43 / 10.48 cm) are the outer envelope; the active Pb/ScFi stays 4 x cellHeight = 76.2 mm as in the epic file.
  static constexpr double svSfilD[5] = {25.718, 26.670, 27.781, 27.305, 28.575}, svSfilL[5] = {51.911, 50.483, 51.752, 50.483, 50.483};
  static constexpr double svSfilU[5] = {56.833, 50.324, 43.815, 37.306, 30.956}, svSfilW[5] = {2.2225, 2.3813, 2.2225, 2.2225, 2.3813};
  static constexpr double svBbcalD = 25.241, svAstroD = 41.910, svAstroU = 59.690;
  static constexpr double svPadD[4] = {49.371, 57.150, 39.688, 48.419}, svPadU[4] = {74.930, 64.929, 76.676, 64.770};   // P3 P4 P5 P6
  static constexpr double svPadX = 2.3813, svPadZ = 1.4287, svPadY = 14.922;
  // light model
  // light attenuation along the fibre: I(d) = alpha exp(-d/lambda1) + (1 - alpha) exp(-d/lambda2); the FTBF Baby BCAL analysis (epic BeamTest
  // branch, analysis/pimHistos.cpp, 'Sylvester's formula') uses alpha = 0.181, lambda1 = 60.9 cm, lambda2 = 418 cm; attenuationAlpha = 0 gives
  // the single exponential with attenuationLength (the 3000 mm used until 2026-10-24)
  double attenuationLength = 4182. * mm, attenuationLength1 = 608.7 * mm, attenuationAlpha = 0.181, vEff = 170. * mm / ns, birksConstant = 0.126 * mm / MeV, pixelThreshold = 1. * keV;
  // derived
  double layerThickness() const { return rowsPerLayer * rowPitchZ; }
  double sfilHeight() const { return nCells * cellHeight; }
  double bbcalHeight() const { return bbcalCells * bbcalCellHeight; }
  double bbcalLayerThickness() const { return bbcalRowsPerLayer * rowPitchZ; }
  double sfilFullThickness() const { return layerThickness() + 2 * edgeThickness; }
  double astroModuleLength() const { return astroChips * chipPitch; }
  // survey frame (cm in, mm out): x from the distance to the beam pipe, z from the distance upstream of the east face
  double svX(double d) const { return (svAstroD - d) * cm + astroModuleLength() / 2; }
  double svZ(double u) const { return (svSfilU[0] + svSfilW[0] / 2 - u) * cm + sfilFullThickness() / 2; }
  // per SFIL: fibre length, x centre, z of the upstream face (incl. the edge skin)
  double sfilLength(int L) const { return layout == 1 && L < 5 ? svSfilL[L] * cm : fiberLength; }
  double sfilXc(int L) const { return layout == 1 && L < 5 ? svX(svSfilD[L] + svSfilL[L] / 2) : 0.; }
  double sfilZfront(int L) const { return layout == 1 && L < 5 ? svZ(svSfilU[L] + svSfilW[L] / 2) - sfilFullThickness() / 2 : L * (sfilFullThickness() + sfilGap); }
  double bbcalXc() const { return layout == 1 ? svX(svBbcalD + bbcalFiberLength / cm / 2) : 0.; }
};
#endif
