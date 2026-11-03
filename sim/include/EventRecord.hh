// what one simulated event records; filled by SteppingAction, written by EventAction/RunAction
#ifndef EVENT_RECORD_HH
#define EVENT_RECORD_HH
#include <vector>
#include <map>
#include <cstring>
struct EventRecord
{
  static const int NSFIL = 8, NCELL = 8, NBBL = 16, NBBC = 8, NPAD = 4;
  double ebeam = 0, x0 = 0, y0 = 0, z0 = 0, ux = 0, uy = 0, uz = 1; int pdg = 0;   // primary vertex and direction (cosmic mode: from above)
  double sfil_efib[NSFIL][NCELL], sfil_evis[NSFIL][NCELL], sfil_epb[NSFIL][NCELL], sfil_ln[NSFIL][NCELL], sfil_ls[NSFIL][NCELL], sfil_x[NSFIL][NCELL], sfil_tn[NSFIL][NCELL], sfil_ts[NSFIL][NCELL], sfil_e[NSFIL];
  double bb_efib[NBBL][NBBC], bb_evis[NBBL][NBBC], bb_epb[NBBL][NBBC], bb_ln[NBBL][NBBC], bb_ls[NBBL][NBBC], bb_e[NBBL], bb_etot;
  double paddle_e[NPAD], e_escape, e_other, e_dep, sfil_etot, astro_etot; int astro_n;
  std::vector<int> astro_chip, astro_row, astro_col; std::vector<double> astro_e, astro_t;
  std::map<long, std::pair<double, double>> pixels;   // key chip*1000000 + row*1000 + col -> (edep, first time)
  void clear()
  {
    memset(sfil_efib, 0, sizeof(sfil_efib)); memset(sfil_evis, 0, sizeof(sfil_evis)); memset(sfil_epb, 0, sizeof(sfil_epb)); memset(sfil_ln, 0, sizeof(sfil_ln)); memset(sfil_ls, 0, sizeof(sfil_ls)); memset(sfil_x, 0, sizeof(sfil_x)); memset(sfil_e, 0, sizeof(sfil_e));
    for (int i = 0; i < NSFIL; i++) for (int j = 0; j < NCELL; j++) sfil_tn[i][j] = sfil_ts[i][j] = 1e9;
    memset(bb_efib, 0, sizeof(bb_efib)); memset(bb_evis, 0, sizeof(bb_evis)); memset(bb_epb, 0, sizeof(bb_epb)); memset(bb_ln, 0, sizeof(bb_ln)); memset(bb_ls, 0, sizeof(bb_ls)); memset(bb_e, 0, sizeof(bb_e)); memset(paddle_e, 0, sizeof(paddle_e));
    bb_etot = e_escape = e_other = e_dep = sfil_etot = astro_etot = 0; astro_n = 0; astro_chip.clear(); astro_row.clear(); astro_col.clear(); astro_e.clear(); astro_t.clear(); pixels.clear();
  }
};
#endif
