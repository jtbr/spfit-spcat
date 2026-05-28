#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include "spinv_internal.h"
#include "SpinvEngine.hpp"
#include "api/InputSchema.hpp"
#include "common/CalError.hpp"

SpinvEngine::SpinvEngine()
{
  // Zero out the entire context first
  memset(&m_context, 0, sizeof(SpinvContext));

  // Apply specific initializations
  m_context.sptzero[0] = 1;
  m_context.sptzero[1] = 0;
  m_context.ssp_head.next = NULL; // ssp_head itself needs init
  m_context.ssp_head.sspt = m_context.sptzero;
  m_context.ssp_head.ssize = 1;
  m_context.ssp_head.nitot = 0;
  m_context.revsym[0] = 0;
  m_context.revsym[1] = 3;
  m_context.revsym[2] = 2;
  m_context.revsym[3] = 1;
  m_context.isoddk[0] = 0;
  m_context.isoddk[1] = 1;
  m_context.isoddk[2] = 1;
  m_context.isoddk[3] = 0;
  m_context.ixphase[0] = 0;
  m_context.ixphase[1] = 1;
  m_context.ixphase[2] = 2;
  m_context.ixphase[3] = 3;
  m_context.ipwr2[0] = 0;
  m_context.ipwr2[1] = 1;
  m_context.ipwr2[2] = 2;
  m_context.ipwr2[3] = 4;

  // Initialize pointers to static fallback members
  m_context.vinfo = &m_context.vinfo1;        // vinfo1 itself needs full init
  memset(&m_context.vinfo1, 0, sizeof(SVIB)); // Zero vinfo1
  m_context.vinfo1.spt = m_context.sptzero;   // Then set its specific pointers

  m_context.dipinfo = &m_context.dipinfo0; // dipinfo0 needs full init
  memset(&m_context.dipinfo0, 0, sizeof(SDIP));

  // Pointing arrays to their single static "zero" element is okay for initial state
  // as setopt/setblk will realloc them.
  m_context.moldv = &m_context.zmoldv;
  m_context.blkptr = &m_context.zblkptr;
  m_context.ipder = &m_context.zipder;
  m_context.ivs = &m_context.zivs;
  m_context.wk = &m_context.zwk;
  m_context.idx = &m_context.zidx;
  m_context.jdx = &m_context.zjdx;
  m_context.iqnsep = &m_context.ziqnsep;
  m_context.ibkptr = &m_context.zibkptr;
  m_context.ikmin = &m_context.zikmin;

  // The single static elements themselves
  m_context.zwk = 0.0;
  m_context.zmoldv = 0;
  m_context.zblkptr = 0;
  m_context.zivs = 0;
  m_context.zipder = 0; // ipder values are significant; -1 is often "not fitted"
  m_context.zidx = 0;
  m_context.zjdx = 0;
  m_context.ziqnsep = 0;
  m_context.zibkptr = 0;
  m_context.zikmin = 0;

  // Initialize glob explicitly too (setopt overwrites many but not all)
  memset(&m_context.glob, 0, sizeof(GLOB));
  // setopt then sets: glob.lsym=TRUE, glob.esym=TRUE, glob.esymdec=100 etc.
}

SpinvEngine::~SpinvEngine()
{
  // Free any dynamically allocated memory in the context
  if (m_context.vinfo != &m_context.vinfo1) {
      free(m_context.vinfo);
  }
  if (m_context.moldv != &m_context.zmoldv) {
      free(m_context.moldv);
  }
  if (m_context.blkptr != &m_context.zblkptr) {
      free(m_context.blkptr);
  }
  if (m_context.ipder != &m_context.zipder) {
      free(m_context.ipder);
  }
  if (m_context.ivs != &m_context.zivs) {
      free(m_context.ivs);
  }
  if (m_context.wk != &m_context.zwk) {
      free(m_context.wk);
  }
  if (m_context.idx != &m_context.zidx) {
      free(m_context.idx);
  }
  if (m_context.jdx != &m_context.zjdx) {
      free(m_context.jdx);
  }
  if (m_context.iqnsep != &m_context.ziqnsep) {
      free(m_context.iqnsep);
  }
  if (m_context.ibkptr != &m_context.zibkptr) {
      free(m_context.ibkptr);
  }
  if (m_context.ikmin != &m_context.zikmin) {
      free(m_context.ikmin);
  }
}

int SpinvEngine::hamx(int iblk, int nsize, int npar, bcd_t *idpar, double *par,
                      double *egy, double *t, double *dedp,
                      double *pmix, const BOOL ifdump)
{
    return ::hamx(&m_context, iblk, nsize, npar, idpar, par, egy, t, dedp, pmix, ifdump);
}

int SpinvEngine::setint(FILE *lu, BOOL *ifdiag, int *nsav, const int ndip,
                        bcd_t *idip, int *isimag)
{
    return ::setint(&m_context, lu, ifdiag, nsav, ndip, idip, isimag);
}

int SpinvEngine::intens(const int iblk, const int isiz, const int jblk,
                           const int jsiz, const int ndip, const bcd_t *idip,
                           const double *dip, double *s)
{
    return ::intens(&m_context, iblk, isiz, jblk, jsiz, ndip, idip, dip, s);
}

int SpinvEngine::getqn(const int iblk, const int indx, const int maxqn,
                        short *iqn, int *idgn)
{
    return ::getqn(&m_context, iblk, indx, maxqn, iqn, idgn);
}

int SpinvEngine::setopt(FILE *lu, int *nfmt, int *itd, int *ndbcd, char *namfil)
{
    return ::setopt(&m_context, lu, nfmt, itd, ndbcd, namfil);
}

// Encode nuclear spins as a decimal string that getbcd will parse into bcdspin.
// getbcd packs digits right-to-left into BCD nibbles, so spins are concatenated
// in reverse order: for spins {s0, s1, ..., sn-1} write "s_{n-1}...s1 s0".
// Each spin value 2*I must be in [1, 9] (the isbig=0 single-nibble format).
// Zero spins outputs "0" which triggers isbig=1 with jj=0, terminating immediately.
static std::string encode_spinv_spin(const VibState &vib)
{
    for (int s : vib.spin_degeneracies) {
        if (s < 1 || s > 9)
            throw ValidationError("spin_degeneracy value must be in [1,9]; "
                                  "use isbig BCD format for larger spins (not yet supported)",
                                  CalErrorCode::InvalidParameter);
    }
    std::string result;
    if (vib.symmetric_rotor_quanta)
        result += '-';
    if (vib.spin_degeneracies.empty()) {
        result += '0';
    } else {
        for (int i = (int)vib.spin_degeneracies.size() - 1; i >= 0; --i)
            result += (char)('0' + vib.spin_degeneracies[i]);
    }
    return result;
}

// Build the text cards that SpinvEngine::setopt can parse from a fmemopen buffer.
// Card format (getbcd reads the leading integer, pcard reads rvec[0..10] after it):
//   <spin_int> <lopt> <knnmin> <knnmax> <ixz> <iax> <iwtpl> <iwtmn> <vsym> <ewt0> <idiag> <phase>
// - Card 1: lopt = ±nvib (sign = oblate flag), ixz/idiag/phase from SpinvOptions globals
// - Card i>1: lopt = vibs[i].index (0-based vib state index to configure)
// - Intermediate cards (not last) use vsym = -1.0 to signal "keep reading"
// - Last card uses vibs.back().vsym (must be >= -0.5; clamped if negative)
static std::string make_spinv_cards(const SpinvOptions &so)
{
    if (so.vibs.empty())
        throw ValidationError("SpinvOptions.vibs must have at least one entry",
                              CalErrorCode::InvalidParameter);

    std::ostringstream oss;
    int nvib = (int)so.vibs.size();

    // Trailing all-default VibStates arise when the original .par had fewer cards
    // than nvib (legacy_parser fills the remainder with VibState() defaults).
    // Emitting a card for them would switch setopt from its single-card path
    // (which handles multi-vib weight distribution via vsym) to the multi-card
    // path, breaking the weight setup for molecules like ch3oh.
    auto vib_is_all_default = [](const VibState &v) {
        return v.spin_degeneracies.empty()
            && v.knmin == 0
            && v.knmax == 359
            && v.iwtpl == 1
            && v.iwtmn == 1
            && std::fabs(v.vsym) < 0.5
            && v.esym_weight == 99
            && v.stat_weight_axis == 1
            && !v.symmetric_rotor_quanta;
    };
    int last_card = nvib - 1;
    while (last_card > 0 && vib_is_all_default(so.vibs[last_card]))
        --last_card;

    for (int i = 0; i <= last_card; ++i) {
        const VibState &vib = so.vibs[i];
        bool is_last = (i == last_card);

        std::string spin = encode_spinv_spin(vib);

        // lopt: card 1 → ±nvib (total vib states, not cards); subsequent cards → vib state index
        int lopt;
        if (i == 0)
            lopt = so.oblate ? -nvib : nvib;
        else
            lopt = i;

        // vsym: intermediate cards must be < -0.5 (loop-continue signal).
        // Last card: user value clamped to >= -0.4 so it doesn't continue the loop.
        double vsym;
        if (!is_last) {
            vsym = -1.0;
        } else {
            vsym = (vib.vsym < -0.5) ? -0.4 : vib.vsym;
        }

        // card 1 carries global options; later cards carry zeros for those fields
        int ixx   = (i == 0) ? so.inclusion_flags : 0;
        int idiag = (i == 0) ? so.diag_order     : 0;
        int phase = (i == 0) ? so.phase_flags     : 0;

        oss << spin << " " << lopt
            << " " << vib.knmin << " " << vib.knmax
            << " " << ixx << " " << vib.stat_weight_axis
            << " " << vib.iwtpl << " " << vib.iwtmn
            << " " << vsym << " " << vib.esym_weight
            << " " << idiag << " " << phase << "\n";
    }
    return oss.str();
}

int SpinvEngine::apply_options(const EngineOptions &opts, int *nfmt, int *itd, int *ndbcd,
                               std::string &namfil)
{
    std::string cards = make_spinv_cards(opts.spinv);
    FILE *f = fmemopen((void *)cards.data(), cards.size(), "r");
    if (!f)
        throw IoError("fmemopen failed in SpinvEngine::apply_options");

    char namfil_buf[256] = {};
    if (!namfil.empty()) {
        strncpy(namfil_buf, namfil.c_str(), sizeof(namfil_buf) - 1);
    }

    int n = ::setopt(&m_context, f, nfmt, itd, ndbcd, namfil_buf);
    fclose(f);

    // nam_file overrides the engine default derived from the card's leading letter
    if (!opts.spinv.nam_file.empty())
        namfil = opts.spinv.nam_file;
    else
        namfil = namfil_buf;

    return n;
}

int SpinvEngine::setfmt(int *iqnfmt, int nfmt)
{
    return ::setfmt(&m_context, iqnfmt, nfmt);
}

int SpinvEngine::setblk(FILE *lu, int npar, bcd_t *idpar, double *par,
                        int *nblkpf, int *negy)
{
    return ::setblk(&m_context, lu, npar, idpar, par, nblkpf, negy);
}