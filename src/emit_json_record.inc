/** @file emit_json_record.inc
 *
 *   JSON record emitter for --emit-json mode. This is #included at the END of
 *   lce_misc.cpp (NOT compiled standalone) because it calls the StatHandler<SH>
 *   template methods (computeStatNow, and the trait/handler accessors) whose
 *   bodies live in stathandler.cpp / stat_rec_base.h, which lce_misc.cpp already
 *   pulls in. That guarantees the concrete TTNeutralSH / TTQuantiSH / MetapopSH
 *   instantiations are visible.
 *
 *   Output is a HEADER + FRAMES stream (schema_version 2):
 *     - ONE `kind:"header"` record is emitted first, declaring the static layout
 *       of the run (patch count, the locus table in canonical order with the
 *       per-locus allele counts, the trait count, and which layer-1 estimators
 *       this run computes).
 *     - then ONE compact `kind:"frame"` record per logged generation, carrying
 *       only positional numeric payloads (allele_freqs, phenotype, stats) read
 *       against the header. No keys, no repeated identity.
 *   A reader must consume the header before it can interpret any frame.
 *
 *   STATS (the metapop-wide scalar layer-1 estimators) are NO LONGER a hardcoded
 *   five. We enumerate the stat recorders the run actually requested (parsed from
 *   the .ini `stat` tokens into StatRecorder objects) and, for each that is a
 *   metapop-wide SCALAR, emit its catalog `stream_key` into header.stats and its
 *   value into frame.stats at the matching position, every logged generation.
 *
 *   The recorder->stream_key mapping is keyed on the recorder's NAME, which the
 *   engine builds to be exactly the .ini `stat` token (e.g. "n.adlt.fst.wc",
 *   "adlt.nbInd", "q.qst"); see ttneutral/ttquanti/metapop setStatRecorders +
 *   stathandler.cpp add(). The token->stream_key table below is derived from
 *   contracts/stats-catalog.json (catalog_version 2) — that catalog is the
 *   SOURCE OF TRUTH; keep this table in sync with it. Only stats with
 *   output:"scalar" (i.e. a `stream_key`) appear here. Per-patch/per-locus/
 *   pairwise recorders carry a name suffix (_p1, _l1, _pair, ...) so they never
 *   match a bare token key and are skipped. A requested scalar recorder with no
 *   catalog mapping is skipped gracefully (no key, no value, no crash).
 *
 *   See emit_json.h for the public API. Distributed under GPL v3 (or later),
 *   like the rest of quantiNemo.
 */

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "emit_json.h"
#include "patch.h"
#include "stat_rec_base.h"
#include "stathandlerbase.h"
#include "tindividual.h"
#include "ttrait.h"
#include "ttneutral.h"
#include "ttquanti.h"
#include "types.h"

namespace {

// Emit a double as a JSON number, or `null` for the engine's my_NAN sentinel /
// non-finite values. A value we cannot reach is null, never a fake number.
void ej_put_num_or_null(std::ostream& os, double v)
{
    if (!std::isfinite(v) || v == (double)my_NAN) { os << "null"; return; }
    std::ostringstream tmp;
    tmp.precision(6);
    tmp << v;
    os << tmp.str();
}

void ej_put_num(std::ostream& os, double v)
{
    std::ostringstream tmp;
    tmp.precision(6);
    tmp << v;
    os << tmp.str();
}

// canonical locus / genotype ordering: all neutral loci first, then quanti
const char* const EJ_LOCUS_ORDER[2] = { "ntrl", "quanti" };

unsigned int ej_count_quanti_traits(TMetapop* pop)
{
    unsigned int n = 0, nbProto = pop->getTraitPrototypeSize();
    for (unsigned int t = 0; t < nbProto; ++t)
        if (pop->getTraitPrototype(t).get_type() == "quanti") ++n;
    return n;
}

// --------------------------------------------------------------------------- //
// token -> stream_key map for the metapop-wide SCALAR stats.
//
// DERIVED FROM contracts/stats-catalog.json (catalog_version 2): every stat with
// output:"scalar" and a `stream_key`. The KEY is the catalog `token`, which the
// engine uses verbatim as the StatRecorder NAME for the scalar form of that stat
// (per-patch/per-locus/pairwise forms append a _p/_l/_pair suffix to the name, so
// they never match here). Keep this in sync with the catalog (the source of
// truth). If a requested scalar recorder is NOT in this table it is skipped.
// --------------------------------------------------------------------------- //
struct EjStatMap { const char* token; const char* stream_key; };
const EjStatMap EJ_STAT_MAP[] = {
    // --- NEU (neutral-marker genetic estimators) ---
    { "n.adlt.fst.wc",    "fst"        },
    { "n.adlt.fst",       "fst_nei"    },
    { "n.adlt.fis.wc",    "fis"        },
    { "n.adlt.fit.wc",    "fit"        },
    { "n.adlt.ho",        "ho"         },
    { "n.adlt.hs",        "hs"         },
    { "n.adlt.ht",        "ht"         },
    { "n.adlt.hs.p2",     "hs_p2"      },
    { "n.adlt.ht.p2",     "ht_p2"      },
    { "n.adlt.nbAll",     "nb_all"     },
    { "n.adlt.nbAll.tot", "nb_all_tot" },
    { "n.adlt.nbFixLoc",  "nb_fix_loc" },
    { "n.adlt.rs",        "rs"         },
    { "n.adlt.rt",        "rt"         },
    { "n.adlt.r",         "a_range"    },
    { "n.adlt.gw",        "gw"         },
    { "n.adlt.chord",     "chord"      },
    // --- KIN (kinship / coancestry) ---
    { "n.adlt.theta",     "theta"      },
    { "n.adlt.alpha",     "alpha"      },
    { "n.adlt.coa.fsib",  "coa_fsib"   },
    { "adlt.fsib",        "fsib_prop"  },
    // --- QT (quantitative-trait estimators) ---
    { "q.qst",            "qst"        },
    { "q.VgB",            "vg_b"       },
    { "q.VgW",            "vg_w"       },
    { "q.VpB",            "vp_b"       },
    { "q.VpW",            "vp_w"       },
    { "q.VaW",            "va_w"       },
    // --- DEM (demographic / census) ---
    { "adlt.nbInd",       "nb_ind"     },
    { "adlt.meanInd",     "mean_ind"   },
    { "adlt.sexRatio",    "sex_ratio"  },
    { "adlt.nbPops",      "nb_pops"    },
    { "ext.rate",         "ext_rate"   },
    // --- DISP (dispersal / migration) ---
    { "emigrants",        "emigrants"  },
    { "immigrate",        "immigrate"  },
    { "colonisers",       "colonisers" },
    // --- FEC (fecundity) ---
    { "fem.meanFec",      "mean_fec_fem" },
    // --- FIT (fitness) ---
    { "VwW",              "vw_w"       },
    // --- COAL (coalescent) ---
    { "mrca.mean",        "mrca_mean"  },
};
const unsigned int EJ_STAT_MAP_N = sizeof(EJ_STAT_MAP) / sizeof(EJ_STAT_MAP[0]);

// stream_key for a recorder name, or NULL if the name is not a catalog scalar.
const char* ej_stream_key_for(const std::string& name)
{
    for (unsigned int k = 0; k < EJ_STAT_MAP_N; ++k)
        if (name == EJ_STAT_MAP[k].token) return EJ_STAT_MAP[k].stream_key;
    return 0;
}

// One streamed scalar stat: its catalog stream_key, plus the handler + recorder
// needed to (re)compute its current value for a frame.
struct EjScalarStat {
    std::string        key;
    StatHandlerBase*   handler;
    StatRecBase*       rec;
};

// Scan one handler's recorders; append the ones that are catalog scalars.
void ej_collect_from_handler(StatHandlerBase* sh, std::vector<EjScalarStat>& out)
{
    if (!sh) return;
    std::list<StatRecBase*>& recs = sh->getStats();
    for (std::list<StatRecBase*>::iterator it = recs.begin(); it != recs.end(); ++it) {
        StatRecBase* r = *it;
        if (!r) continue;
        const char* sk = ej_stream_key_for(r->getName());
        if (!sk) continue;                 // not a streamable scalar: skip
        EjScalarStat s;
        s.key = sk; s.handler = sh; s.rec = r;
        out.push_back(s);
    }
}

// Build the ordered list of streamed scalar stats this run computes. Order is
// STABLE (neutral handlers, then quanti handlers, then the metapop handler, in
// prototype/recorder registration order) so header.stats and frame.stats stay
// positionally aligned across the whole stream.
std::vector<EjScalarStat> ej_collect_scalar_stats(TMetapop* pop)
{
    std::vector<EjScalarStat> out;
    unsigned int nbProto = pop->getTraitPrototypeSize();

    // neutral trait handlers first, then quanti, matching the locus/stat ordering
    for (int pass = 0; pass < 2; ++pass) {
        for (unsigned int t = 0; t < nbProto; ++t) {
            TTraitProto& proto = pop->getTraitPrototype(t);
            if (proto.get_type() != EJ_LOCUS_ORDER[pass]) continue;
            if (proto.get_type() == "ntrl") {
                TTraitNeutralProto* np = dynamic_cast<TTraitNeutralProto*>(&proto);
                if (np && np->_stats) ej_collect_from_handler(np->_stats, out);
            } else { // quanti
                TTraitQuantiProto* qp = dynamic_cast<TTraitQuantiProto*>(&proto);
                if (qp && qp->_stats) ej_collect_from_handler(qp->_stats, out);
            }
        }
    }

    // demographic / dispersal / fecundity / fitness / kinship metapop scalars
    ej_collect_from_handler(pop->get_popStatHandler(), out);

    return out;
}

// ---- header: static run layout, emitted ONCE, first ---------------------- //
void ej_emit_header(TMetapop* pop)
{
    std::ostringstream rec;
    rec << "{\"schema_version\":\"2\",\"kind\":\"header\""
        << ",\"replicate\":1"
        << ",\"patches\":" << pop->get_nbPatch();

    // loci: ntrl first then quanti; each declares its (declared) allele count,
    // which is the RAGGED inner length of that locus's frame allele_freqs.
    rec << ",\"loci\":[";
    {
        bool first = true;
        unsigned int nbProto = pop->getTraitPrototypeSize();
        for (int pass = 0; pass < 2; ++pass) {
            for (unsigned int t = 0; t < nbProto; ++t) {
                TTraitProto& proto = pop->getTraitPrototype(t);
                if (proto.get_type() != EJ_LOCUS_ORDER[pass]) continue;
                unsigned int nbLocus = proto.get_nb_locus();
                for (unsigned int l = 0; l < nbLocus; ++l) {
                    if (!first) rec << ",";
                    first = false;
                    rec << "{\"type\":\"" << EJ_LOCUS_ORDER[pass]
                        << "\",\"n_alleles\":" << proto.get_nb_allele(l) << "}";
                }
            }
        }
    }
    rec << "]";

    // traits: one per quantitative trait prototype (phenotype outer axis).
    rec << ",\"traits\":" << ej_count_quanti_traits(pop);

    // stats: the stream_key of every metapop-wide SCALAR estimator this run
    // computes, in the stable order ej_collect_scalar_stats() produces. Frames
    // carry a parallel `stats` array of values in this same order.
    std::vector<EjScalarStat> stats = ej_collect_scalar_stats(pop);
    rec << ",\"stats\":[";
    for (unsigned int s = 0; s < stats.size(); ++s) {
        if (s) rec << ",";
        rec << "\"" << stats[s].key << "\"";
    }
    rec << "]";

    rec << "}";
    std::cout << rec.str() << "\n" << std::flush;
}

// ---- frame: one logged generation, compact positional payload ------------ //
void ej_emit_frame(TMetapop* pop)
{
    const age_idx AGE = ADLTx;            // allele/phenotype payloads use ADULTS
    const sex_t   SEXES[2] = { FEM, MAL };

    unsigned int nbPatch = pop->get_nbPatch();
    unsigned int nbProto = pop->getTraitPrototypeSize();

    std::ostringstream rec;
    rec << "{\"schema_version\":\"2\",\"kind\":\"frame\",\"generation\":"
        << pop->getCurrentGeneration();

    // ---- allele_freqs[locus][patch][allele] (ntrl loci first, then quanti) ----
    rec << ",\"allele_freqs\":[";
    {
        bool firstLocus = true;
        for (int pass = 0; pass < 2; ++pass) {
            for (unsigned int t = 0; t < nbProto; ++t) {
                TTraitProto& proto = pop->getTraitPrototype(t);
                if (proto.get_type() != EJ_LOCUS_ORDER[pass]) continue;
                unsigned int absIdx  = proto.get_absolute_index();
                unsigned int nbLocus = proto.get_nb_locus();

                for (unsigned int l = 0; l < nbLocus; ++l) {
                    unsigned int nAll = proto.get_nb_allele(l);
                    if (!firstLocus) rec << ",";
                    firstLocus = false;

                    rec << "[";                      // patch axis
                    for (unsigned int p = 0; p < nbPatch; ++p) {
                        TPatch* patch = pop->get_vPatch(p);
                        std::vector<unsigned int> counts(nAll, 0);
                        unsigned int total = 0;

                        for (int s = 0; s < 2; ++s) {
                            sex_t SEX = SEXES[s];
                            unsigned int n = patch->size(SEX, AGE);
                            for (unsigned int i = 0; i < n; ++i) {
                                TIndividual* ind = patch->get(SEX, AGE, i);
                                if (!ind) continue;
                                ALLELE** seq = (ALLELE**)ind->getTrait((int)absIdx)->get_sequence();
                                if (!seq) continue;
                                for (unsigned int c = 0; c < ploidy; ++c) {
                                    ALLELE a = seq[l][c];
                                    if ((unsigned int)a < nAll) { counts[a]++; ++total; }
                                }
                            }
                        }

                        if (p) rec << ",";
                        rec << "[";                  // allele axis (length nAll)
                        for (unsigned int a = 0; a < nAll; ++a) {
                            if (a) rec << ",";
                            double f = (total > 0) ? (double)counts[a] / (double)total : 0.0;
                            ej_put_num(rec, f);
                        }
                        rec << "]";
                    }
                    rec << "]";
                }
            }
        }
    }
    rec << "]";

    // ---- phenotype[trait][patch][mean,var] (omitted when no quanti traits) ----
    {
        std::ostringstream ph;
        ph << ",\"phenotype\":[";
        bool anyTrait = false, firstTrait = true;
        for (unsigned int t = 0; t < nbProto; ++t) {
            TTraitProto& proto = pop->getTraitPrototype(t);
            if (proto.get_type() != "quanti") continue;
            anyTrait = true;
            unsigned int absIdx = proto.get_absolute_index();
            if (!firstTrait) ph << ",";
            firstTrait = false;

            ph << "[";                               // patch axis
            for (unsigned int p = 0; p < nbPatch; ++p) {
                TPatch* patch = pop->get_vPatch(p);
                double sum = 0.0, sum2 = 0.0;
                unsigned int n = 0;
                for (int s = 0; s < 2; ++s) {
                    sex_t SEX = SEXES[s];
                    unsigned int ns = patch->size(SEX, AGE);
                    for (unsigned int i = 0; i < ns; ++i) {
                        TIndividual* ind = patch->get(SEX, AGE, i);
                        if (!ind) continue;
                        double v = ind->getTrait((int)absIdx)->get_phenotype();
                        if (!std::isfinite(v) || v == (double)my_NAN) continue;
                        sum += v; sum2 += v * v; ++n;
                    }
                }
                // dense per-patch: an empty patch reports [0,0] (no individuals).
                double mean = n ? sum / n : 0.0;
                double var  = n ? sum2 / n - mean * mean : 0.0;
                if (var < 0.0) var = 0.0;            // guard tiny negative drift
                if (p) ph << ",";
                ph << "["; ej_put_num(ph, mean); ph << ","; ej_put_num(ph, var); ph << "]";
            }
            ph << "]";
        }
        ph << "]";
        if (anyTrait) rec << ph.str();
    }

    // ---- stats: bare numbers aligned positionally to header.stats ----
    // Recomputed from the SAME ordered recorder list the header declared, so the
    // i-th value here is the i-th stream_key there. Any scalar stat the run did
    // not request is simply absent from both (no nulls for missing stats).
    {
        std::vector<EjScalarStat> stats = ej_collect_scalar_stats(pop);
        if (!stats.empty()) {
            rec << ",\"stats\":[";
            for (unsigned int s = 0; s < stats.size(); ++s) {
                if (s) rec << ",";
                double v = stats[s].handler->computeStatNow(stats[s].rec);
                ej_put_num_or_null(rec, v);
            }
            rec << "]";
        }
    }

    rec << "}";
    std::cout << rec.str() << "\n" << std::flush;
}

} // namespace

// Driven once per generation by the FileServices (LCE_FileServicesNotifier, which
// runs AFTER the stat notifier, so the estimators are already up to date for this
// generation). We emit on the stat cadence and always on the final generation.
void EmitJsonFH::update()
{
    if (!_pop || !_statDb) return;
    unsigned int gen = _pop->getCurrentGeneration();
    unsigned int occ = _statDb->get_occurrence();
    if ((occ && !(gen % occ)) || gen == _pop->getGenerations()) FHwrite();
}

void EmitJsonFH::FHwrite()
{
    if (!_pop) return;
    // The header declares the run's static layout ONCE, before any frame.
    if (!_headerDone) {
        ej_emit_header(_pop);
        _headerDone = true;
    }
    ej_emit_frame(_pop);
}
