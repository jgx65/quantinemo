/** @file emit_json_record.inc
 *
 *   JSON record emitter for --emit-json mode. This is #included at the END of
 *   lce_misc.cpp (NOT compiled standalone) because it calls the StatHandler<SH>
 *   template methods (getFst_WC, getHsnei, getQst, ...) whose bodies live in
 *   stathandler.cpp, which lce_misc.cpp already pulls in. That guarantees the
 *   concrete TTNeutralSH / TTQuantiSH instantiations are visible.
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
 *   See emit_json.h for the public API. Distributed under GPL v3 (or later),
 *   like the rest of quantiNemo.
 */

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "emit_json.h"
#include "patch.h"
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
    tmp.precision(10);
    tmp << v;
    os << tmp.str();
}

void ej_put_num(std::ostream& os, double v)
{
    std::ostringstream tmp;
    tmp.precision(10);
    tmp << v;
    os << tmp.str();
}

// Layout decided once when the header is emitted, then reused by every frame so
// the positional `stats` payload stays aligned to the header's `stats` list.
bool ej_header_done = false;
bool ej_have_fst = false, ej_have_qst = false,
     ej_have_ho  = false, ej_have_hs  = false, ej_have_ht = false;

// canonical locus / genotype ordering: all neutral loci first, then quanti
const char* const EJ_LOCUS_ORDER[2] = { "ntrl", "quanti" };

unsigned int ej_count_quanti_traits(TMetapop* pop)
{
    unsigned int n = 0, nbProto = pop->getTraitPrototypeSize();
    for (unsigned int t = 0; t < nbProto; ++t)
        if (pop->getTraitPrototype(t).get_type() == "quanti") ++n;
    return n;
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

    // stats: which layer-1 estimators this run computes, in canonical order.
    // A neutral stat handler provides fst/ho/hs/ht; a quanti one provides qst.
    TTraitNeutralProto* ntrlProto =
        dynamic_cast<TTraitNeutralProto*>(pop->getFirstPrototype("ntrl"));
    TTraitQuantiProto* quantiProto =
        dynamic_cast<TTraitQuantiProto*>(pop->getFirstPrototype("quanti"));
    ej_have_fst = ej_have_ho = ej_have_hs = ej_have_ht =
        (ntrlProto && ntrlProto->_stats);
    ej_have_qst = (quantiProto && quantiProto->_stats);

    rec << ",\"stats\":[";
    {
        bool first = true;
        if (ej_have_fst) {                rec << "\"fst\""; first = false; }
        if (ej_have_qst) { if (!first) rec << ","; rec << "\"qst\""; first = false; }
        if (ej_have_ho)  { if (!first) rec << ","; rec << "\"ho\"";  first = false; }
        if (ej_have_hs)  { if (!first) rec << ","; rec << "\"hs\"";  first = false; }
        if (ej_have_ht)  { if (!first) rec << ","; rec << "\"ht\"";  first = false; }
    }
    rec << "]";

    rec << "}";
    std::cout << rec.str() << "\n" << std::flush;
}

// ---- frame: one logged generation, compact positional payload ------------ //
void ej_emit_frame(TMetapop* pop)
{
    const age_idx AGE = ADLTx;            // stats default to ADULTS
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
    if (ej_have_fst || ej_have_qst || ej_have_ho || ej_have_hs || ej_have_ht) {
        TTraitNeutralProto* ntrlProto =
            dynamic_cast<TTraitNeutralProto*>(pop->getFirstPrototype("ntrl"));
        TTraitQuantiProto* quantiProto =
            dynamic_cast<TTraitQuantiProto*>(pop->getFirstPrototype("quanti"));

        double fst = (double)my_NAN, ho = (double)my_NAN,
               hs = (double)my_NAN, ht = (double)my_NAN, qst = (double)my_NAN;
        if (ntrlProto && ntrlProto->_stats) {
            TTNeutralSH* sh = ntrlProto->_stats;
            fst = sh->getFst_WC(ADULTS);    // Weir & Cockerham theta
            ho  = sh->getHo(ADULTS);        // observed heterozygosity
            hs  = sh->getHsnei(ADULTS);     // Nei within-patch expected het
            ht  = sh->getHtnei(ADULTS);     // Nei total expected het
        }
        if (quantiProto && quantiProto->_stats)
            qst = quantiProto->_stats->getQst(ADULTS);

        rec << ",\"stats\":[";
        bool first = true;
        if (ej_have_fst) {                ej_put_num_or_null(rec, fst); first = false; }
        if (ej_have_qst) { if (!first) rec << ","; ej_put_num_or_null(rec, qst); first = false; }
        if (ej_have_ho)  { if (!first) rec << ","; ej_put_num_or_null(rec, ho);  first = false; }
        if (ej_have_hs)  { if (!first) rec << ","; ej_put_num_or_null(rec, hs);  first = false; }
        if (ej_have_ht)  { if (!first) rec << ","; ej_put_num_or_null(rec, ht);  first = false; }
        rec << "]";
    }

    rec << "}";
    std::cout << rec.str() << "\n" << std::flush;
}

} // namespace

void emit_json_record(TMetapop* pop)
{
    if (!g_emit_json || !pop) return;

    // The header declares the run's static layout ONCE, before any frame.
    if (!ej_header_done) { ej_emit_header(pop); ej_header_done = true; }
    ej_emit_frame(pop);
}
