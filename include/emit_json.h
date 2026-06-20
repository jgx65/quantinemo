/** @file emit_json.h
 *
 *   JSON output mode for quantiNemo.
 *
 *   ADDITIVE, opt-in. When the engine is invoked with --emit-json it:
 *     - requires exactly ONE simulation and ONE replicate;
 *     - silences the normal console chatter so stdout carries ONLY JSON;
 *     - emits a header + frames stream of line-delimited JSON.
 *
 *   The stream is ONE `kind:"header"` record first — declaring the run's static
 *   layout once (patch count, the locus table with per-locus allele counts, the
 *   trait count, and which layer-1 estimators are computed) — followed by ONE
 *   compact `kind:"frame"` record per logged generation (every stat_log_time,
 *   plus the final generation) carrying only positional numeric payloads
 *   (allele_freqs, phenotype, stats) read against the header.
 *
 *   This file is part of quantiNemo and is distributed under the GNU General
 *   Public License v3 (or later), like the rest of the program.
 */

#ifndef emit_jsonH
#define emit_jsonH

class TMetapop;

/** Global switch: true when --emit-json was passed. Defined in emit_json.cpp,
 *  defaulted to false so absent-flag behaviour is byte-for-byte unchanged. */
extern bool g_emit_json;

/** Scan argv for --emit-json; if present, set g_emit_json, silence the normal
 *  stdout chatter (verbose_message=0) and remove the flag token from argv so the
 *  existing .ini/CLI parser never sees it. Returns the new argc. Safe to call
 *  once from main() before constructing TSimManager. */
int emit_json_parse_cli(int argc, char** argv);

/** Emit one JSON record (single line, flushed) for the current generation of
 *  \a pop to stdout. No-op when g_emit_json is false. Called from the
 *  stat-notifier hook (LCE_StatServiceNotifier::execute) at the stat_log_time
 *  cadence and on the final generation. */
void emit_json_record(TMetapop* pop);

#endif
