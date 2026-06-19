/** @file emit_json.cpp
 *
 *   JSON output mode for quantiNemo: CLI handling and the global mode switch.
 *   The per-generation JSON emitter lives in emit_json_record.inc, which is
 *   #included into a translation unit that already pulls in the StatHandler
 *   template definitions (lce_misc.cpp). See emit_json.h for the public API.
 *
 *   This file is part of quantiNemo and is distributed under the GNU General
 *   Public License v3 (or later), like the rest of the program.
 */

#include <cstring>

#include "emit_json.h"

// verbose_* are defined in main.cpp; we silence stdout chatter with them.
extern unsigned int verbose_message;
extern unsigned int verbose_warning;

bool g_emit_json = false;

// --------------------------------------------------------------------------- //
// CLI handling
// --------------------------------------------------------------------------- //
int emit_json_parse_cli(int argc, char** argv)
{
    int out = 0;
    char** dst = argv;
    for (int i = 0; i < argc; ++i) {
        const char* a = argv[i];
        if (i > 0 && strcmp(a, "--emit-json") == 0) {
            g_emit_json = true;       // recognised: drop the token
            continue;
        }
        *dst++ = argv[i];
        ++out;
    }
    if (g_emit_json) {
        // stdout must carry ONLY JSON records: silence message() AND warning()
        // (both print to stdout via vprintf/printf). error() already goes to
        // stderr and throws, so it never pollutes the JSON output.
        verbose_message = 0;
        verbose_warning = 0;
    }
    return out;
}
