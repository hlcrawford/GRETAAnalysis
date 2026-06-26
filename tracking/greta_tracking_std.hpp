/**
 * Prepended to all vendored tracking translation units (see SConstruct).
 * Replaces the old C-only tracking_compat.h forced-include.
 */
#pragma once

#include <cassert>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#if !defined(_WIN32)
#include <strings.h> /* bzero (POSIX); omit on Windows */
#endif

#ifdef GRETA_TRACKING_USE_WORKSPACE
#include "GretaTrackingWorkspace.h"

/* Pair-production scratch (per active workspace). */
#define GRETA_WS (GretaTrackingActiveWorkspace())
#define npp (GRETA_WS->nPairProd)
#define ppe (GRETA_WS->ppE)
#define ppx (GRETA_WS->ppX)
#define ppy (GRETA_WS->ppY)
#define ppz (GRETA_WS->ppZ)
#define ppfom (GRETA_WS->ppFom)
#define ppTS (GRETA_WS->ppTs)
#endif

#ifdef GRETA_TRACK_DEBUG
#include "GretaTrackingDebug.h"
#endif
