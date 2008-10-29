/*
 *
 *
 *
 */

#include <lpmd/util.h>
#include "plugincommon.h"
#include "config.h"
#include "version.h"

using namespace lpmd;

const char * PluginVersion()
{
 std::string lver;
 lver = VERSION;
 #ifndef NUMBERED_RELEASE
 lver += " (from ";
 lver += SVNBRANCH;
 lver += (", revision "+ToString<int>(SVNREVISION)+")");
 #endif
 return lver.c_str();
}

