// -*- Mode: C; fill-column: 79 -*-

#ifndef __getpath_h
#define __getpath_h

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

void               setArgv(int argc, const char* const argv[]);
const char*        pathToExecutable(void);

int                getArgc(void);
const char* const* getArgv(void);
static void        pteDirname(char *filepath);
static void        pteJoinPath(char* const buffer, const char* const stuff);
bool               pteIsFile(const char* filepath);
bool               pteIsAnExecutableFile(const char *filepath);
bool               pteIsDir(const char* filepath);
void               pteCopyAbsolute(char* destStr, const char* filepath);
void               pteAbsolutize(char* filepath);


#ifdef __cplusplus
} // extern "C"
#endif

#endif // __getpath_h
