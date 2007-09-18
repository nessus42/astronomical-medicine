// -*- Mode: C; fill-column: 79 -*-

#ifndef __getpath_h
#define __getpath_h

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

void               setArgv(int argc, const char* const argv[]);
const char*        pathToExecutable(void);

const char*        pathToExecutableDir(void);
int                getArgc(void);
const char* const* getArgv(void);
void               pteDirname(char *filepath);
void               pteJoinPathM(char* buffer, const char* stuff);
const char*        pteJoinPath(const char* buffer, const char* stuff);
bool               pteIsFile(const char* filepath);
bool               pteIsAnExecutableFile(const char *filepath);
bool               pteIsDir(const char* filepath);
void               pteCopyAbsolute(char* destStr, const char* filepath);
void               pteAbsolutize(char* filepath);


#ifdef __cplusplus
} // extern "C"
#endif

#endif // __getpath_h
