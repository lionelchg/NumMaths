#ifndef _DEFS_
#define _DEFS_

// Length of strings
const int lenstr = 100;

// Explicit schemes (Euler time integration implied)
const int EXPL_SCHEME = 0;
const int LW_SCHEME = 0;
const int WB_SCHEME = 1;

// Semi-discretized schemes (more general for multi-step integrations)
const int SEMI_SCHEME = 1;
const int FOU_SCHEME = 0;

#endif