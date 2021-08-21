#include <stdio.h>
#include <stdlib.h>
#include "argp.h"

const char *argp_program_version =
    "argex 1.0";

const char *argp_program_bug_address =
    "<bug-gnu-utils@gnu.org>";

/* This structure is used by main to communicate with parse_opt. */
struct arguments
{
    // char **args;           /* ARG1 and ARG2 */
    char *args[10];           /* ARG1 and ARG2 */
    int nschemes;             /* for Number of schemes */
    int verbose;             /* The -v flag */
    char *outfile;           /* Argument for -o */
};

/*
   OPTIONS.  Field 1 in ARGP.
   Order of fields: {NAME, KEY, ARG, FLAGS, DOC}.
*/
static struct argp_option options[] =
    {
        {"verbose", 'v', 0, 0, "Produce verbose output"},
        {"nschemes", 'n', "NSCHEMES", 0, "Number of schemes"},
        {"output", 'o', "OUTFILE", 0,
         "Output to OUTFILE instead of to standard output"},
        {0}};

/*
   PARSER. Field 2 in ARGP.
   Order of parameters: KEY, ARG, STATE.
*/
static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
    struct arguments *arguments = state->input;

    switch (key)
    {
    case 'v':
        arguments->verbose = 1;
        break;
    case 'o':
        arguments->outfile = arg;
        break;
    case 'n':
        arguments->nschemes = atoi(arg);
    case ARGP_KEY_ARG:
        arguments->args[state->arg_num] = arg;
        break;
    default:
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

/*
   ARGS_DOC. Field 3 in ARGP.
   A description of the non-option command-line arguments
     that we accept.
*/
static char args_doc[] = "ARG1 ARG2";

/*
  DOC.  Field 4 in ARGP.
  Program documentation.
*/
static char doc[] =
    "argex -- A program to demonstrate how to code command-line options \
        and arguments.\vFrom the GNU C Tutorial.";

/*
   The ARGP structure itself.
*/
static struct argp argp = {options, parse_opt, args_doc, doc};

