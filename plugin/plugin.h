#ifndef H_PLUGIN
#define H_PLUGIN

#include <stdio.h>
#include <stdlib.h>

#include <molfile_plugin.h>
#include <Gromacs.h>

extern "C"
static void *open_trr_read(const char *filename, const char *filetype, int *natom);
extern "C"
static int read_trr_timestep(void *v, int natoms, molfile_timestep_t *ts);
extern "C"
static void close_trr_read(void *v);
extern "C"
static void *open_trr_write(const char *filename, const char *filetype,int natoms);
extern "C"
static int write_trr_timestep(void *mydata, const molfile_timestep_t *ts);
extern "C"
static void close_trr_write(void *v);
#endif
