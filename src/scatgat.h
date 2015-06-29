/*
 * scatgat.h
 * Jun 25, 2015 20:54:32 EDT
 * Copyright 2015
 *        
 * Andre Young <andre.young@cfa.harvard.edu>
 * Harvard-Smithsonian Center for Astrophysics
 * 60 Garden Street, Cambridge
 * MA 02138
 * 
 * Changelog:
 * 	AY: Created 2015-06-25
 */

#ifndef SCATGAT_H
#define SCATGAT_H

#define _GNU_SOURCE 

#include <fcntl.h>
#include <limits.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>

#include "sg_access.h"
#include "dplane_proxy.h"

#define LAST_VDIF_SECS_INRE(a) ((VDIFHeader *)(&(a->data_buf[(a->n_frames-1)*(a->sgi->pkt_size)/sizeof(uint32_t)])))->w1.secs_inre
#define FIRST_VDIF_SECS_INRE(a) ((VDIFHeader *)(a->data_buf))->w1.secs_inre
#define LAST_VDIF_DF_NUM_INSEC(a) ((VDIFHeader *)(&(a->data_buf[(a->n_frames-1)*(a->sgi->pkt_size)/sizeof(uint32_t)])))->w2.df_num_insec
#define FIRST_VDIF_DF_NUM_INSEC(a) ((VDIFHeader *)(a->data_buf))->w2.df_num_insec

/* Set SGPlan to read / write mode */
enum scatgat_mode {
	SCATGAT_MODE_READ,
	SCATGAT_MODE_WRITE
};

/* Encapsulates single SG file */
typedef struct sg_part {
	SGInfo *sgi;														// points to SGInfo for single SG file
	off_t iblock; 														// next block to read from / write to in SG file
	uint32_t *data_buf; 												// points to start VDIF buffer from previous read / for pending write
	uint32_t n_frames; 													// number of VDIF frames in buffer
} SGPart;

/* Encapsulates group of SG files */
typedef struct sg_plan {
	int sgm;															// scatgat_mode: read / write
	int n_sgprt; 														// number of SGPart elements
	SGPart *sgprt; 														// array of SGPart elements (one per SG file)
} SGPlan;

/*
 * Allocate memory and fill it with SGInfo instances.
 * Arguments:
 *   SGPlan **sgplan -- Address of SGPlan pointer to allocate memory.
 *   const char *pattern -- Filename pattern to search for.
 *   const char *fmtstr -- Format string used to compile the file full
 *     path. It should have the form <..>%d<..>%d<..>%s where the first
 *     %d is replaced with an element from mod_list, the second %d 
 *     replaced with an element from disk_list, and the %s replaced with
 *     pattern.
 *   int *mod_list -- Array of module numbers to use.
 *   int n_mod -- Number of modules to use.
 *   int *disk_list -- Array of disk numbers to use.
 *   int n_disk -- Number of disks to use.
 * Returns:
 *   int -- Number of SGInfo instances (SG files found mathcing pattern)
 * Notes:
 *   All names that match the pattern /mnt/disks/MOD/DISK/data/PATTERN
 *     where MOD and DISK are elements of mod_list and disk_list, 
 *     respectively, and PATTERN is the string in pattern, are given
 *     to sg_access. For each valid SG file found an SGInfo entry is 
 *     allocated in the buffer pointed to by *sgi.
 *   The SGInfo entries stored in sgplan are sorted in ascending order
 *     according to the timestamp on the first VDIF frame in each SG 
 *     file.
 *   For each valid SG file encountered an SGPart element is stored in
 *     SGPlan.
 */
int make_sg_read_plan(SGPlan **sgpln, const char *pattern, 
					const char *fmtstr, int *mod_list, int n_mod, 
					int *disk_list, int n_disk);

/*
 * Read the next block of VDIF frames.
 * Arguments:
 *   SGPlan *sgpln -- The SGPlan created for a given filename pattern.
 *   uint32_t **data_buf -- Address of pointer which can be used to 
 *     store the location of the data buffer created and filled by 
 *     reading the next block.
 * Returns:
 *   int -- The number of VDIF frames contained in the buffer, zero if
 *     end of all files reached, and -1 if the data is no longer 
 *     contiguous.
 * Notes:
 *   This method attempts to read a contiguous set of VDIF frames that
 *     is the equivalent of one SG block per SG file contained in the SG
 *     plan. Blocks from different files are stitched together such that
 *     the first frame in one block directly follows the last frame of
 *     another block. Blocks with data that do not flow contiguously 
 *     from the first frame for the current block are stored in the 
 *     buffer of the associated SGPart, inside SGPlan. Upon subsequent 
 *     calls to this method no further blocks of data is read from that 
 *     particular SG file until its block can be stitched togther with 
 *     the contiguous flow.
 *   Block counter for each SGPart is updated if frames where read from
 *     that file.
 */
int read_next_block_vdif_frames(SGPlan *sgpln, uint32_t **vdif_buf);

/*
 * Read one block's worth of VDIF frames from a group of SG files.
 * Arguments:
 *   SGInfo *sgi -- Array of valid SGInfo instances (i.e. had to have
 *     been accessed prior to calling this method).
 *   int n_sgi -- Number of SGInfo elements in the array.
 *   off_t iblock -- The block index.
 *   uint32_t **vdif_buf -- Address of a pointer which can be used to
 *     store the location of the VDIF buffer created and filled.
 * Returns:
 *   int -- The number of VDIF frames contained in the buffer.
 * Notes:
 *   This method creates as many threads as there are SGInfo elements in 
 *     the array, using pthread_create with sgthread_read_block as the
 *     thread start method.
 *   The VDIF buffer size is determined by counting the packets per 
 *     block total for all SGInfo instances, although the actual used
 *     size may be smaller if some of the blocks are short.
 */
int read_block_vdif_frames(SGPlan *sgpln, off_t iblock, 
							uint32_t **vdif_buf);

/*
 * Close scatter gather reader plan
 */
void close_sg_read_plan(SGPlan *sgplan);

/*
 * Make scatter gather write plan. 
 */
int make_sg_write_plan(SGPlan **sgpln, const char *pattern, 
					const char *fmtstr, int *mod_list, int n_mod, 
					int *disk_list, int n_disk);

/*
 * Write VDIF buffer.
 */
int write_vdif_frames(SGPlan *sgpln, uint32_t *vdif_buf, int n_frames);

/*
 * Close scatter gather write plan
 */
void close_sg_write_plan(SGPlan *sgpln);

/*
 * Free the resources allocated for an SGPlan structure.
 * Arguments:
 *   SGPlan *sgpln -- Pointer to allocated SGPlan structure.
 */
void free_sg_plan(SGPlan *sgpln);

#endif // SCATGAT_H
