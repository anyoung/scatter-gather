/*
 * scatgat.c
 * Jun 25, 2015 21:14:25 EDT
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

#include "scatgat.h"

/* Measures the number of VDIF frames per second and determines where
 * the VDIFHeader.df_num_insec wraps to zero.
 */
#define VDIF_FRAMES_PER_SECOND 125000

/* File permissions with which scatter-gather files are created. */ 
#define SG_FILE_PERMISSIONS (S_IWUSR | S_IRUSR | S_IWGRP | S_IRGRP | S_IROTH)
#define SG_FILE_WRITE_OPEN_MODE (O_WRONLY|O_APPEND)
#define SG_MMAP_WRITE_OPEN_PROTO (PROT_WRITE)
#define SG_MMAP_WRITE_OPEN_MODE (MAP_SHARED)

/* Debugging utilities */
#define DEBUG_LEVEL_DEBUG 40
#define DEBUG_LEVEL_INFO 30
#define DEBUG_LEVEL_WARNING 20
#define DEBUG_LEVEL_ERROR 10
//~ #define DEBUG_LEVEL DEBUG_LEVEL_DEBUG
//~ #define DEBUG_LEVEL DEBUG_LEVEL_INFO
//~ #define DEBUG_LEVEL DEBUG_LEVEL_WARNING
//~ #define DEBUG_LEVEL DEBUG_LEVEL_ERROR
#ifdef DEBUG_LEVEL
#define _DBGMSGLEN 0x400
void debug_msg(const char *msg, const char *filename, const char *funcname, int linenum);
#define DEBUGMSG(m) debug_msg(m,__FILE__,__FUNCTION__,__LINE__)
void error_msg(const char *msg, const char *filename, const char *funcname, int linenum);
#define ERRORMSG(m) error_msg(m,__FILE__,__FUNCTION__,__LINE__)
void warning_msg(const char *msg, const char *filename, const char *funcname, int linenum);
#define WARNINGMSG(m) warning_msg(m,__FILE__,__FUNCTION__,__LINE__)
void info_msg(const char *msg, const char *filename, const char *funcname, int linenum);
#define INFOMSG(m) info_msg(m,__FILE__,__FUNCTION__,__LINE__)
void print_sg_part(SGPart *sgprt, const char *label);
void print_sg_plan(SGPlan *sgpln, const char *label);
#define DEBUGMSG_ENTERFUNC snprintf(_dbgmsg,_DBGMSGLEN,"Enter %s.",__FUNCTION__); DEBUGMSG(_dbgmsg)
#define DEBUGMSG_LEAVEFUNC snprintf(_dbgmsg,_DBGMSGLEN,"Leave %s.",__FUNCTION__); DEBUGMSG(_dbgmsg)
#endif

/* Comparison methods, may be used with qsort */
int compare_int_descend(const void *a, const void *b);
int compare_sg_info(const void *a, const void *b);
int compare_sg_part(const void *a, const void *b);

/* For sorting and continuity testing */
int map_sg_parts_contiguous(SGPlan *sgpln, int *mapping);
int test_sg_parts_contiguous(SGPart *a, SGPart *b);

/* Threaded implementations compatible with pthread */
static void * sgthread_read_block(void *arg);
static void * sgthread_fill_read_sgi(void *arg);

/* Memory management */
void clear_sg_part_buffer(SGPart *sgprt);
void free_sg_info(SGInfo *sgi);

/*
 * Description needed.
 */
int compare_int_descend(const void *a, const void *b)
{
	int *int_a = (int *)a;
	int *int_b = (int *)b;
	return *int_b < *int_a ? -1 : *int_b > *int_a;
}

//////////////////////////////////////////////////////////////////////// SCATTER GATHER READING
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
						int *disk_list, int n_disk)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	int itmp; // just a counter
	int idisk, imod; // disk, module counters
	char filename[n_mod*n_disk][PATH_MAX]; // full filename searched for
	int ithread; // thread counter
	int thread_result; // return result for pthread methods
	pthread_t sg_threads[n_mod*n_disk]; // pthreads to do filling
	int valid_sgi = 0; // number of valid SG files found
	/* Allocate temporary buffer to store maximum possible SGInfo 
	 * instances.
	 */
	SGInfo *sgi_buf = (SGInfo *)calloc(sizeof(SGInfo),n_mod*n_disk);
	/* And allocate temporary single SGInfo. */
	SGInfo *sgi_tmp = (SGInfo *)calloc(sizeof(SGInfo),1); 
	/* Step through all modules and disks, and access files that 
	 * match the pattern.
	 */
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG("\tLaunching threads.");
	#endif
	for (imod=0; imod<n_mod; imod++)
	{
		#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
			snprintf(_dbgmsg,_DBGMSGLEN,"\t\tmod[%d] = %d",imod,mod_list[imod]);
			DEBUGMSG(_dbgmsg);
		#endif
		for (idisk=0; idisk<n_disk; idisk++)
		{
			#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
				snprintf(_dbgmsg,_DBGMSGLEN,"\t\t\tdisk[%d] = %d",idisk,disk_list[idisk]);
				DEBUGMSG(_dbgmsg);
			#endif
			ithread = imod*n_disk + idisk;
			snprintf(filename[ithread],PATH_MAX,fmtstr,mod_list[imod],disk_list[idisk],pattern);
			#if DEBUG_LEVEL >= DEBUG_LEVEL_INFO
				snprintf(_dbgmsg,_DBGMSGLEN,"\t\t\tAccessing file '%s'.",filename[ithread]);
				INFOMSG(_dbgmsg);
			#endif
			thread_result = pthread_create(&(sg_threads[ithread]), NULL, &sgthread_fill_read_sgi, filename[ithread]);
			if (thread_result != 0)
			{
				perror("Unable to launch thread.");
				exit(EXIT_FAILURE);
			}
		}
	}
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG("\tJoining threads.");
	#endif
	for (imod=0; imod<n_mod; imod++)
	{
		#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
			snprintf(_dbgmsg,_DBGMSGLEN,"\t\tmod[%d] = %d",imod,mod_list[imod]);
			DEBUGMSG(_dbgmsg);
		#endif
		for (idisk=0; idisk<n_disk; idisk++)
		{
			#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
				snprintf(_dbgmsg,_DBGMSGLEN,"\t\t\tdisk[%d] = %d",idisk,disk_list[idisk]);
				DEBUGMSG(_dbgmsg);
			#endif
			ithread = imod*n_disk + idisk;
			thread_result = pthread_join(sg_threads[ithread],(void *)&sgi_tmp);
			if (thread_result != 0)
			{
				perror("Unable to join thread.");
				exit(EXIT_FAILURE);
			}
			if (sgi_tmp->smi.mmfd > 0)
			{
				memcpy(sgi_buf+valid_sgi++, sgi_tmp, sizeof(SGInfo));
				#if DEBUG_LEVEL >= DEBUG_LEVEL_INFO
					sg_report(&(sgi_buf[valid_sgi-1]),"\tSG report (sgi_buf):");
					DEBUGMSG("\t\t\tClosing SGInfo.");
				#endif
				//~ sg_close(sgi_tmp);
			}
			/* Free the temporary SGInfo resources, but DO NOT free
			 * the malloc'ed NAME to which we still keep a pointer.
			 */
			free(sgi_tmp);
		}
	}
	if (valid_sgi == 0)
	{
		/* Done with this, free it. */
		free(sgi_buf);
		#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
			DEBUGMSG_LEAVEFUNC;
		#endif
		return 0;
	}
	/* Allocate space for storing valid SGInfos. */
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG("\tAllocating buffer space and copying SGInfo.");
	#endif
	
	/* Sort SGInfo array according to second / frames */
	qsort((void *)sgi_buf, valid_sgi, sizeof(SGInfo), compare_sg_info);
	/* Allocate memory for SGPlan */
	*sgpln = (SGPlan *)malloc(sizeof(SGPlan));
	(*sgpln)->sgprt = (SGPart *)malloc(sizeof(SGPart)*valid_sgi);
	for (itmp=0; itmp<valid_sgi; itmp++)
	{
		// Allocate memory for SGInfo and copy
		(*sgpln)->sgprt[itmp].sgi = (SGInfo *)malloc(sizeof(SGInfo));
		memcpy((*sgpln)->sgprt[itmp].sgi, &(sgi_buf[itmp]), sizeof(SGInfo));
		// Initialize block counter, VDIF buffer, and frame counter
		(*sgpln)->sgprt[itmp].iblock = 0;
		(*sgpln)->sgprt[itmp].data_buf = NULL;
		(*sgpln)->sgprt[itmp].n_frames = 0;
	}
	(*sgpln)->n_sgprt = valid_sgi;
	/* Done with the temporary buffer, free it. */
	free(sgi_buf);
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		print_sg_plan(*sgpln,"\t");
		DEBUGMSG_LEAVEFUNC;
	#endif
	return valid_sgi;
}

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
int read_next_block_vdif_frames(SGPlan *sgpln, uint32_t **vdif_buf)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
		print_sg_plan(sgpln,"\t");
	#endif
	int ithread; // thread counter
	int thread_result; // result of calls to pthread methods
	pthread_t sg_threads[sgpln->n_sgprt]; // the pthreads used
	int sg_threads_mask[sgpln->n_sgprt];
	
	int frames_estimate = 0; // estimate the size of buffer to create
	int frames_read = 0; // count the number of frames received
	int frame_size = sgpln->sgprt[0].sgi->pkt_size; // size of a frame
	
	int isgprt;
	int n_contiguous_blocks = 0;
	int mapping[sgpln->n_sgprt];
	
	/* Launch threads to read data */
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG("\tLaunching threads.");
	#endif
	for (ithread=0; ithread<sgpln->n_sgprt; ithread++)
	{
		/* For each SGPart, check if its data buffer is empty, which 
		 * indicates that the next block of data should be read.
		 */
		sg_threads_mask[ithread] = 0;
		if (sgpln->sgprt[ithread].n_frames == 0 && sgpln->sgprt[ithread].iblock < sgpln->sgprt[ithread].sgi->sg_total_blks)
		{
			sg_threads_mask[ithread] = 1;
			thread_result = pthread_create(&(sg_threads[ithread]),NULL,&sgthread_read_block,&(sgpln->sgprt[ithread]));
			if (thread_result != 0)
			{
				perror("Unable to create thread.");
				exit(EXIT_FAILURE);
			}
			frames_estimate += sgpln->sgprt[ithread].sgi->sg_wr_pkts;
		}
	}
	/* Create storage buffer. Assume that the number of frames read
	 * is always smaller than or equal to the number of estimated frames
	 */
	*vdif_buf = (uint32_t *)malloc(frames_estimate*frame_size);
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG("\tJoining threads.");
	#endif
	/* Join the threads */
	for (ithread=0; ithread<sgpln->n_sgprt; ithread++)
	{
		/* Only join threads that have been started. */
		if (sg_threads_mask[ithread] == 1)
		{
			thread_result = pthread_join(sg_threads[ithread],NULL);
			//printf("Thread %d read %d frames.\n",ithread,msg_in->num_frames);
			if (thread_result != 0)
			{
				perror("Unable to join thread.");
				exit(EXIT_FAILURE);
			}
			/* If we read frames from this SG file, update the block 
			 * counter.
			 */
			if (sgpln->sgprt[ithread].n_frames > 0)
			{
				sgpln->sgprt[ithread].iblock++;
			}
		}
	}
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		print_sg_plan(sgpln,"\t");
	#endif
	
/***********************************************************************
 * For now, ignore discontinuity in the data.
 *
 */
	n_contiguous_blocks = map_sg_parts_contiguous(sgpln, mapping);
	if (n_contiguous_blocks == 0)
	{
		printf("No contiguous blocks found.\n");
		return 0;
	}
	for (isgprt=0; isgprt<n_contiguous_blocks; isgprt++)
	{
		memcpy((void *)(*vdif_buf + frames_read*frame_size/sizeof(uint32_t)),
				(void *)(sgpln->sgprt[mapping[isgprt]-1].data_buf),sgpln->sgprt[mapping[isgprt]-1].n_frames*frame_size);
		frames_read += sgpln->sgprt[mapping[isgprt]-1].n_frames;
		clear_sg_part_buffer(&(sgpln->sgprt[mapping[isgprt]-1]));
	}
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		snprintf(_dbgmsg,_DBGMSGLEN,"Found %d contiguous blocks\n",n_contiguous_blocks);
		DEBUGMSG(_dbgmsg);
	#endif
 /*
 *
 **********************************************************************/
	//~ for (ithread=0; ithread<sgpln->n_sgprt; ithread++)
	//~ {
		//~ if (sg_threads_mask[ithread] == 1 && sgpln->sgprt[ithread].n_frames > 0)
		//~ {
			//~ memcpy((void *)(*vdif_buf + frames_read*frame_size/sizeof(uint32_t)),
					//~ (void *)(sgpln->sgprt[ithread].data_buf),sgpln->sgprt[ithread].n_frames*frame_size);
			//~ frames_read += sgpln->sgprt[ithread].n_frames;
			//~ clear_sg_part_buffer(&(sgpln->sgprt[ithread]));
		//~ }
	//~ } 
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		print_sg_plan(sgpln,"\t");
		DEBUGMSG_LEAVEFUNC;
	#endif
	/* Return the frame count. */
	return frames_read;
}

// Reading VDIF from SG files
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
 *     size may be smaller if one of the blocks is short.
 */
int read_block_vdif_frames(SGPlan *sgpln, off_t iblock, uint32_t **vdif_buf)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	int ithread; // thread counter
	int thread_result; // result of calls to pthread methods
	pthread_t sg_threads[sgpln->n_sgprt]; // the pthreads used
	
	int frames_estimate = 0; // estimate the size of buffer to create
	int frames_read = 0; // count the number of frames received
	int frame_size = sgpln->sgprt[0].sgi->pkt_size; // size of a frame
	
	/* Launch threads to read data */
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG("\tLaunching threads.");
	#endif
	for (ithread=0; ithread<sgpln->n_sgprt; ithread++)
	{
		thread_result = pthread_create(&(sg_threads[ithread]),NULL,&sgthread_read_block,&(sgpln->sgprt[ithread]));
		if (thread_result != 0)
		{
			perror("Unable to create thread.");
			exit(EXIT_FAILURE);
		}
		frames_estimate += sgpln->sgprt[ithread].sgi->sg_wr_pkts;
	}
	/* Create storage buffer. */
	*vdif_buf = (uint32_t *)malloc(frames_estimate*frame_size);
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG("\tJoining threads.");
	#endif
	/* Join the threads and copy data. */
	for (ithread=0; ithread<sgpln->n_sgprt; ithread++)
	{
		thread_result = pthread_join(sg_threads[ithread],NULL);
		if (thread_result != 0)
		{
			perror("Unable to join thread.");
			exit(EXIT_FAILURE);
		}
		if (sgpln->sgprt[ithread].n_frames > 0)
		{
			memcpy((void *)(*vdif_buf + frames_read*frame_size/sizeof(uint32_t)),(void *)(sgpln->sgprt[ithread].data_buf),sgpln->sgprt[ithread].n_frames*frame_size);
			frames_read += sgpln->sgprt[ithread].n_frames;
		}
	}
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
	/* Return the frame count. */
	return frames_read;
}

//////////////////////////////////////////////////////////////////////// SCATTER GATHER WRITING
/*
 * Make scatter gather write plan. 
 */
int make_sg_write_plan(SGPlan **sgpln, const char *pattern, 
					const char *fmtstr, int *mod_list, int n_mod, 
					int *disk_list, int n_disk, int pkt_size)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	int idisk, imod; // disk, module counters
	char filename[n_mod*n_disk][PATH_MAX]; // full filename searched for
	int ithread; // thread counter
	int thread_result; // return result for pthread methods
	pthread_t sg_threads[n_mod*n_disk]; // pthreads to do filling
	int valid_sgi = 0; // number of valid SG files created
	/* Step through all modules and disks, and access files that 
	 * match the pattern.
	 */
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG("\tLaunching threads.");
	#endif
	for (imod=0; imod<n_mod; imod++)
	{
		#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
			snprintf(_dbgmsg,_DBGMSGLEN,"\t\tmod[%d] = %d",imod,mod_list[imod]);
			DEBUGMSG(_dbgmsg);
		#endif
		for (idisk=0; idisk<n_disk; idisk++)
		{
			#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
				snprintf(_dbgmsg,_DBGMSGLEN,"\t\t\tdisk[%d] = %d",idisk,disk_list[idisk]);
				DEBUGMSG(_dbgmsg);
			#endif
			ithread = imod*n_disk + idisk;
			snprintf(filename[ithread],PATH_MAX,fmtstr,mod_list[imod],disk_list[idisk],pattern);
			#if DEBUG_LEVEL >= DEBUG_LEVEL_INFO
				snprintf(_dbgmsg,_DBGMSGLEN,"\t\t\tCreating file '%s'.",filename[ithread]);
				INFOMSG(_dbgmsg);
			#endif
			//~ thread_result = pthread_create(&(sg_threads[ithread]), NULL, &sgthread_fill_write_sgi, filename[ithread]);
			//~ if (thread_result != 0)
			//~ {
				//~ perror("Unable to launch thread.");
				//~ exit(EXIT_FAILURE);
			//~ }
		}
	}
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG("\tLaunching threads.");
	#endif
	for (imod=0; imod<n_mod; imod++)
	{
		#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
			snprintf(_dbgmsg,_DBGMSGLEN,"\t\tmod[%d] = %d",imod,mod_list[imod]);
			DEBUGMSG(_dbgmsg);
		#endif
		for (idisk=0; idisk<n_disk; idisk++)
		{
			#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
				snprintf(_dbgmsg,_DBGMSGLEN,"\t\t\tdisk[%d] = %d",idisk,disk_list[idisk]);
				DEBUGMSG(_dbgmsg);
			#endif
			ithread = imod*n_disk + idisk;
			snprintf(filename[ithread],PATH_MAX,fmtstr,mod_list[imod],disk_list[idisk],pattern);
			#if DEBUG_LEVEL >= DEBUG_LEVEL_INFO
				snprintf(_dbgmsg,_DBGMSGLEN,"\t\t\tCreating file '%s'.",filename[ithread]);
				INFOMSG(_dbgmsg);
			#endif
			//~ thread_result = pthread_join(&(sg_threads[ithread]), NULL);
			//~ if (thread_result != 0)
			//~ {
				//~ perror("Unable to launch thread.");
				//~ exit(EXIT_FAILURE);
			//~ }
		}
	}
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
}

/*
 * Write VDIF buffer.
 */
int write_next_block_vdif_frames(SGPlan *sgpln, uint32_t *vdif_buf, 
									int n_frames)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	int ithread; // thread counter
	int thread_result; // result of calls to pthread methods
	pthread_t sg_threads[sgpln->n_sgprt]; // the pthreads used
	/* TODO:
	 * -Subdivide VDIF buffer into blocks.
	 * -Assign thread per block to write
	 * -Return total number of frames written
	/* Launch threads to read data */
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG("\tLaunching threads.");
	#endif
	for (ithread=0; ithread<sgpln->n_sgprt; ithread++)
	{
		//~ thread_result = pthread_create(&(sg_threads[ithread]),NULL,&sgthread_read_block,&(sgpln->sgprt[ithread]));
		//~ if (thread_result != 0)
		//~ {
			//~ perror("Unable to create thread.");
			//~ exit(EXIT_FAILURE);
		//~ }
	}
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG("\tJoining threads.");
	#endif
	/* Join the threads and copy data. */
	for (ithread=0; ithread<sgpln->n_sgprt; ithread++)
	{
		//~ thread_result = pthread_join(sg_threads[ithread],NULL);
		//~ if (thread_result != 0)
		//~ {
			//~ perror("Unable to join thread.");
			//~ exit(EXIT_FAILURE);
		//~ }
	}
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
}

/*
 * Close scatter gather write plan
 */
int close_sg_write_plan(SGPlan *sgplan)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
}

//////////////////////////////////////////////////////////////////////// THREAD IMPLEMENTATIONS
/* 
 * Create an SGInfo instance for the given filename.
 * Arguments:
 *   void *arg -- Pointer to filename string.
 * Returns:
 *   void * -- Pointer to SGInfo instance if the filename produced a 
 *     valid sg_access result (test on SGInfo.smi) or NULL if not.
 * Notes:
 *   This method is suitable for a call via pthread_create.
 */
static void * sgthread_fill_read_sgi(void *arg)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	char *filename = (char *)arg; // filename to try to access
	SGInfo *sgi = (SGInfo *)calloc(sizeof(SGInfo), 1); // SGInfo pointer to return
	sgi->name = NULL;
	sgi->verbose = 0;
	sg_open(filename,sgi);
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		snprintf(_dbgmsg,_DBGMSGLEN,"\tsgi->smi.mmfd = %d",sgi->smi.mmfd);
		DEBUGMSG(_dbgmsg);
		sg_report(sgi,"\tSG Report:");
		DEBUGMSG_LEAVEFUNC;
	#endif
	return (void *)sgi;
}

/*
 * Create SG file and set SGInfo properties accordingly.
 */
static void * sgthread_fill_write_sgi(void *arg)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	char *filename = (char *)arg; // filename to try to access
	SGInfo *sgi = (SGInfo *)malloc(sizeof(SGInfo));
	// Fill in the fields for SGInfo
	sgi->verbose = 0;        /* diagnostic on internal processing */
	sgi->name = (char *)malloc(strlen(filename)*sizeof(char));
	strcpy(sgi->name,filename);
	sgi->total_pkts = 0;     /* total number of packets */
	sgi->sg_version = FILE_VERSION;
	sgi->sg_fht_size = sizeof(struct file_header_tag);
	sgi->sg_wbht_size = sizeof(struct wb_header_tag);
	// Fill in the fields for SGMMInfo
	/* TODO:
	 * -Open file
	umask((mode_t)0);
	sgi->smi->mmfd = open(filename, SG_FILE_WRITE_OPEN_MODE, SG_FILE_PERMISSIONS);
	if (sgi->smi->mmfd == -1)
	{
		perror("Unable to open file.");
		exit(EXIT_FAILURE);
	}
	sgi->smi->size = ???;
	* -Truncate file to initial size
	* -Create mmap
	sgi->smi->start = mmap(NULL, sgi->smi->size, SG_MMAP_WRITE_OPEN_PROTO, SG_MMAP_WRITE_OPEN_MODE, sgi->smi->mmfd, 0);
	if (sgi->smi->start == MAP_FAILED)
	{
		perror("Unable to map file.");
		exit(EXIT_FAILURE);
	}
	sgi->smi->eomem = sgi->smi->start + sgi->smi->size
	sgi->smi->users = 1;
	*/
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
	return (void *)sgi;
}

/*
 * Read one block's worth of VDIF packets from the given SG file.
 * Arguments:
 *   void *arg -- MsgToSGThread by reference that contains a pointer to 
 *     a valid SGInfo instance, and an index specifying which block to
 *     read.
 * Returns:
 *   void *arg -- SGPart pointer associated with the SG file to be read 
 *     from.
 * Notes:
 *   This method is compatible with pthread.
 */
static void * sgthread_read_block(void *arg)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	SGPart *sgprt = (SGPart *)arg;
	uint32_t *start = NULL;
	uint32_t *end = NULL;
	int ii;
	uint32_t *data_buf;// pointer to VDIF data buffer
	uint32_t n_frames;
	
	// check if this is a valid block number
	if (sgprt->iblock < sgprt->sgi->sg_total_blks) 
	{
		//~ start = sg_pkt_by_blk(sgprt->sgi,0,&(sgprt->n_frames),&end);
		start = sg_pkt_by_blk(sgprt->sgi,sgprt->iblock,&(sgprt->n_frames),&end);
		// allocate data storage and copy data to memory
		sgprt->data_buf = (uint32_t *)malloc(sgprt->n_frames*sgprt->sgi->pkt_size);
		if (sgprt->data_buf != NULL)
		{
			memcpy(sgprt->data_buf,start,sgprt->n_frames*sgprt->sgi->pkt_size);
		}
	}
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
	return NULL;
}

/* 
 * Write a block
 */
static void * sgthread_write_block(void *arg)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	SGPart *sgprt = (SGPart *)arg;
	/* TODO: 
	 * -Check number of VDIF frames against current file size, and
	 * resize if necessary, with update to mmap.
	 * -Copy VDIF buffer data to mmap.
	 * -Update the SGInfo fields
	 * -Reset frame counter to zero (or subtract number of written 
	 * frames)
	// To be filled upon first write:
	//~ uint32_t    pkt_size;       /* VDIF packet size */
	//~ uint32_t    pkt_offset;     /* offset into packet */
	//~ uint32_t    read_size;      /* total packet + overhead */
	//~ uint32_t    ref_epoch;      /* reference epoch */
	//~ uint32_t    first_secs;     /* seconds of epoch of first packet */
	//~ uint32_t    first_frame;    /* frame counter of first packet */
	// To be filled with every write
	//~ uint32_t    final_secs;     /* seconds of epoch of final packet */
	//~ uint32_t    final_frame;    /* frame counter of final packet */
	//~ uint32_t    sg_wr_block;    /* standard write block size */
	//~ uint32_t    sg_wr_pkts;     /* pkts in standard write block */
	//~ uint32_t    sg_se_block;    /* ending write block size */
	//~ uint32_t    sg_se_pkts;     /* pkts in ending write block */
	//~ uint32_t    sg_total_blks;  /* total number of blocks */
	// To be left unfilled
	//~ uint32_t    frame_cnt_max;  /* maximum frame counter value seen */
	//~ VDIFsigu    vdif_signature; /* header signature */
	//~ uint32_t    sg_wr_blks_bs;  /* number of write blocks before sb */
	//~ uint32_t    sg_wr_blks_as;  /* number of write blocks after sb */
	//~ uint32_t    sg_wr_pkts_bs;  /* pkts in normal wbs before sb */
	//~ uint32_t    sg_wr_pkts_as;  /* pkts in normal wbs after sb */
	//~ uint32_t    sg_sh_block;    /* starting write block size */
	//~ off_t       sg_sh_blk_off;  /* offset in the file */
	//~ off_t       sg_se_blk_off;  /* offset in the file */
	//~ uint32_t    sg_sh_pkts;     /* pkts in starting write block */
	//~ double      eval_time;      /* diagnostic on file access time */
    #if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
	return NULL;
} 

//////////////////////////////////////////////////////////////////////// TIME ORDERING UTILITIES
/*
 * Comparison method to sort an array of SGInfo elements.
 * Arguments:
 *   const void *a -- SGInfo by reference.
 *   const void *b -- SGInfo by reference.
 * Return:
 *   int - Returns -1 if a < b, 0 if a == b, and 1 if a > b.
 * Notes:
 *   The comparison is done by comparing the timestamp on the first VDIF
 *     frame in the file associated with *a and *b.
 */
int compare_sg_info(const void *a, const void *b)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	SGInfo *sgi_a = (SGInfo *)a;
	SGInfo *sgi_b = (SGInfo *)b;
	int result = sgi_a->first_secs < sgi_b->first_secs ? -1 : sgi_a->first_secs > sgi_b->first_secs;
	if (result == 0)
	{
		result = sgi_a->first_frame < sgi_b->first_frame ? -1 : sgi_a->first_frame > sgi_b->first_frame;
	}
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
	return result;
}

/*
 * Comparison method to sort an array of SGPart elements.
 * Arguments
 *   const void *a -- SGPart by reference.
 *   const void *b -- SGPart by reference.
 * Return:
 *   int -- Returns -1 if a < b, 0 if a == b, and 1 if a > b.
 * Notes:
 *   The comparison is based on the timestamp on the first VDIF frame in
 *     the a->data_buf and b->data_buf.
 */
int compare_sg_part(const void *a, const void *b)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	SGPart *sgprt_a = (SGPart *)a;
	SGPart *sgprt_b = (SGPart *)b;
	/* Seconds since reference epoch. */
	uint32_t secs_inre_a = FIRST_VDIF_SECS_INRE(sgprt_a);
	uint32_t secs_inre_b = FIRST_VDIF_SECS_INRE(sgprt_b);
	/* Data frame number within second */
	uint32_t df_num_insec_a = FIRST_VDIF_DF_NUM_INSEC(sgprt_a);
	uint32_t df_num_insec_b = FIRST_VDIF_DF_NUM_INSEC(sgprt_b);
	
	int result = secs_inre_a < secs_inre_b ? -1 : secs_inre_a > secs_inre_b;
	#if DEBUG_LEVEL >= DEBUG_LEVEL_INFO
		snprintf(_dbgmsg, _DBGMSGLEN, "%d = %d ? %d : %d\n",result,secs_inre_a < secs_inre_b,-1,secs_inre_a > secs_inre_b);
		INFOMSG(_dbgmsg);
	#endif
	if (result == 0)
	{
		result = df_num_insec_a < df_num_insec_b ? -1 : df_num_insec_a > df_num_insec_b;
		#if DEBUG_LEVEL >= DEBUG_LEVEL_INFO
			snprintf(_dbgmsg, _DBGMSGLEN, "%d = %d ? %d : %d\n",result,df_num_insec_a < df_num_insec_b,-1,df_num_insec_a > df_num_insec_b);
			INFOMSG(_dbgmsg);
		#endif
	}
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		snprintf(_dbgmsg,_DBGMSGLEN,"Result = %d\n",result);
		DEBUGMSG(_dbgmsg);
		DEBUGMSG_LEAVEFUNC;
	#endif
	return result;
}

/*
 * Find a contiguous mapping of SGParts in the given SGPlan.
 * Arguments:
 *   SGPlan *sgpln -- SGPlan that contains the SGParts array to order.
 *   int *mapping -- Allocated integer array that will contain the 
 *     ordered mapping so that the first M entries will list the 
 *     contiguous blocks from start to end. This is followed by 
 *     sgpln->n_sgprt-M negative indecies that list blocks that are not
 *     contiguous with this block set.
 * Returns:
 *   int -- The number of contiguous blocks found.
 */
int map_sg_parts_contiguous(SGPlan *sgpln, int *mapping)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	int ii, jj;
	int idx_min_new = 0;
	int tmp_map = 0;
	int return_value = 0;
	int dead_nodes = sgpln->n_sgprt;
	/* Initialize mapping array. Just put indecies in increasing 
	 * magnitude, and make entries for dead nodes negative.
	 */
	for (ii=0; ii<sgpln->n_sgprt; ii++)
	{
		if (sgpln->sgprt[ii].n_frames > 0)
		{
			dead_nodes--;
			mapping[ii] = ii+1;
		}
		else
		{
			mapping[ii] = -(ii+1);
		}
	}
	/* If all nodes dead, just return zero. */
	if (dead_nodes == sgpln->n_sgprt)
	{
		return 0;
	}
	/* Put all dead nodes at the end */
	qsort((void *)mapping, sgpln->n_sgprt, sizeof(int), compare_int_descend);
	/* Sort according to timestamps */
	for (ii=0; ii<sgpln->n_sgprt-dead_nodes; ii++)
	{
		idx_min_new = ii;
		for (jj=ii; jj<sgpln->n_sgprt-dead_nodes; jj++)
		{
			if (compare_sg_part((void *)&(sgpln->sgprt[mapping[jj]-1]),(void *)&(sgpln->sgprt[mapping[idx_min_new]-1])) < 0)
			{
				idx_min_new = jj;
			}
		}
		tmp_map = mapping[ii];
		mapping[ii] = mapping[idx_min_new];
		mapping[idx_min_new] = tmp_map;
	}
	/* Check data continuity, and set index negative if not. */
	for (ii = 0;ii<sgpln->n_sgprt-dead_nodes-1; ii++)
	{
		if (!test_sg_parts_contiguous(&(sgpln->sgprt[mapping[ii]-1]),&(sgpln->sgprt[mapping[ii+1]-1])))
		{
			break;
		}
	}
	return_value = ii+1;
	for (jj=return_value; jj<sgpln->n_sgprt-dead_nodes; jj++)
	{
		mapping[jj] = -mapping[jj];
	}
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		
		printf("Mapping = [");
		for (ii=0; ii<sgpln->n_sgprt; ii++)
		{
			printf("%5d",mapping[ii]);
		}
		printf("]\n");
		DEBUGMSG_LEAVEFUNC;
	#endif
	return return_value;
}

/*
 * Test whether two SG parts are contiguous.
 * Arguments:
 *   SGPart *a -- Pointer to SGPart assumed to contain first data.
 *   SGPart *b -- Pointer to SGPart assumed to contain last data.
 * Returns:
 *   int -- 1 if contiguous, 0 if not.
 * Notes:
 *   Continuinity means that the last VDIF frame in the a->data_buf
 *     and the first VDIF frame in b->data_buf are adjacent in time, 
 *     according to seconds-since-reference-epoch and data-frame-within-
 *     second counters.
 */
int test_sg_parts_contiguous(SGPart *a, SGPart *b)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	/* Seconds since reference epoch. */
	uint32_t secs_inre_a = LAST_VDIF_SECS_INRE(a);
	uint32_t secs_inre_b = FIRST_VDIF_SECS_INRE(b);
	/* Data frame number within second */
	uint32_t df_num_insec_a = LAST_VDIF_DF_NUM_INSEC(a);
	uint32_t df_num_insec_b = FIRST_VDIF_DF_NUM_INSEC(b);
	if (secs_inre_b == secs_inre_a)
	{
		if (df_num_insec_b == df_num_insec_a+1)
		{
			return 1;
		}
	}
	else if (secs_inre_b == secs_inre_a+1)
	{
		if (df_num_insec_a == VDIF_FRAMES_PER_SECOND-1 && df_num_insec_b == 0)
		{
		#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
			DEBUGMSG_LEAVEFUNC;
		#endif
			return 1;
		}
	}
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
	return 0;
}

//////////////////////////////////////////////////////////////////////// MEMORY MANAGEMENT
/*
 * Clear the data buffer in SGPart.
 * Arguments:
 *   SGPart *sgprt -- Pointer to SGPart instance.
 * Return:
 *   void
 */
void clear_sg_part_buffer(SGPart *sgprt)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	sgprt->n_frames = 0;
	if (sgprt->data_buf != NULL)
	{
		free(sgprt->data_buf);
		sgprt->data_buf = NULL;
	}
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
}

/*
 * Free the resources allocated for an SGInfo structure.
 * Arguments:
 *   SGInfo *sgi -- Pointer to allocated SGInfo structure.
 */
void free_sg_info(SGInfo *sgi)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	if (sgi != NULL)
	{
		if (sgi->name != NULL)
		{
			free(sgi->name);
		}
		free(sgi);
	}
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
}

/*
 * Free the resources allocated for an SGPlan structure.
 * Arguments:
 *   SGPlan *sgpln -- Pointer to allocated SGPlan structure.
 */
void free_sg_plan(SGPlan *sgpln)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	int ii;
	for (ii=0; ii<sgpln->n_sgprt; ii++)
	{
		clear_sg_part_buffer(&(sgpln->sgprt[ii]));
		free_sg_info(sgpln->sgprt[ii].sgi);
	}
	free(sgpln->sgprt);
	free(sgpln);
	#if DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
}

//////////////////////////////////////////////////////////////////////// DEBUG UTILITIES
#ifdef DEBUG_LEVEL
void debug_msg(const char *msg, const char *filename, const char *funcname, int linenum)
{
	printf("DEBUG:%s:%d:%s:%s\n",filename,linenum,funcname,msg);
}

void error_msg(const char *msg, const char *filename, const char *funcname, int linenum)
{
	printf("ERROR:%s:%d:%s:%s\n",filename,linenum,funcname,msg);
}

void warning_msg(const char *msg, const char *filename, const char *funcname, int linenum)
{
	printf("WARNING:%s:%d:%s:%s\n",filename,linenum,funcname,msg);
}

void info_msg(const char *msg, const char *filename, const char *funcname, int linenum)
{
	printf("INFO:%s:%d:%s:%s\n",filename,linenum,funcname,msg);
}

void print_sg_part(SGPart *sgprt, const char *label)
{
	printf("%sSGPart 0x%lx:",label,(unsigned long int)sgprt);
	if (sgprt->data_buf != NULL)
	{
		printf(" %u.%u -->> %u.%u",(uint32_t)FIRST_VDIF_SECS_INRE(sgprt),(uint32_t)FIRST_VDIF_DF_NUM_INSEC(sgprt),
			(uint32_t)LAST_VDIF_SECS_INRE(sgprt),(uint32_t)LAST_VDIF_DF_NUM_INSEC(sgprt));
	}
	printf("\n");
	printf("%s\t.iblock = %lu\n",label,(unsigned long int)(sgprt->iblock));
	printf("%s\t.data_buf = 0x%lx\n",label,(unsigned long int)(sgprt->data_buf));
	printf("%s\t.n_frames = %d\n",label,sgprt->n_frames);
}

void print_sg_plan(SGPlan *sgpln, const char *label)
{
	int ii;
	char new_label[strlen(label)+2];
	snprintf(new_label, strlen(label)+2, "\t\t%s",label);
	printf("%sSGPlan 0x%lx:\n",label,(unsigned long int)sgpln);
	for (ii=0; ii<sgpln->n_sgprt; ii++)
	{
		print_sg_part(&(sgpln->sgprt[ii]),new_label);
	}
}
#endif
