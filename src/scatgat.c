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

/* Reliance on the frame rate is ignored to make the code more portable.
 * This means we cannot check data continuity across 1-second 
 * boundaries. */
//~ /* Measures the number of VDIF frames per second and determines where
 //~ * the VDIFHeader.df_num_insec wraps to zero.
 //~ */
//~ #define VDIF_FRAMES_PER_SECOND 125000

/* This defines the initial memory mapped output file sizes. Guess that
 * we are recording at 1GB/s for 300s. Number needs to be divided by 
 * the number of scatter-gather files created to get the size per 
 * scatter-gather file. */
#define INITIAL_SIZE_IN_BLOCKS 1000
/* Defines by how many blocks an SG file size is incremented whenever 
 * resize is necessary */
#define GROWTH_SIZE_IN_BLOCKS 1000

#define LAST_VDIF_SECS_INRE(a) ((VDIFHeader *)(&(a->data_buf[(a->n_frames-1)*(a->sgi->pkt_size)/sizeof(uint32_t)])))->w1.secs_inre
#define FIRST_VDIF_SECS_INRE(a) ((VDIFHeader *)(a->data_buf))->w1.secs_inre
#define LAST_VDIF_DF_NUM_INSEC(a) ((VDIFHeader *)(&(a->data_buf[(a->n_frames-1)*(a->sgi->pkt_size)/sizeof(uint32_t)])))->w2.df_num_insec
#define FIRST_VDIF_DF_NUM_INSEC(a) ((VDIFHeader *)(a->data_buf))->w2.df_num_insec

/* File permissions with which scatter-gather files are created. */ 
#define SG_FILE_PERMISSIONS (S_IWUSR | S_IRUSR | S_IWGRP | S_IRGRP | S_IROTH)
#define SG_FILE_WRITE_OPEN_MODE (O_RDWR|O_TRUNC|O_CREAT)
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
static void * sgthread_fill_write_sgi(void *arg);
static void * sgthread_write_block(void *arg);

/* Handles writing and resizing of files through mmap */
int first_write_sg_plan(SGPlan *sgpln);
int write_to_sg(SGInfo *sgi, const void *src, size_t n);
int resize_to_sg(SGInfo *sgi, off_t new_size);

/* Memory management */
void clear_sg_part_buffer(SGPart *sgprt);
void free_sg_info(SGInfo *sgi);
void init_sg_part(SGPart *sgprt, const SGInfo *sgi);
void init_sg_info(SGInfo *sgi, const char *filename);

/* Misc checks */
int first_write_sgplan(SGPlan *sgpln);

//////////////////////////////////////////////////////////////////////// SCATTER GATHER READING
/*
 * Create an SGPlan instance in read-mode.
 * Arguments:
 *   SGPlan **sgpln -- Address of SGPlan pointer to allocate memory.
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
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
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
	SGInfo *sgi_tmp;// = (SGInfo *)calloc(sizeof(SGInfo),1); 
	/* Step through all modules and disks, and access files that 
	 * match the pattern.
	 */
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG("\tLaunching threads.");
	#endif
	for (imod=0; imod<n_mod; imod++)
	{
		#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
			snprintf(_dbgmsg,_DBGMSGLEN,"\t\tmod[%d] = %d",imod,mod_list[imod]);
			DEBUGMSG(_dbgmsg);
		#endif
		for (idisk=0; idisk<n_disk; idisk++)
		{
			#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
				snprintf(_dbgmsg,_DBGMSGLEN,"\t\t\tdisk[%d] = %d",idisk,disk_list[idisk]);
				DEBUGMSG(_dbgmsg);
			#endif
			ithread = imod*n_disk + idisk;
			snprintf(filename[ithread],PATH_MAX,fmtstr,mod_list[imod],disk_list[idisk],pattern);
			#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_INFO
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
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG("\tJoining threads.");
	#endif
	for (imod=0; imod<n_mod; imod++)
	{
		#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
			snprintf(_dbgmsg,_DBGMSGLEN,"\t\tmod[%d] = %d",imod,mod_list[imod]);
			DEBUGMSG(_dbgmsg);
		#endif
		for (idisk=0; idisk<n_disk; idisk++)
		{
			#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
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
				#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_INFO
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
		#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
			DEBUGMSG_LEAVEFUNC;
		#endif
		return 0;
	}
	/* Allocate space for storing valid SGInfos. */
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG("\tAllocating buffer space and copying SGInfo.");
	#endif
	
	/* Sort SGInfo array according to second / frames */
	qsort((void *)sgi_buf, valid_sgi, sizeof(SGInfo), compare_sg_info);
	/* Allocate memory for SGPlan */
	*sgpln = (SGPlan *)malloc(sizeof(SGPlan));
	(*sgpln)->sgm = SCATGAT_MODE_READ;
	(*sgpln)->sgprt = (SGPart *)malloc(sizeof(SGPart)*valid_sgi);
	for (itmp=0; itmp<valid_sgi; itmp++)
	{
		// Replaced with initialization function
		//~ // Allocate memory for SGInfo and copy
		//~ (*sgpln)->sgprt[itmp].sgi = (SGInfo *)malloc(sizeof(SGInfo));
		//~ memcpy((*sgpln)->sgprt[itmp].sgi, &(sgi_buf[itmp]), sizeof(SGInfo));
		//~ // Initialize block counter, VDIF buffer, and frame counter
		//~ (*sgpln)->sgprt[itmp].iblock = 0;
		//~ (*sgpln)->sgprt[itmp].data_buf = NULL;
		//~ (*sgpln)->sgprt[itmp].n_frames = 0;
		init_sg_part(&((*sgpln)->sgprt[itmp]),&(sgi_buf[itmp]));
	}
	(*sgpln)->n_sgprt = valid_sgi;
	/* Done with the temporary buffer, free it. */
	free(sgi_buf);
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
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
 *     end no frames could be read, and -1 on error.
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
 *     the contiguous flow. However, data continuity is NOT checked 
 *     between consecutive calls to this method.
 *   Block counter for each SGPart is updated if frames where read from
 *     that file.
 */
int read_next_block_vdif_frames(SGPlan *sgpln, uint32_t **vdif_buf)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
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
	
	/* Check if read mode */
	if (sgpln->sgm != SCATGAT_MODE_READ)
	{
		fprintf(stderr,"Trying to read from non-read-mode SGPlan.\n");
		return -1;
	}
	
	/* Launch threads to read data */
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
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
			/* Remove the frames_estimate increment here ...
			frames_estimate += sgpln->sgprt[ithread].sgi->sg_wr_pkts;
			... and put it outside the thread launch loop. If 
			newly read data maps continuously to already read
			data, then we underestimate the number of frames
			that we'll need to put in the output buffer. */
		}
		/* As per above, increment frames_estimate here. */
		frames_estimate += sgpln->sgprt[ithread].sgi->sg_wr_pkts;

	}
	/* Create storage buffer. Assume that the number of frames read
	 * is always smaller than or equal to the number of estimated frames
	 */
	*vdif_buf = (uint32_t *)malloc(frames_estimate*frame_size);
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
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
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		print_sg_plan(sgpln,"\t");
	#endif
	
/***********************************************************************
 * This part of the code checks continuity of data across block 
 * boundaries.
 */
	n_contiguous_blocks = map_sg_parts_contiguous(sgpln, mapping);
	//~ printf("n_contiguous_blocks = %d\n",n_contiguous_blocks);
	if (n_contiguous_blocks == 0)
	{
		printf("No contiguous blocks found.\n");
		return 0;
	}
	for (isgprt=0; isgprt<n_contiguous_blocks; isgprt++)
	{
		//~ printf("memcpy %d\n",isgprt);
		memcpy((void *)(*vdif_buf + frames_read*frame_size/sizeof(uint32_t)),
				(void *)(sgpln->sgprt[mapping[isgprt]-1].data_buf),sgpln->sgprt[mapping[isgprt]-1].n_frames*frame_size);
		frames_read += sgpln->sgprt[mapping[isgprt]-1].n_frames;
		clear_sg_part_buffer(&(sgpln->sgprt[mapping[isgprt]-1]));
	}
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		snprintf(_dbgmsg,_DBGMSGLEN,"Found %d contiguous blocks\n",n_contiguous_blocks);
		DEBUGMSG(_dbgmsg);
	#endif
/*
 */
//~ /***********************************************************************
 //~ * This part of the code ignores discontinuities of data across block 
 //~ * boundaries.
 //~ */
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
//~ /*
//~  ***********************************************************************/
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		print_sg_plan(sgpln,"\t");
		DEBUGMSG_LEAVEFUNC;
	#endif
	/* Return the frame count. */
	return frames_read;
}

/*
 * Read one block's worth of VDIF frames from a group of SG files.
 * Arguments:
 *   SGPlan *sgpln -- SGPlan instance created in read-mode.
 *   int n_sgi -- Number of SGInfo elements in the array.
 *   off_t iblock -- The block index.
 *   uint32_t **vdif_buf -- Address of a pointer which can be used to
 *     store the location of the VDIF buffer created and filled.
 * Returns:
 *   int -- The number of VDIF frames contained in the buffer, zero if
 *     no frames read, and -1 on error.
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
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	int ithread; // thread counter
	int thread_result; // result of calls to pthread methods
	pthread_t sg_threads[sgpln->n_sgprt]; // the pthreads used
	
	int frames_estimate = 0; // estimate the size of buffer to create
	int frames_read = 0; // count the number of frames received
	int frame_size = sgpln->sgprt[0].sgi->pkt_size; // size of a frame
	
	/* Check if read mode */
	if (sgpln->sgm != SCATGAT_MODE_READ)
	{
		fprintf(stderr,"Trying to read from non-read-mode SGPlan.\n");
		return -1;
	}
	
	/* Launch threads to read data */
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
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
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
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
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
	/* Return the frame count. */
	return frames_read;
}

/*
 * Close scatter gather read plan.
 * Arguments:
 *   SGPlan *sgpln -- Pointer to SGPlan opened in read mode.
 * Return:
 *   void
 */
void close_sg_read_plan(SGPlan *sgpln)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	/* Check if read mode */
	if (sgpln->sgm != SCATGAT_MODE_READ)
	{
		fprintf(stderr,"Cannot close non-read-mode SGPlan as read-mode.\n");
	}
	int ii;
	for (ii=0; ii<sgpln->n_sgprt; ii++)
	{
		sg_close(sgpln->sgprt[ii].sgi);
	}
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
}


//////////////////////////////////////////////////////////////////////// SCATTER GATHER WRITING
/*
 * Create a write-mode SGPlan instance.
 * Arguments:
 *   SGPlan **sgpln -- Address of SGPlan pointer to allocate memory.
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
 *   int -- Number of SGInfo instances (SG files created)
 * Notes:
 *   For each SG file created an SGPart element is stored in SGPlan.
 */
int make_sg_write_plan(SGPlan **sgpln, const char *pattern, 
					const char *fmtstr, int *mod_list, int n_mod, 
					int *disk_list, int n_disk)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	int itmp; // just a counter
	int idisk, imod; // disk, module counters
	char filename[n_mod*n_disk][PATH_MAX]; // full filename searched for
	int ithread; // thread counter
	int thread_result; // return result for pthread methods
	pthread_t sg_threads[n_mod*n_disk]; // pthreads to do filling
	int valid_sgi = 0; // number of valid SG files created
	/* Temporary store for allocated SGInfo instances. */
	SGInfo *sgi_tmp;
	/* Temporary store for SGPart instances. */
	SGPart sgprt_tmp[n_mod*n_disk];
	/* Step through all modules and disks, and access files that 
	 * match the pattern.
	 */
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG("\tLaunching threads.");
	#endif
	for (imod=0; imod<n_mod; imod++)
	{
		#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
			snprintf(_dbgmsg,_DBGMSGLEN,"\t\tmod[%d] = %d",imod,mod_list[imod]);
			DEBUGMSG(_dbgmsg);
		#endif
		for (idisk=0; idisk<n_disk; idisk++)
		{
			#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
				snprintf(_dbgmsg,_DBGMSGLEN,"\t\t\tdisk[%d] = %d",idisk,disk_list[idisk]);
				DEBUGMSG(_dbgmsg);
			#endif
			ithread = imod*n_disk + idisk;
			snprintf(filename[ithread],PATH_MAX,fmtstr,mod_list[imod],disk_list[idisk],pattern);
			#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_INFO
				snprintf(_dbgmsg,_DBGMSGLEN,"\t\t\tCreating file '%s'.",filename[ithread]);
				INFOMSG(_dbgmsg);
			#endif
			thread_result = pthread_create(&(sg_threads[ithread]), NULL, &sgthread_fill_write_sgi, filename[ithread]);
			if (thread_result != 0)
			{
				perror("Unable to launch thread.");
				exit(EXIT_FAILURE);
			}
		}
	}
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG("\tLaunching threads.");
	#endif
	for (imod=0; imod<n_mod; imod++)
	{
		#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
			snprintf(_dbgmsg,_DBGMSGLEN,"\t\tmod[%d] = %d",imod,mod_list[imod]);
			DEBUGMSG(_dbgmsg);
		#endif
		for (idisk=0; idisk<n_disk; idisk++)
		{
			#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
				snprintf(_dbgmsg,_DBGMSGLEN,"\t\t\tdisk[%d] = %d",idisk,disk_list[idisk]);
				DEBUGMSG(_dbgmsg);
			#endif
			ithread = imod*n_disk + idisk;
			snprintf(filename[ithread],PATH_MAX,fmtstr,mod_list[imod],disk_list[idisk],pattern);
			#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_INFO
				snprintf(_dbgmsg,_DBGMSGLEN,"\t\t\tCreating file '%s'.",filename[ithread]);
				INFOMSG(_dbgmsg);
			#endif
			thread_result = pthread_join(sg_threads[ithread], (void *)&sgi_tmp);
			if (thread_result != 0)
			{
				perror("Unable to launch thread.");
				exit(EXIT_FAILURE);
			}
			if (sgi_tmp != NULL)
			{
				init_sg_part(&(sgprt_tmp[valid_sgi++]), sgi_tmp);
			}
			else
			{
				#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_INFO
					snprintf(_dbgmsg,_DBGMSGLEN,"\t\t\tUnable to create file '%s'.",filename[ithread]);
					WARNINGMSG(_dbgmsg);
				#endif
			}
		}
	}
	*sgpln = (SGPlan *)malloc(sizeof(SGPlan));
	(*sgpln)->sgm = SCATGAT_MODE_WRITE;
	(*sgpln)->n_sgprt = valid_sgi;
	(*sgpln)->sgprt = (SGPart *)malloc(sizeof(SGPart)*valid_sgi);
	memcpy((*sgpln)->sgprt, sgprt_tmp, sizeof(SGPart)*valid_sgi);
	/* Free the temporary SGInfo resources, but DO NOT free
	 * the malloc'ed NAME to which we still keep a pointer.
	 */
	free(sgi_tmp);
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
	return valid_sgi;
}

/*
 * Write VDIF frames to SG files.
 * Arguments:
 *   SGPlan *sgpln -- Pointer to write-mode SGPlan instance.
 *   uint32_t *vdif_buf -- Buffer that contains VDIF data to write
 *   int n_frames -- Number of frames to write from buffer.
 * Returns
 *   int -- The number of frames written if all writes were successful.
 */
int write_vdif_frames(SGPlan *sgpln, uint32_t *vdif_buf, int n_frames)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	int ii;
	int ithread; // thread counter
	int ipass; // passes through data
	int thread_result; // result of calls to pthread methods
	pthread_t sg_threads[sgpln->n_sgprt]; // the pthreads used
	int frames_per_block;
	int frames_written = 0;
	int first_sg_idx = 0;
	int this_sg_idx;
	int thread_count;
	
	/* Check if read mode */
	if (sgpln->sgm != SCATGAT_MODE_WRITE)
	{
		fprintf(stderr,"Trying to write to non-write-mode SGPlan.\n");
		return -1;
	}
	
	/* If first write, set some properties */
	if (first_write_sg_plan(sgpln))
	{
		// To be filled upon first write:
		VDIFHeader *vdif_header = (VDIFHeader *)vdif_buf;
		for (ii=0; ii<sgpln->n_sgprt; ii++)
		{
			sgpln->sgprt[ii].sgi->pkt_size = vdif_header->w3.df_len * 8;
			sgpln->sgprt[ii].sgi->pkt_offset = sizeof(VDIFHeader);
			sgpln->sgprt[ii].sgi->first_secs = vdif_header->w1.secs_inre;
			sgpln->sgprt[ii].sgi->first_frame = vdif_header->w2.df_num_insec;
			sgpln->sgprt[ii].sgi->ref_epoch = vdif_header->w2.ref_epoch;
		}
	}
	else
	{
		// TODO: Check if incoming packets are valid?
	}
	printf("WBLOCK_SIZE = %ld, pkt_size = %ld\n",(long int)WBLOCK_SIZE,(long int)sgpln->sgprt[0].sgi->pkt_size);
	frames_per_block = WBLOCK_SIZE/sgpln->sgprt[0].sgi->pkt_size;
	/* Find the first SG file that is short */
	for (ithread=1; ithread<sgpln->n_sgprt; ithread++)
	{
		if (sgpln->sgprt[ithread].iblock < sgpln->sgprt[first_sg_idx].iblock)
		{
			first_sg_idx = ithread;
		}
	}
	/* While there is data unwritten, pass through write-cycle */
	while (frames_written < n_frames)
	{
		#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
			DEBUGMSG("\tLaunching threads.");
		#endif
		thread_count = 0;
		for (ithread=0; ithread<sgpln->n_sgprt && frames_written<n_frames; ithread++)
		{
			thread_count++;
			this_sg_idx = (ithread+first_sg_idx) % sgpln->n_sgprt;
			sgpln->sgprt[this_sg_idx].data_buf = vdif_buf + frames_written*(sgpln->sgprt[0].sgi->pkt_size)/sizeof(uint32_t);
			sgpln->sgprt[this_sg_idx].n_frames = (n_frames-frames_written)<frames_per_block ? (n_frames-frames_written) : frames_per_block;
			thread_result = pthread_create(&(sg_threads[ithread]),NULL,&sgthread_write_block,&(sgpln->sgprt[this_sg_idx]));
			if (thread_result != 0)
			{
				perror("Unable to create thread.");
				exit(EXIT_FAILURE);
			}
			frames_written += sgpln->sgprt[this_sg_idx].n_frames;
			#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
				snprintf(_dbgmsg,_DBGMSGLEN,"%d / %d frames written.",frames_written,n_frames);
				DEBUGMSG(_dbgmsg);
			#endif
		}
		#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
			DEBUGMSG("\tJoining threads.");
		#endif
		for (ithread=0; ithread<thread_count; ithread++)
		{
			this_sg_idx = (ithread+first_sg_idx) % sgpln->n_sgprt;
			thread_result = pthread_join(sg_threads[ithread],NULL);
			if (thread_result != 0)
			{
				perror("Unable to join thread.");
				exit(EXIT_FAILURE);
			}
		}
		#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
			DEBUGMSG("Unwritten frames left, repeat write-cycle.");
		#endif
	}
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
	return frames_written;
}

/*
 * Close scatter gather write plan
 * Arguments:
 *   SGPlan *sgpln -- Pointer to SGPlan opened in write-mode.
 * Return:
 *   void
 */
void close_sg_write_plan(SGPlan *sgpln)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	/* Check if write mode plan */
	if (sgpln->sgm != SCATGAT_MODE_WRITE)
	{
		fprintf(stderr,"Cannot close non-write-mode SGPlan as write-mode\n");
	}
	int ii;
	for (ii=0; ii<sgpln->n_sgprt; ii++)
	{
		if (sgpln->sgprt[ii].sgi->smi.size != (sgpln->sgprt[ii].sgi->smi.eomem - sgpln->sgprt[ii].sgi->smi.start))
		{
			if (sgpln->sgprt[ii].sgi->smi.size == 0)
			{
				/* Reset size to original value, to fool sg_close() */
				sgpln->sgprt[ii].sgi->smi.size = sgpln->sgprt[ii].sgi->smi.eomem - sgpln->sgprt[ii].sgi->smi.start;
				sg_close(sgpln->sgprt[ii].sgi);
				/* Then delete the file */
				if (unlink(sgpln->sgprt[ii].sgi->name) == -1)
				{
					perror("Unable to remove empty file.");
				}
			}
			else
			{
				resize_to_sg(sgpln->sgprt[ii].sgi, sgpln->sgprt[ii].sgi->smi.size);
				sg_close(sgpln->sgprt[ii].sgi);
			}
		}
		else
		{
			sg_close(sgpln->sgprt[ii].sgi);
		}
	}
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
}

//////////////////////////////////////////////////////////////////////// THREAD IMPLEMENTATIONS
/* 
 * Create an SGInfo instance for reading for the given filename.
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
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	char *filename = (char *)arg; // filename to try to access
	SGInfo *sgi = (SGInfo *)calloc(sizeof(SGInfo), 1); // SGInfo pointer to return
	sgi->name = NULL;
	sgi->verbose = 0;
	sg_open(filename,sgi);
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		snprintf(_dbgmsg,_DBGMSGLEN,"\tsgi->smi.mmfd = %d",sgi->smi.mmfd);
		DEBUGMSG(_dbgmsg);
		sg_report(sgi,"\tSG Report:");
		DEBUGMSG_LEAVEFUNC;
	#endif
	return (void *)sgi;
}

/* Create an SGInfo instance for writing for the given filename.
 * Arguments:
 *   void *arg -- Pointer to filename string.
 * Returns:
 *   void * -- Pointer to SGInfo instance if the the file could be 
 *     opened for writing and mapped to memory (i.e. valid SGInfo 
 *     instance could be created), or NULL otherwise.
 * Notes:
 *   This method is suitable for a call via pthread_create.
 */
static void * sgthread_fill_write_sgi(void *arg)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	char *filename = (char *)arg; // filename to try to access
	SGInfo *sgi = (SGInfo *)malloc(sizeof(SGInfo));
	// Fill in the fields for SGInfo
	init_sg_info(sgi, filename);
	// Fill in the fields for SGMMInfo
	umask((mode_t)0);
	sgi->smi.mmfd = open(filename, SG_FILE_WRITE_OPEN_MODE, SG_FILE_PERMISSIONS);
	if (sgi->smi.mmfd == -1)
	{
		perror("Unable to open / create file.");
		free_sg_info(sgi);
		return (void *)NULL;
	}
	/* Guess initial file size equal to ideal single write block, 
	 * defined in dplane_proxy.h. This size setting is just temporary,
	 * and needs to be reset to zero to indicate the actual number of
	 * bytes written to file.
	 */
	sgi->smi.size = (off_t)INITIAL_SIZE_IN_BLOCKS*WBLOCK_SIZE; //~ WBLOCK_SIZE + sizeof(struct file_header_tag);
	/* Truncate file to initial size */
	if (ftruncate(sgi->smi.mmfd, sgi->smi.size) == -1)
	{
		perror("Unable to reset file size.");
		free_sg_info(sgi);
		return (void *)NULL;
	}
	/* Create mmap */
	sgi->smi.start = mmap(NULL, sgi->smi.size, SG_MMAP_WRITE_OPEN_PROTO, SG_MMAP_WRITE_OPEN_MODE, sgi->smi.mmfd, 0);
	if (sgi->smi.start == MAP_FAILED)
	{
		perror("Unable to map file.");
		free_sg_info(sgi);
		return (void *)NULL;
	}
	sgi->smi.eomem = sgi->smi.start + sgi->smi.size;
	sgi->smi.users = 1;
	/* And reset file size to number of bytes written, zero */
	sgi->smi.size = 0;
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
	return (void *)sgi;
}

/*
 * Read one block's worth of VDIF packets from the given SG file.
 * Arguments:
 *   void *arg -- SGPart by reference that contains a pointer to 
 *     a valid SGInfo instance, and an index specifying which block to
 *     read.
 * Returns:
 *   void *arg -- NULL
 * Notes:
 *   Upon successful read, the data_buf and n_frames fields of the 
 *     received SGPart instance are filled with VDIF data and a frame
 *     count, respectively.
 *   This method is compatible with pthread.
 */
static void * sgthread_read_block(void *arg)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	SGPart *sgprt = (SGPart *)arg;
	uint32_t *start = NULL;
	uint32_t *end = NULL;
	
	// check if this is a valid block number
	if (sgprt->iblock < sgprt->sgi->sg_total_blks) 
	{
		//~ start = sg_pkt_by_blk(sgprt->sgi,0,&(sgprt->n_frames),&end);
		start = sg_pkt_by_blk(sgprt->sgi,sgprt->iblock,(int *)&(sgprt->n_frames),&end);
		// allocate data storage and copy data to memory
		sgprt->data_buf = (uint32_t *)malloc(sgprt->n_frames*sgprt->sgi->pkt_size);
		if (sgprt->data_buf != NULL)
		{
			memcpy(sgprt->data_buf,start,sgprt->n_frames*sgprt->sgi->pkt_size);
		}
	}
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
	return NULL;
}

/* 
 * Write one block's worht of VDIF packets to the given SG file.
 * Arguments:
 *   void *arg -- SGPart by reference that contains the data to be 
 *     written as well as a reference to the SGInfo for the file to 
 *     write to.
 * Return:
 *   void *arg -- NULL
 * Notes:
 *   The data_buf and n_frames fields in the SGPart instance contain the
 *     VDIF data and frame count, respectively, that should be written 
 *     as the single next block in the SG file.
 */
static void * sgthread_write_block(void *arg)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	SGPart *sgprt = (SGPart *)arg;
	struct file_header_tag fht = { 
					.sync_word = SYNC_WORD, 
					.version = FILE_VERSION, 
					.packet_format = VDIF };
	struct wb_header_tag wbht = { 
					.blocknum = sgprt->iblock,
					.wb_size = sgprt->sgi->pkt_size*sgprt->n_frames + sizeof(struct wb_header_tag)};
	/* If first block, write file header */
	if (sgprt->iblock == 0)
	{
		fht.packet_size = sgprt->sgi->pkt_size;
		fht.block_size = fht.packet_size*(WBLOCK_SIZE/fht.packet_size) + sizeof(struct wb_header_tag);
		if (write_to_sg(sgprt->sgi, (void *)&fht, sizeof(struct file_header_tag)) == -1)
		{
			fprintf(stderr,"Unable to write file header tag to SG in thread.\n");
			return NULL;
		}
	}
	/* Write block header */
	if (write_to_sg(sgprt->sgi, (void *)&wbht, sizeof(struct wb_header_tag)) == -1)
	{
		fprintf(stderr,"Unable to write block header tag to SG in thread.\n");
		return NULL;
	}
	/* Write data */
	if (write_to_sg(sgprt->sgi, (void *)(sgprt->data_buf), sgprt->sgi->pkt_size*sgprt->n_frames) == -1)
	{
		fprintf(stderr,"Unable to write data block to SG in thread.\n");
		return NULL;
	}
	/* Update block counter for this SG file */
	sgprt->iblock++;
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
	return NULL;
}

/*
 * Write to SG file and resize if necessary
 * Arguments:
 *   SGInfo *sgi -- Pointer to SGInfo instance for the file to be 
 *     written to.
 *   const void *src -- Pointer to buffer containing source data.
 *   size_t n -- Number of bytes to be written from the source data.
 * Return:
 *   int -- 0 on success, -1 on failure
 * Notes:
 *   The size of the file is increased as necessary, currently in steps
 *     of WBLOCK_SIZE (see dplane_proxy.h).
 */
int write_to_sg(SGInfo *sgi, const void *src, size_t n)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	/* Check if resize is necessary */
	if (sgi->smi.size+n > (off_t)(sgi->smi.eomem-sgi->smi.start))
	{
		if (resize_to_sg(sgi, (off_t)(sgi->smi.eomem-sgi->smi.start)+(off_t)GROWTH_SIZE_IN_BLOCKS*WBLOCK_SIZE) == -1)
		{
			return -1;
		}
	}
	/* Memcopy */
	memcpy(sgi->smi.start+sgi->smi.size, src, n);
	sgi->smi.size += n;
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
	return 0;
}

/*
 * Resize SG file.
 * Arguments:
 *   SGInfo *sgi -- Pointer to SGInfo instance for the file to be 
 *     resized.
 *   off_t new_size -- The new file size, in bytes.
 * Return:
 *   int -- 0 on success, -1 on failure.
 */
int resize_to_sg(SGInfo *sgi, off_t new_size)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	/* Truncate file to new size */
	if (ftruncate(sgi->smi.mmfd, new_size) == -1)
	{
		perror("Unable to reset file size.");
		return -1;
	}
	/* Re-mmap for non-zero new_size, and unmap if size is zero. */
	if (new_size != (off_t)0)
	{
		sgi->smi.start = mremap(sgi->smi.start, (sgi->smi.eomem - sgi->smi.start), new_size, MREMAP_MAYMOVE);
		if (sgi->smi.start == MAP_FAILED)
		{
			perror("Unable to map file.");
			return -1;
		}
		sgi->smi.eomem = sgi->smi.start + new_size;
	}
	else
	{
		if (munmap(sgi->smi.start, (sgi->smi.eomem - sgi->smi.start)) == -1)
		{
			perror("Unable to unmap file for resize to zero.");
			return -1;
		}
		sgi->smi.start = NULL;
		sgi->smi.eomem = NULL;
	}
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
	return 0;
}


//////////////////////////////////////////////////////////////////////// TIME ORDERING UTILITIES
/*
 * Comparison method to sort an array of integers in reverse order, i.e.
 *   from largest to smallest.
 * Arguments:
 *   const void *a -- Integer by reference
 *   const void *b -- Integer by reference
 * Return:
 *   int -- Returns -1 if a > b, 0 if a == b, and 1 if a < b.
 */
int compare_int_descend(const void *a, const void *b)
{
	int *int_a = (int *)a;
	int *int_b = (int *)b;
	return *int_b < *int_a ? -1 : *int_b > *int_a;
}

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
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	SGInfo *sgi_a = (SGInfo *)a;
	SGInfo *sgi_b = (SGInfo *)b;
	int result = sgi_a->first_secs < sgi_b->first_secs ? -1 : sgi_a->first_secs > sgi_b->first_secs;
	if (result == 0)
	{
		result = sgi_a->first_frame < sgi_b->first_frame ? -1 : sgi_a->first_frame > sgi_b->first_frame;
	}
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
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
	//~ #if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		//~ DEBUGMSG_ENTERFUNC;
	//~ #endif
	SGPart *sgprt_a = (SGPart *)a;
	SGPart *sgprt_b = (SGPart *)b;
	/* Seconds since reference epoch. */
	uint32_t secs_inre_a = FIRST_VDIF_SECS_INRE(sgprt_a);
	uint32_t secs_inre_b = FIRST_VDIF_SECS_INRE(sgprt_b);
	/* Data frame number within second */
	uint32_t df_num_insec_a = FIRST_VDIF_DF_NUM_INSEC(sgprt_a);
	uint32_t df_num_insec_b = FIRST_VDIF_DF_NUM_INSEC(sgprt_b);
	
	int result = secs_inre_a < secs_inre_b ? -1 : secs_inre_a > secs_inre_b;
	//~ #if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		//~ snprintf(_dbgmsg, _DBGMSGLEN, "%d = %d ? %d : %d",result,secs_inre_a < secs_inre_b,-1,secs_inre_a > secs_inre_b);
		//~ DEBUGMSG(_dbgmsg);
	//~ #endif
	if (result == 0)
	{
		result = df_num_insec_a < df_num_insec_b ? -1 : df_num_insec_a > df_num_insec_b;
		//~ #if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
			//~ snprintf(_dbgmsg, _DBGMSGLEN, "%d = %d ? %d : %d",result,df_num_insec_a < df_num_insec_b,-1,df_num_insec_a > df_num_insec_b);
			//~ DEBUGMSG(_dbgmsg);
		//~ #endif
	}
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		snprintf(_dbgmsg,_DBGMSGLEN,"Result = %d (%u.%u ? %u.%u)",result,secs_inre_a,df_num_insec_a,secs_inre_b,df_num_insec_b);
		DEBUGMSG(_dbgmsg);
		//~ DEBUGMSG_LEAVEFUNC;
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
	//~ printf(" -->> map_sg_parts_contiguous\n");
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
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
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		
		printf("Mapping(init) = [");
		for (ii=0; ii<sgpln->n_sgprt; ii++)
		{
			printf("%5d",mapping[ii]);
		}
		printf("]\n");
		DEBUGMSG_LEAVEFUNC;
	#endif
	/* If all nodes dead, just return zero. */
	if (dead_nodes == sgpln->n_sgprt)
	{
		return 0;
	}
	/* Put all dead nodes at the end */
	qsort((void *)mapping, sgpln->n_sgprt, sizeof(int), compare_int_descend);
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		
		printf("Mapping(qsort) = [");
		for (ii=0; ii<sgpln->n_sgprt; ii++)
		{
			printf("%5d",mapping[ii]);
		}
		printf("]\n");
		DEBUGMSG_LEAVEFUNC;
	#endif
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
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		
		printf("Mapping(bsort) = [");
		for (ii=0; ii<sgpln->n_sgprt; ii++)
		{
			printf("%5d",mapping[ii]);
		}
		printf("]\n");
		DEBUGMSG_LEAVEFUNC;
	#endif
	
	/* Check data continuity, and set index negative if not. */
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		snprintf(_dbgmsg,_DBGMSGLEN,"%d / %d dead nodes",dead_nodes,sgpln->n_sgprt);
		DEBUGMSG(_dbgmsg);
	#endif
	for (ii = 0;ii<sgpln->n_sgprt-dead_nodes-1; ii++)
	{
		#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
			snprintf(_dbgmsg,_DBGMSGLEN,"ii = %d",ii);
			DEBUGMSG(_dbgmsg);
		#endif
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
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		
		printf("Mapping(final) = [");
		for (ii=0; ii<sgpln->n_sgprt; ii++)
		{
			printf("%5d",mapping[ii]);
		}
		printf("]\n");
		DEBUGMSG_LEAVEFUNC;
	#endif
	//~ printf(" <<-- map_sg_parts_contiguous\n");
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
 *     and the first VDIF frame in b->data_buf are adjacent or aligned
 *     in time, according to seconds-since-reference-epoch and 
 *     data-frame-within-second counters. The aligned case is needed if
 *     two or more parallel streams are processed simultaneously, in 
 *     which case packets may have duplicate timestamps.
 *   In order to make the code portable for different frame rates, the
 *     reliance on VDIF_FRAMES_PER_SECOND is removed. This means that 
 *     continuity across a 1-second boundary is ignored.
 */
int test_sg_parts_contiguous(SGPart *a, SGPart *b)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	
	// Test if either pointer is NULL
	if (a == NULL || b == NULL) {
		#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
			DEBUGMSG_LEAVEFUNC;
		#endif
		return 0;
	}
	
	/* Seconds since reference epoch. */
	uint32_t secs_inre_a_last = LAST_VDIF_SECS_INRE(a);
	uint32_t secs_inre_a_first = FIRST_VDIF_SECS_INRE(a);
	uint32_t secs_inre_b_first = FIRST_VDIF_SECS_INRE(b);
	uint32_t secs_inre_b_last = LAST_VDIF_SECS_INRE(b);
	/* Data frame number within second */
	uint32_t df_num_insec_a_last = LAST_VDIF_DF_NUM_INSEC(a);
	uint32_t df_num_insec_a_first = FIRST_VDIF_DF_NUM_INSEC(a);
	uint32_t df_num_insec_b_first = FIRST_VDIF_DF_NUM_INSEC(b);
	uint32_t df_num_insec_b_last = LAST_VDIF_DF_NUM_INSEC(b);
	
	/* Set valid range for b start */
	uint32_t min_secs_inre_b_first = secs_inre_a_first;
	uint32_t max_secs_inre_b_first = secs_inre_a_last;
	uint32_t min_df_num_insec_b_first;
	uint32_t max_df_num_insec_b_first;
	
	
	if (secs_inre_a_first == secs_inre_a_last) // If a is contained in single second
	{
		if (secs_inre_b_first == secs_inre_a_last) // If b starts in the same second as that containing a
		{
			if (df_num_insec_b_first >= df_num_insec_a_first && // b cannot start before a
				df_num_insec_b_first <= df_num_insec_a_last+1) // and b should start at least directly after a ends
			{
				#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
					snprintf(_dbgmsg,_DBGMSGLEN,"Contiguous (%u.%u -> %u.%u ? %u.%u -> %u.%u)",
								secs_inre_a_first,df_num_insec_a_first,secs_inre_a_last,df_num_insec_a_last,
								secs_inre_b_first,df_num_insec_b_first,secs_inre_b_last,df_num_insec_b_last);
					DEBUGMSG(_dbgmsg);
					DEBUGMSG_LEAVEFUNC;
				#endif
				return 1;
			}
		} // b cannot start outside the second in which a is contained
	}
	else // a is split across multiple seconds
	{
		if (secs_inre_b_first == secs_inre_a_first) // if b starts in same second as a starts
		{
			if (df_num_insec_b_first >= df_num_insec_a_first) // since a ends in future second, only need to check that b did not start before a
			{
				#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
					snprintf(_dbgmsg,_DBGMSGLEN,"Contiguous (%u.%u -> %u.%u ? %u.%u -> %u.%u)",
								secs_inre_a_first,df_num_insec_a_first,secs_inre_a_last,df_num_insec_a_last,
								secs_inre_b_first,df_num_insec_b_first,secs_inre_b_last,df_num_insec_b_last);
					DEBUGMSG(_dbgmsg);
					DEBUGMSG_LEAVEFUNC;
				#endif
				return 1;
			}
		}
		else if (secs_inre_b_first == secs_inre_a_last) // if b starts in the same second as a ends
		{
			if (df_num_insec_b_first <= df_num_insec_a_last+1) // since a started in past second, only need to check that b starts directly after a ends or earlier
			{
				#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
					snprintf(_dbgmsg,_DBGMSGLEN,"Contiguous (%u.%u -> %u.%u ? %u.%u -> %u.%u)",
								secs_inre_a_first,df_num_insec_a_first,secs_inre_a_last,df_num_insec_a_last,
								secs_inre_b_first,df_num_insec_b_first,secs_inre_b_last,df_num_insec_b_last);
					DEBUGMSG(_dbgmsg);
					DEBUGMSG_LEAVEFUNC;
				#endif
				return 1;
			}
		}
		else if (secs_inre_b_first > secs_inre_a_first && // if b starts in some second between where a starts and ends, that's okay too
					secs_inre_b_first < secs_inre_a_last)
		{
			#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
				snprintf(_dbgmsg,_DBGMSGLEN,"Contiguous (%u.%u -> %u.%u ? %u.%u -> %u.%u)",
							secs_inre_a_first,df_num_insec_a_first,secs_inre_a_last,df_num_insec_a_last,
							secs_inre_b_first,df_num_insec_b_first,secs_inre_b_last,df_num_insec_b_last);
				DEBUGMSG(_dbgmsg);
				DEBUGMSG_LEAVEFUNC;
			#endif
			return 1;
		}
	}
	/* All other cases fail */
	//~ if (secs_inre_b_first == secs_inre_a_last)
	//~ {
		//~ if ((df_num_insec_b_first == df_num_insec_a_last+1) || (df_num_insec_b_first == df_num_insec_a_first))
		//~ {
			//~ #if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
				//~ snprintf(_dbgmsg,_DBGMSGLEN,"Contiguous (%u.%u -> %u.%u ? %u.%u -> %u.%u)",secs_inre_a_first,df_num_insec_a_first,secs_inre_a_last,df_num_insec_a_last,secs_inre_b_first,df_num_insec_b_first,secs_inre_b_last,df_num_insec_b_last);
				//~ DEBUGMSG(_dbgmsg);
				//~ DEBUGMSG_LEAVEFUNC;
			//~ #endif
			//~ #if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
				//~ DEBUGMSG_LEAVEFUNC;
			//~ #endif
			//~ return 1;
		//~ }
	//~ }
	//~ else if (secs_inre_b_first == secs_inre_a_last+1)
	//~ {
		//~ if (df_num_insec_a_last == VDIF_FRAMES_PER_SECOND-1 && df_num_insec_b_first == 0)
		//~ {
		//~ #if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
			//~ DEBUGMSG_LEAVEFUNC;
		//~ #endif
			//~ return 1;
		//~ }
	//~ }
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		snprintf(_dbgmsg,_DBGMSGLEN,"Non-Contiguous (%u.%u -> %u.%u ? %u.%u -> %u.%u)",
						secs_inre_a_first,df_num_insec_a_first,secs_inre_a_last,df_num_insec_a_last,
						secs_inre_b_first,df_num_insec_b_first,secs_inre_b_last,df_num_insec_b_last);
			DEBUGMSG(_dbgmsg);
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
 * Notes:
 *   Free memory pointed to by data_buf field (if not NULL), and set it 
 *     to NULL. Reset frame counter to zero.
 */
void clear_sg_part_buffer(SGPart *sgprt)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	sgprt->n_frames = 0;
	if (sgprt->data_buf != NULL)
	{
		free(sgprt->data_buf);
		sgprt->data_buf = NULL;
	}
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
}

/*
 * Free the resources allocated for an SGInfo structure.
 * Arguments:
 *   SGInfo *sgi -- Pointer to allocated SGInfo structure.
 * Return:
 *   void
 * Notes:
 *   Free memory pointed to by name field (if not NULL), and free the 
 *     memory pointed to by sgi.
 */
void free_sg_info(SGInfo *sgi)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
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
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
}

/*
 * Free the resources allocated for an SGPlan structure.
 * Arguments:
 *   SGPlan *sgpln -- Pointer to allocated SGPlan structure.
 * Return:
 *   void
 * Notes:
 *   Free each SGPart instance associated with this SGPlan, and finally
 *     free the SGPlan.
 *   If this is a read-mode SGPlan, then clear_sg_part_buffer is also
 *     called on each SGPart instance.
 */
void free_sg_plan(SGPlan *sgpln)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	int ii;
	for (ii=0; ii<sgpln->n_sgprt; ii++)
	{
		if (sgpln->sgm == SCATGAT_MODE_READ)
		{
			clear_sg_part_buffer(&(sgpln->sgprt[ii]));
		}
		free_sg_info(sgpln->sgprt[ii].sgi);
	}
	free(sgpln->sgprt);
	free(sgpln);
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
}

/* 
 * Set default values for new SGPart instance, and deep copy given
 * SGInfo.
 * Arguments:
 *   SGPart *sgprt -- Pointer to SGPart instance to initialize.
 *   SGInfo *sgi -- Pointer to SGInfo to use for initialization.
 * Return:
 *   void
 * Notes:
 *   The sgprt->sgi member is allocated new memory and the contents of 
 *     the passed argument sgi are copied to that location.
 */
void init_sg_part(SGPart *sgprt, const SGInfo *sgi)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	sgprt->sgi = (SGInfo *)malloc(sizeof(SGInfo));
	memcpy(sgprt->sgi, sgi, sizeof(SGInfo));
	sgprt->iblock = 0;
	sgprt->data_buf = NULL;
	sgprt->n_frames = 0;
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif	
}

/*
 * Set default values for new SGInfo instance.
 * Arguments:
 *   SGInfo *sgi -- Pointer to SGInfo instance.
 *   const char *filename -- Filename associated with new SGInfo.
 * Return:
 *   void
 * Notes:
 *   The sgi->name member is allocated new memory and the contents of
 *     the passed filename are copied into that location.
 */
void init_sg_info(SGInfo *sgi, const char *filename)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	sgi->name = (char *)malloc(strlen(filename)*sizeof(char));
	strcpy(sgi->name,filename);
	sgi->verbose = 0;
	sgi->total_pkts = 0;
	sgi->sg_version = FILE_VERSION;
	sgi->sg_fht_size = sizeof(struct file_header_tag);
	sgi->sg_wbht_size = sizeof(struct wb_header_tag);
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif	
} 

//////////////////////////////////////////////////////////////////////// MISC CHECKS
/*
 * Check if this is the first write operation to a write-mode SGPlan.
 * Arguments:
 *   SGPlan *sgpln -- Pointer to allocated SGPlan opened in write-mode.
 * Return:
 *   int -- 1 if true, 0 if false.
 * Notes:
 *   First write assumed when the block counters for all SGPart 
 *     instances are equal to zero.
 */
int first_write_sg_plan(SGPlan *sgpln)
{
	#ifdef DEBUG_LEVEL
		char _dbgmsg[_DBGMSGLEN];
	#endif
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_ENTERFUNC;
	#endif
	int ii;
	for (ii=0; ii<sgpln->n_sgprt; ii++)
	{
		if (sgpln->sgprt[ii].iblock > 0)
		{
			return 0;
		}
	}
	#if defined(DEBUG_LEVEL) && DEBUG_LEVEL >= DEBUG_LEVEL_DEBUG
		DEBUGMSG_LEAVEFUNC;
	#endif
	return 1;
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
