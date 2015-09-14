


#include "launcher.h"


 
 int launch (const char* script, const char* arg_vec[], int num_tasks) {
    char error[DRMAA_ERROR_STRING_BUFFER];
    int errnum = 0;
    drmaa_job_template_t *jt = NULL;
 
   errnum = drmaa_init (NULL, error, DRMAA_ERROR_STRING_BUFFER);
 
    if (errnum != DRMAA_ERRNO_SUCCESS) {
       fprintf (stderr, "Could not initialize the DRMAA library: %s\n", error);
       return 1;
    }
 
    errnum = drmaa_allocate_job_template (&jt, error, DRMAA_ERROR_STRING_BUFFER);
 
    if (errnum != DRMAA_ERRNO_SUCCESS) {
       fprintf (stderr, "Could not create job template: %s\n", error);
    }
    else {
       errnum = drmaa_set_attribute (jt, DRMAA_REMOTE_COMMAND, script,
                                     error, DRMAA_ERROR_STRING_BUFFER);
 
       if (errnum != DRMAA_ERRNO_SUCCESS) {
          fprintf (stderr, "Could not set attribute \"%s\": %s\n",
                   DRMAA_REMOTE_COMMAND, error);
       }
       else {
        
          
          errnum = drmaa_set_vector_attribute (jt, DRMAA_V_ARGV, arg_vec, error,
                                               DRMAA_ERROR_STRING_BUFFER);
       }
       
       if (errnum != DRMAA_ERRNO_SUCCESS) {
          fprintf (stderr, "Could not set attribute \"%s\": %s\n",
                   DRMAA_REMOTE_COMMAND, error);
       }else {
         drmaa_job_ids_t *ids = NULL;
 
          errnum = drmaa_run_bulk_jobs (&ids, jt, 1, num_tasks, 1, error, DRMAA_ERROR_STRING_BUFFER);
 
          if (errnum != DRMAA_ERRNO_SUCCESS) {
             fprintf (stderr, "Could not submit job: %s\n", error);
          }else {
             char jobid[DRMAA_JOBNAME_BUFFER];
		const char *jobids[2] = {DRMAA_JOB_IDS_SESSION_ALL, NULL};
             while (drmaa_get_next_job_id (ids, jobid, DRMAA_JOBNAME_BUFFER) == DRMAA_ERRNO_SUCCESS) {
                printf ("A job task has been submitted with id %s\n", jobid);
             }
          
 	     errnum = drmaa_synchronize (jobids, DRMAA_TIMEOUT_WAIT_FOREVER,
                                         1, error, DRMAA_ERROR_STRING_BUFFER);
             
             if (errnum != DRMAA_ERRNO_SUCCESS) {
                fprintf (stderr, "Could not wait for jobs: %s\n", error);
             }else {
                printf ("All job tasks have finished.\n");
             }
         } 
           drmaa_release_job_ids (ids);
       } /* else */
 
       errnum = drmaa_delete_job_template (jt, error, DRMAA_ERROR_STRING_BUFFER);
 
       if (errnum != DRMAA_ERRNO_SUCCESS) {
          fprintf (stderr, "Could not delete job template: %s\n", error);
       }
    } /* else */
 
    errnum = drmaa_exit (error, DRMAA_ERROR_STRING_BUFFER);
 
    if (errnum != DRMAA_ERRNO_SUCCESS) {
       fprintf (stderr, "Could not shut down the DRMAA library: %s\n", error);
      return 1;
   }

    return 0;
 }
