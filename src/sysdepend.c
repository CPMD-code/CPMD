/*-- IBM ------------------------------------------*/
#if defined(__IBM) && !defined(__OSX) && !defined(__OSX_IFC) && !defined(__PWRLinux)
/*-- NOMEMINFO ------------------------------------------*/
#elif defined(__NOMEMINFO)
#include <stdlib.h>
int  printmemsize(long *value)
{
  int st =  sbrk(0);
  return((int) st);
}

/*-- SUN ------------------------------------------*/
#elif defined(__SUN)
#define PRINTMEMSIZE
/*-- SUN ------------------------------------------*/


/*-- SGI ------------------------------------------*/
#elif defined(__SGI)
#define PRINTMEMSIZE
/*-- SGI ------------------------------------------*/

/*--Hewlet Packard --------------------------------*/
#elif defined(__HP)
#include <string.h>
#include <sys/param.h>
#include <sys/pstat.h>
#include <unistd.h>

#define ERROR { *vdata=0; *vstack=0; return; };

/* Give memory used by process */
void memme_(double *vdata, double *vstack)
{
  int pstat_getstatic(struct pst_static *buf, size_t elemsize, 
                      size_t elemcount,int index);
  int pstat_getproc(struct pst_status *buf, size_t elemsize, 
                      size_t elemcount,int index);     
  pid_t getpgid ();

  size_t elemcount;
  int    index;
  int    pagesize, chk;
  struct pst_static buf_static;
  struct pst_status buf_status;

  elemcount=1;
  index=0;
  /* get pagesize */
  chk=pstat_getstatic(&buf_static,(size_t) sizeof(buf_static),elemcount,index);
  if (chk == -1) ERROR
  pagesize = (int) buf_static.page_size;

  elemcount=0;
  index=(int) getpid();
  if (index <= 0) ERROR
  chk=pstat_getproc(&buf_status,(size_t) sizeof(buf_status),elemcount,index); 
  if (chk == -1) ERROR

  *vdata  = (double) buf_status.pst_dsize;
  *vstack = (double) buf_status.pst_ssize;

  *vdata  = *vdata /1024*(pagesize/1024);
  *vstack = *vstack/1024*(pagesize/1024);

  return;
}
/*-- Hewlet Packard -------------------------------*/


/*-- Fujitsu VPP ----------------------------------*/
#elif defined(_vpp_)
#define PRINTMEMSIZE
/*-- Linux ----------------------------------------*/
#elif (defined(__Linux) || defined(__ALPHALINUX) || defined(__PWRLinux) || defined(__WINNT))

/*
 * redirect stdout from c to /dev/null. needed for path integrals
 */

#include <stdio.h>

void silentstdout_(void)
{
#if !defined(__WINNT)
        fopen("/dev/null", "w+");
#endif
        return;
}

/* 
 * get the current resident set size from the /proc filesystem.  this is
 * _much_ faster and far more stable and portable across linux kernel
 * versions than parsing some obscure ps output which sometimes changes
 * with every new moon.
 * Cray XT3 does not have /proc, so we cannot use this.
 * 
 * added 07/2001 by axel.kohlmeyer@theochem.rub.de
 * changed from RSS to VSIZE 04/2004 by axel.kohlmeyer@theochem.rub.de
 * 10/2005 disabled for Cray XT3, alternative at the top of the file
 * by akohlmey@cmm.upenn.edu.
 */

#include <stdio.h>
#include <string.h>

#if defined(__ALPHALINUX)
# define PGSIZE 8192U
#else
# define PGSIZE 4096U
#endif
#define MYBUFSIZE 1024

#define UNKNOWN "(unknown)"
void printmemsize_(void)
{
    FILE *pf;
    unsigned int sz, rss;
    char buf[MYBUFSIZE], *ptr;

    rss = sz = 0;
    memset(buf, 0, MYBUFSIZE);
    strncpy(buf,UNKNOWN,strlen(UNKNOWN));

    /* newer linux kernels support the more convenient /proc/self/status
     * file, so we try read that one first. */
    pf = fopen("/proc/self/status", "r");
    if (pf) {
        fread(buf, MYBUFSIZE, 1, pf);
        buf[MYBUFSIZE - 1] = 0;
        fclose(pf);

        /* get the total size */
        ptr = strstr(buf, "VmSize");
        if (ptr != NULL) {
            sscanf(ptr, "%*s %u", &sz);
        }
        
        /* get the resident set size */
        ptr = strstr(buf, "VmRSS");
        if (ptr != NULL) {
            sscanf(ptr, "%*s %u", &rss);
        }
        
        printf("%7u/%7u", rss, sz);
    } else {
        pf = fopen("/proc/self/stat", "r");
        if (pf) {
            fread(buf, MYBUFSIZE, 1, pf);
            buf[MYBUFSIZE - 1] = 0;
            fclose(pf);
            sscanf(buf, "%*d %*s %*c %*d %*d %*d %*d %*d %*u %*u "
                   "%*u %*u %*u %*d %*d %*d %*d %*d %*d %*u %*u "
                   "%*d %u %u", &sz, &rss);
            printf("%7u/%7u", PGSIZE / 1024U * rss, sz / 1024U);
        } else {
            printf("(unknown)");
        }
    }
    (void)fflush(stdout);
    return;
}
#undef MYBUFSIZE
#undef PGSIZE

/*-- Other computers ------------------------------*/
#elif !defined(__NEC) && !defined (__SGI) && !defined(__SUN) && !defined(__OSX) && !defined(__OSX_IFC)

#endif
/*-- Other computers ------------------------------*/

/*-- PRINTMEMSIZE -----------------------------------*/
#if defined(PRINTMEMSIZE)

#include <sys/fcntl.h>   /* define O_RDONLY */
#include <sys/procfs.h>  /* define PIOCPSINFO and PIOCGETPR */

#include <limits.h>      /* PID_MAX */
#include <string.h>      /* strlen */     
#include <math.h>        /* log10 */
#include <stdio.h>       /* printf */

#include <stdlib.h>      /* malloc */

/* Print Memory used by process */
void printmemsize_( void ) 
{ 
  int  fid,imax,i,j,namelen,n;
  size_t l;
  long vrtsize;
  char sproc[] = "/proc/";
  char *pname;
  prpsinfo_t info; 
  
  imax = (int)log10( (double)PID_MAX)+1;
  l = strlen(sproc);
  pname = malloc( (imax+(int)l)*sizeof(char) );
  if( ! pname) {
    (void)fprintf(stdout,"\n printmemsize! pname not allocated!\n");
    (void)fflush(stdout);
    return;
  }
  n = getpid();
  (void)sprintf (pname,"/proc/%d",n);
  /* First, we try %d */
  fid = open (pname, O_RDONLY);
  if (fid == -1) {
    /* Failed, we try %n.nd */
    (void)sprintf(pname,"%s",sproc);
    j = (int)log10( (double)n)+1;
    for (i=0; i<imax-j; i++) {
      (void)sprintf(&pname[l+i],"0");
    }
    (void)sprintf(&pname[l+imax-j],"%d",n);
    fid = open (pname, O_RDONLY);
  }
  if( ioctl(fid, PIOCPSINFO, (char *)&info) ) {
    (void)fprintf(stdout,"\n printmemsize! error in ioctl\n");
    (void)fflush(stdout);
    return;    
  }
  close (fid);
  /* We use the image size */
  vrtsize = info.pr_size*getpagesize()/1000;
  (void)fprintf(stdout,"%7ld",vrtsize);
  (void)fflush(stdout);
  /* Deallocation */
  free(pname);
} /* printmemsize_ */
#endif

/*-- END PRINTMEMSIZE --*/

