#include "tools.h"
#include <stdlib.h>
#include <stdio.h>
#include <pwd.h>
#include <unistd.h>
#include <sys/types.h>
#include <string.h>
#include <math.h>

char *tilde(const char *file)
{
    struct passwd *pass ;
    char buf[512] ;
    const char *remains ;
    char *newfile ;
    
    if(*file != '~') {
        newfile = (char *) malloc(strlen(file)+1) ;
        strcpy(newfile,file) ;
        return(newfile) ;
    }    
    file++ ;
    if(*file == '/') {
        pass = getpwuid(getuid()) ;
        if(pass == NULL) {
            fprintf(stderr,
                    "warning: Can't get password entry for your userid!\n") ;
            --file ;
            newfile = (char *) malloc(strlen(file)+1) ;
            strcpy(newfile,file) ;
            return newfile ;
        }
    } else {
        char *t ;
        
        remains = file - 1 ;
        t = buf ;
        while(*file != '/' && *file != '\0')
          *t++ = *file++ ;
        *t = '\0' ;
        pass = getpwnam(buf) ;
        if(pass == NULL) {
            fprintf(stderr,
                  "warning: Can't get password entry for userid '%s'\n",buf) ;
            newfile = (char *) malloc(strlen(remains)+1) ;
            strcpy(newfile,remains) ;
            return newfile ;
        }
    }
    newfile = (char *)malloc(strlen(file)+strlen(pass->pw_dir)+1) ;
    strcpy(newfile,pass->pw_dir) ;
    strcat(newfile,file) ;
    return(newfile) ;
}

char * fourSigDigs(double num)
{
    static char buf[256] ;
    
    if(num == 0.0)
      snprintf(buf,256,"0.0") ;
    else if(fabs(num) < 1.0|| fabs(num)> 99999.9) {
        if(fabs(num) > 0.009 && fabs(num)<1.0)
          snprintf(buf,256,"%3.3f",num) ;
        else
          snprintf(buf,256,"%3.2e",num) ;
    } else if(fabs(num) > 999.9)
      snprintf(buf,256,"%g",num) ;
    else
      snprintf(buf,256,"%3.1f",num) ;
    return buf ;

}
