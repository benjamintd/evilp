#include "exact.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(void){char line[32768];while(fgets(line,sizeof line,stdin)){line[strcspn(line,"\r\n")]=0;char converted[32768];strcpy(converted,line);for(char *s=converted;*s;s++)if(*s=='d'||*s=='D')*s='e';double rounded=strtod(converted,NULL);LabInterval q=lab_source_number_interval(line,rounded);printf("%a %a %a\n",rounded,q.lower,q.upper);}return 0;}
