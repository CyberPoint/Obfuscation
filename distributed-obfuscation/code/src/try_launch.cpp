#include "launcher.h"
#include <stdlib.h>
#include<stdio.h>


int main(int argc, char* argv[])
{
if (argc<2){
fprintf(stdout,"Gimme something!\n");
}
else{ 
const char* arg_vec[2];
char* home = getenv("HOME");
char* dirstring = (char*)malloc(101);
snprintf(dirstring,100,"%s/obfuscation_grid/code/circuits/point-4.acirc.obf.8",home);
arg_vec[0] = dirstring;
arg_vec[1] = NULL;

launch("/home/mroberts/obfuscation_grid/code/src/launch_zen.sh",arg_vec,4);
free(dirstring);}
return 1;

}
