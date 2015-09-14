#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]){

	if (argc<2){ fprintf(stdout,"Gimme something!\n");}
	else {fprintf(stdout,"Your string: %s\n",argv[1]);}


}
