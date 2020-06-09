#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "main.h"

static char REF[1024];
static char BED[1024];
int ftools_NoNregion(int argc, char * argv[]){
	char Name[1024];
	int c = 0;
	int start = 0;
	int end = 0;
	int nonNlen = 0;
	int totallen = 0;
	while ((c = getopt(argc, argv, "f:b:")) >= 0) {
		switch (c) {
			case 'f': sscanf (optarg, "%s",REF);break;
			case 'b': sscanf (optarg, "%s",BED);break;
		}
	}
	if (optind == 1){
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   f-tools noNregion [options]\n\n");
		fprintf(stderr, "        -f FILE     fasta file\n");
		fprintf(stderr, "        -b FILE     bed file\n");
		fprintf(stderr, "\n");
		return 1;
	}
	FILE *ref = fopen(REF, "r");
	FILE *bed = fopen(BED, "w+");
	if (ref == NULL){
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
	if (bed == NULL){
		fprintf(stderr, "[E::%s] failed to open the output file/stream.\n", __func__);
		return 1;
	}
	while(!feof(ref)){
		char tmp[4869447];
		fscanf(ref,"%s",tmp);
		if(!feof(ref))
		{
			if (tmp[0] == '>' )
			{	
				if(strlen(Name) != 0 && nonNlen >= 1)
				{
					end = start + nonNlen;
					fprintf(bed, "%s\t%d\t%d\n", Name, start, end);
				}
				int a;
				for(a=0;tmp[a]!='\0';a++){
					tmp[a]=tmp[a+1];
					tmp[a]=tmp[a+1];
				}
				strcpy(Name, tmp);
				start = 0;
				end = 0;
				nonNlen = 0;
				totallen = 0;
			}
			else if (tmp[0] == 'A' || tmp[0] == 'G' || tmp[0] == 'C' || tmp[0] == 'T')
			{
				int i;
				for(i = 0; i < strlen(tmp); i++)
				{
					if(tmp[i] == 'N' || tmp[i] == 'n')
					{
						if(nonNlen >= 1)
						{
							end = start + nonNlen;
							fprintf(bed, "%s\t%d\t%d\n", Name, start, end);
						}
						nonNlen = 0;
					}
					else if(tmp[i] == 'A' || tmp[i] == 'C' || tmp[i] == 'G' || tmp[i] == 'T')
					{
						if(nonNlen == 0)
						{
							start = i + totallen;
						}
						nonNlen++;
					}
					else
					{
						if(nonNlen >= 1)
						{
							end = start + nonNlen;
							fprintf(bed, "%s\t%d\t%d\n", Name, start, end);
						}
						nonNlen = 0;
					}
				}
				totallen += strlen(tmp);
			}
		}
	}
	if(strlen(Name) != 0 && nonNlen >= 1)
	{
		end = start + nonNlen;
		fprintf(bed, "%s\t%d\t%d\n", Name, start, end);
	}
	fclose(ref);
	fclose(bed);
	return 0;
}
