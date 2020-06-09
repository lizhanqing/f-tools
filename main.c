#include <stdio.h>
#include <string.h>
#include "main.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.1 (2020)"
#endif

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: fqtools ()\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Zhanqing LI <peaches.li@foxmail.com>\n\n");
	fprintf(stderr, "Usage:   f-tools <command> [options]\n\n");
	fprintf(stderr, "Command:fqstat           fastq statistics\n");
	fprintf(stderr, "        noNregion        fasta noN-region\n");
	fprintf(stderr, "        seq              common transformation of FASTA/Q\n");
	fprintf(stderr, "        comp             get the nucleotide composition of FASTA/Q\n");
	fprintf(stderr, "        sample           subsample sequences\n");
	fprintf(stderr, "        subseq           extract subsequences from FASTA/Q\n");
	fprintf(stderr, "        fqchk            fastq QC (base/quality summary)\n");
	fprintf(stderr, "        mergepe          interleave two PE FASTA/Q files\n");
	fprintf(stderr, "        trimfq           trim FASTQ using the Phred algorithm\n\n");
	fprintf(stderr, "        hety             regional heterozygosity\n");
	fprintf(stderr, "        gc               identify high- or low-GC regions\n");
	fprintf(stderr, "        mutfa            point mutate FASTA at specified positions\n");
	fprintf(stderr, "        mergefa          merge two FASTA/Q files\n");
	fprintf(stderr, "        famask           apply a X-coded FASTA to a source FASTA\n");
	fprintf(stderr, "        dropse           drop unpaired from interleaved PE FASTA/Q\n");
	fprintf(stderr, "        rename           rename sequence names\n");
	fprintf(stderr, "        randbase         choose a random base from hets\n");
	fprintf(stderr, "        cutN             cut sequence at long N\n");
	fprintf(stderr, "        gap              get the gap locations\n");
	fprintf(stderr, "        listhet          extract the position of each het\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage();
	if (strcmp(argv[1], "comp") == 0) return stk_comp(argc-1, argv+1);
	else if (strcmp(argv[1], "fqstat") == 0) return ftools_fqstat(argc-1, argv+1);
	else if (strcmp(argv[1], "noNregion") == 0) return ftools_NoNregion(argc-1, argv+1);
	else if (strcmp(argv[1], "fqchk") == 0) return stk_fqchk(argc-1, argv+1);
	else if (strcmp(argv[1], "hety") == 0) return stk_hety(argc-1, argv+1);
	else if (strcmp(argv[1], "gc") == 0) return stk_gc(argc-1, argv+1);
	else if (strcmp(argv[1], "subseq") == 0) return stk_subseq(argc-1, argv+1);
	else if (strcmp(argv[1], "mutfa") == 0) return stk_mutfa(argc-1, argv+1);
	else if (strcmp(argv[1], "mergefa") == 0) return stk_mergefa(argc-1, argv+1);
	else if (strcmp(argv[1], "mergepe") == 0) return stk_mergepe(argc-1, argv+1);
	else if (strcmp(argv[1], "dropse") == 0) return stk_dropse(argc-1, argv+1);
	else if (strcmp(argv[1], "randbase") == 0) return stk_randbase(argc-1, argv+1);
	else if (strcmp(argv[1], "cutN") == 0) return stk_cutN(argc-1, argv+1);
	else if (strcmp(argv[1], "gap") == 0) return stk_gap(argc-1, argv+1);
	else if (strcmp(argv[1], "listhet") == 0) return stk_listhet(argc-1, argv+1);
	else if (strcmp(argv[1], "famask") == 0) return stk_famask(argc-1, argv+1);
	else if (strcmp(argv[1], "trimfq") == 0) return stk_trimfq(argc-1, argv+1);
	else if (strcmp(argv[1], "hrun") == 0) return stk_hrun(argc-1, argv+1);
	else if (strcmp(argv[1], "sample") == 0) return stk_sample(argc-1, argv+1);
	else if (strcmp(argv[1], "seq") == 0) return stk_seq(argc-1, argv+1);
	else if (strcmp(argv[1], "kfreq") == 0) return stk_kfreq(argc-1, argv+1);
	else if (strcmp(argv[1], "rename") == 0) return stk_rename(argc-1, argv+1);
	else if (strcmp(argv[1], "split") == 0) return stk_split(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;
}
