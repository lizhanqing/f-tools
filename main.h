#ifndef FQTOOLS_MAIN_H
#define FQTOOLS_MAIN_H

#ifdef __cplusplus
extern "C" {
#endif

	int stk_comp(int argc, char *argv[]);
	int ftools_fqstat(int argc, char *argv[]);
	int ftools_NoNregion(int argc, char *argv[]);
	int stk_fqchk(int argc, char *argv[]);
	int stk_hety(int argc, char *argv[]);
	int stk_gc(int argc, char *argv[]);
	int stk_subseq(int argc, char *argv[]);
	int stk_mutfa(int argc, char *argv[]);
	int stk_mergefa(int argc, char *argv[]);
	int stk_mergepe(int argc, char *argv[]);
	int stk_dropse(int argc, char *argv[]);
	int stk_randbase(int argc, char *argv[]);
	int stk_cutN(int argc, char *argv[]);
	int stk_gap(int argc, char *argv[]);
	int stk_listhet(int argc, char *argv[]);
	int stk_famask(int argc, char *argv[]);
	int stk_trimfq(int argc, char *argv[]);
	int stk_hrun(int argc, char *argv[]);
	int stk_sample(int argc, char *argv[]);
	int stk_seq(int argc, char *argv[]);
	int stk_kfreq(int argc, char *argv[]);
	int stk_rename(int argc, char *argv[]);
	int stk_split(int argc, char *argv[]);
#ifdef __cplusplus
}
#endif

#endif
