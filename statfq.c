#include <zlib.h>  
#include <stdio.h>
#include <string.h>  
#include "main.h" 
#include "kseq.h"  
KSEQ_INIT(gzFile, gzread)
int ftools_fqstat(int argc, char *argv[])  
{
	gzFile fp;  
	kseq_t *seq;
	int l;
	if (argc <= 1 || argc >= 4) {  
		fprintf(stderr, "Usage: %s <fq1> \n", argv[0]);
		fprintf(stderr, " or  : %s <fq1> <fq2>\n", argv[0]);
		return 1;  
	}else if (argc  == 2 ){
		int se_read1_length;
		long se_read1_reads = 0; long se_read1_bases  = 0; long se_read1_Q10_base = 0; long se_read1_Q20_base = 0; long se_read1_Q30_base = 0; long se_read1_gc_base = 0; long se_read1_A_total_base = 0; long se_read1_G_total_base = 0; long se_read1_C_total_base = 0; long se_read1_T_total_base = 0; long se_read1_N_total_base = 0; long se_read1_q10_more_80 = 0; long se_read1_q20_more_80 = 0; long se_read1_q30_more_80 = 0; 
		float se_read1_read_float; float se_read1_base_float; float se_read1_GC_rate; float se_read1_q10_rate; float se_read1_q20_rate; float se_read1_q30_rate; float se_read1_A_T_diference_rate; float se_read1_G_C_diference_rate; float se_read1_A_T_diference_base_float; float se_read1_G_C_diference_base_float; float se_read1_A_total_base_rate; float se_read1_G_total_base_rate; float se_read1_C_total_base_rate; float se_read1_T_total_base_rate; float se_read1_N_total_base_rate; float se_read1_q10_more_80_rate; float se_read1_q20_more_80_rate; float se_read1_q30_more_80_rate;
		fp = gzopen(argv[1], "r"); 
		seq = kseq_init(fp);
		while ((l = kseq_read(seq)) >= 0) { 
			char *q = seq->qual.s;
			int c = 0;
			int se_read1_q10_per_read = 0;
			int se_read1_q20_per_read = 0;
			int se_read1_q30_per_read = 0;
			while (c < strlen(seq->qual.s)) {
				if (*q - 33 >= 10) { se_read1_Q10_base++; se_read1_q10_per_read++;}
				if (*q - 33 >= 20) { se_read1_Q20_base++; se_read1_q20_per_read++;}
				if (*q - 33 >= 30) { se_read1_Q30_base++; se_read1_q30_per_read++;}
				q++;
				c++;
			}
			if (se_read1_q10_per_read > c * 0.8){se_read1_q10_more_80++;}
			if (se_read1_q20_per_read > c * 0.8){se_read1_q20_more_80++;}
			if (se_read1_q30_per_read > c * 0.8){se_read1_q30_more_80++;}
			char *s = seq->seq.s;
			int   d = 0;
			while (d < strlen(seq->seq.s)) {
				if (*s == 'A' ){se_read1_A_total_base ++;}
				else if (*s == 'G'){se_read1_gc_base++; se_read1_G_total_base ++;}
				else if (*s == 'C'){se_read1_gc_base++; se_read1_C_total_base ++;}
				else if (*s == 'T'){se_read1_T_total_base ++;}
				else{ se_read1_N_total_base ++;}
				s++;
				d++;
			}
		se_read1_bases += strlen(seq->seq.s);
		se_read1_reads += 1;  
		}
		se_read1_A_T_diference_base_float = abs(se_read1_A_total_base - se_read1_T_total_base);
		se_read1_G_C_diference_base_float = abs(se_read1_G_total_base - se_read1_C_total_base);
		se_read1_read_float = se_read1_reads;
		se_read1_base_float = se_read1_bases;
		se_read1_length = (se_read1_bases/se_read1_reads);
		se_read1_q10_rate = (se_read1_Q10_base / se_read1_base_float) * 100;
		se_read1_q20_rate = (se_read1_Q20_base / se_read1_base_float) * 100;
		se_read1_q30_rate = (se_read1_Q30_base / se_read1_base_float) * 100;
		se_read1_A_T_diference_rate = (se_read1_A_T_diference_base_float / se_read1_base_float) * 100;
		se_read1_G_C_diference_rate = (se_read1_G_C_diference_base_float / se_read1_base_float) * 100;
		se_read1_GC_rate = (se_read1_gc_base / se_read1_base_float) * 100;
		se_read1_A_total_base_rate = (se_read1_A_total_base / se_read1_base_float) * 100;
		se_read1_G_total_base_rate = (se_read1_G_total_base / se_read1_base_float) * 100;
		se_read1_C_total_base_rate = (se_read1_C_total_base / se_read1_base_float) * 100;
		se_read1_T_total_base_rate = (se_read1_T_total_base / se_read1_base_float) * 100;
		se_read1_N_total_base_rate = (se_read1_N_total_base / se_read1_base_float) * 100;
		se_read1_q10_more_80_rate = (se_read1_q10_more_80 / se_read1_read_float) * 100;
		se_read1_q20_more_80_rate = (se_read1_q20_more_80 / se_read1_read_float) * 100;
		se_read1_q30_more_80_rate = (se_read1_q30_more_80 / se_read1_read_float) * 100;
		printf("Sequence Type: SE\n");
		printf("Item\tfq1\n");
		printf("Read length\t%.0d\n",se_read1_length);
		printf("Total number of reads\t%ld\n",se_read1_reads);
		printf("Total number of bases\t%ld\n",se_read1_bases);
		printf("Q10 number of bases\t%ld ( %.3lf )\n",se_read1_Q20_base, se_read1_q10_rate);
		printf("Q20 number of bases\t%ld ( %.3lf )\n",se_read1_Q20_base, se_read1_q20_rate);
		printf("Q30 number of bases\t%ld ( %.3lf )\n",se_read1_Q30_base, se_read1_q30_rate);
		printf("Q10 than 80%% of per read\t%ld ( %.3lf )\n",se_read1_q10_more_80, se_read1_q10_more_80_rate);
		printf("Q20 than 80%% of per read\t%ld ( %.3lf )\n",se_read1_q20_more_80, se_read1_q20_more_80_rate);
		printf("Q30 than 80%% of per read\t%ld ( %.3lf )\n",se_read1_q30_more_80, se_read1_q30_more_80_rate);
		printf("Number of base A\t%ld ( %.3lf )\n",se_read1_A_total_base, se_read1_A_total_base_rate);
		printf("Number of base G\t%ld ( %.3lf )\n",se_read1_G_total_base, se_read1_G_total_base_rate);
		printf("Number of base C\t%ld ( %.3lf )\n",se_read1_C_total_base, se_read1_C_total_base_rate);
		printf("Number of base T\t%ld ( %.3lf )\n",se_read1_T_total_base, se_read1_T_total_base_rate);
		printf("Number of base N\t%ld ( %.3lf )\n",se_read1_N_total_base, se_read1_N_total_base_rate);
		printf("A/T Difference ratio %% : %.3lf\n",se_read1_A_T_diference_rate);
		printf("G/C Difference ratio %% : %.3lf\n",se_read1_G_C_diference_rate);
		printf("GC rate : %.3lf\n",se_read1_GC_rate);
		kseq_destroy(seq); 
		gzclose(fp); 
		return 0;  
	}else{
		int pe_read1_length; int pe_read2_length;
		long pe_read1_reads = 0; long pe_read2_reads = 0; long pe_read1_bases  = 0; long pe_read2_bases  = 0; long pe_read1_Q10_base = 0; long pe_read2_Q10_base = 0; long pe_read1_Q20_base = 0; long pe_read2_Q20_base = 0; long pe_read1_Q30_base = 0; long pe_read2_Q30_base = 0; long pe_read1_gc_base = 0; long pe_read2_gc_base = 0; long pe_read1_A_total_base = 0; long pe_read2_A_total_base = 0; long pe_read1_G_total_base = 0; long pe_read2_G_total_base = 0; long pe_read1_C_total_base = 0; long pe_read2_C_total_base = 0; long pe_read1_T_total_base = 0; long pe_read2_T_total_base = 0; long pe_read1_N_total_base = 0; long pe_read2_N_total_base = 0; long pe_read1_q10_more_80 = 0; long pe_read2_q10_more_80 = 0; long pe_read1_q20_more_80 = 0; long pe_read2_q20_more_80 = 0; long pe_read1_q30_more_80 = 0; long pe_read2_q30_more_80 = 0;
		float pe_read1_read_float; float pe_read2_read_float; float pe_read1_base_float; float pe_read2_base_float; float pe_read1_GC_rate; float pe_read2_GC_rate; float pe_read1_q10_rate; float pe_read2_q10_rate; float pe_read1_q20_rate; float pe_read2_q20_rate; float pe_read1_q30_rate; float pe_read2_q30_rate; float pe_read1_A_T_diference_rate; float pe_read2_A_T_diference_rate; float pe_read1_G_C_diference_rate; float pe_read2_G_C_diference_rate; float pe_read1_A_T_diference_base_float; float pe_read2_A_T_diference_base_float; float pe_read1_G_C_diference_base_float; float pe_read2_G_C_diference_base_float; float pe_read1_A_total_base_rate; float pe_read2_A_total_base_rate; float pe_read1_G_total_base_rate; float pe_read2_G_total_base_rate; float pe_read1_C_total_base_rate; float pe_read2_C_total_base_rate; float pe_read1_T_total_base_rate; float pe_read2_T_total_base_rate; float pe_read1_N_total_base_rate; float pe_read2_N_total_base_rate; float pe_read1_q10_more_80_rate; float pe_read2_q10_more_80_rate; float pe_read1_q20_more_80_rate; float pe_read2_q20_more_80_rate; float pe_read1_q30_more_80_rate; float pe_read2_q30_more_80_rate;
		fp = gzopen(argv[1], "r");
		seq = kseq_init(fp); 
		while ((l = kseq_read(seq)) >= 0) {
			char *q = seq->qual.s;
			int c = 0;
			int pe_read1_q10_per_read = 0;
			int pe_read1_q20_per_read = 0;
			int pe_read1_q30_per_read = 0;
			while (c < strlen(seq->qual.s)) {
				if (*q - 33 >= 10) { pe_read1_Q10_base++; pe_read1_q10_per_read++;}
				if (*q - 33 >= 20) { pe_read1_Q20_base++; pe_read1_q20_per_read++;}
				if (*q - 33 >= 30) { pe_read1_Q30_base++; pe_read1_q30_per_read++;}
				q++;
				c++;
			}
			if (pe_read1_q10_per_read > c * 0.8){pe_read1_q10_more_80++;}
			if (pe_read1_q20_per_read > c * 0.8){pe_read1_q20_more_80++;}
			if (pe_read1_q30_per_read > c * 0.8){pe_read1_q30_more_80++;}
			char *s = seq->seq.s;
			int   d = 0;
			while (d < strlen(seq->seq.s)) {
				if (*s == 'A' ){pe_read1_A_total_base ++;}
				else if (*s == 'G'){pe_read1_gc_base++; pe_read1_G_total_base ++;}
				else if (*s == 'C'){pe_read1_gc_base++; pe_read1_C_total_base ++;}
				else if (*s == 'T'){pe_read1_T_total_base ++;}
				else{ pe_read1_N_total_base ++;}
				s++;
				d++;
			}
		pe_read1_bases += strlen(seq->seq.s);
		pe_read1_reads += 1;
		}
		kseq_destroy(seq);
		gzclose(fp);
		fp = gzopen(argv[2], "r");
		seq = kseq_init(fp);
		while ((l = kseq_read(seq)) >= 0) {
			char *q = seq->qual.s;
			int c = 0;
			int pe_read2_q10_per_read = 0;
			int pe_read2_q20_per_read = 0;
			int pe_read2_q30_per_read = 0;
			while (c < strlen(seq->qual.s)) {
				if (*q - 33 >= 10) { pe_read2_Q10_base++; pe_read2_q10_per_read++;}
				if (*q - 33 >= 20) { pe_read2_Q20_base++; pe_read2_q20_per_read++;}
				if (*q - 33 >= 30) { pe_read2_Q30_base++; pe_read2_q30_per_read++;}
				q++;
				c++;
			}
			if (pe_read2_q10_per_read > c * 0.8){pe_read2_q10_more_80++;}
			if (pe_read2_q20_per_read > c * 0.8){pe_read2_q20_more_80++;}
			if (pe_read2_q30_per_read > c * 0.8){pe_read2_q30_more_80++;}
			char *s = seq->seq.s;
			int   d = 0;
			while (d < strlen(seq->seq.s)) {
				if (*s == 'A' ){pe_read2_A_total_base ++;}
				else if (*s == 'G'){pe_read2_gc_base++; pe_read2_G_total_base ++;}
				else if (*s == 'C'){pe_read2_gc_base++; pe_read2_C_total_base ++;}
				else if (*s == 'T'){pe_read2_T_total_base ++;}
				else{ pe_read2_N_total_base ++;}
				s++;
				d++;
			}
		pe_read2_bases += strlen(seq->seq.s);
		pe_read2_reads += 1;
		}
		pe_read1_A_T_diference_base_float = abs(pe_read1_A_total_base - pe_read1_T_total_base);
		pe_read2_A_T_diference_base_float = abs(pe_read2_A_total_base - pe_read2_T_total_base);
		pe_read1_G_C_diference_base_float = abs(pe_read1_G_total_base - pe_read1_C_total_base);
		pe_read2_G_C_diference_base_float = abs(pe_read2_G_total_base - pe_read2_C_total_base);
		pe_read1_read_float = pe_read1_reads;
		pe_read2_read_float = pe_read2_reads;
		pe_read1_base_float = pe_read1_bases;
		pe_read2_base_float = pe_read2_bases;
		pe_read1_length = (pe_read1_bases/pe_read1_reads);
		pe_read2_length = (pe_read2_bases/pe_read2_reads);
		pe_read1_q10_rate = (pe_read1_Q10_base / pe_read1_base_float) * 100;
		pe_read2_q10_rate = (pe_read2_Q10_base / pe_read2_base_float) * 100;
		pe_read1_q20_rate = (pe_read1_Q20_base / pe_read1_base_float) * 100;
		pe_read2_q20_rate = (pe_read2_Q20_base / pe_read2_base_float) * 100;
		pe_read1_q30_rate = (pe_read1_Q30_base / pe_read1_base_float) * 100;
		pe_read2_q30_rate = (pe_read2_Q30_base / pe_read2_base_float) * 100;
		pe_read1_A_T_diference_rate = (pe_read1_A_T_diference_base_float / pe_read1_base_float) * 100;
		pe_read2_A_T_diference_rate = (pe_read2_A_T_diference_base_float / pe_read2_base_float) * 100;
		pe_read1_G_C_diference_rate = (pe_read1_G_C_diference_base_float / pe_read1_base_float) * 100;
		pe_read2_G_C_diference_rate = (pe_read2_G_C_diference_base_float / pe_read2_base_float) * 100;
		pe_read1_GC_rate = (pe_read1_gc_base / pe_read1_base_float) * 100;
		pe_read2_GC_rate = (pe_read2_gc_base / pe_read2_base_float) * 100;
		pe_read1_A_total_base_rate = (pe_read1_A_total_base / pe_read1_base_float) * 100;
		pe_read2_A_total_base_rate = (pe_read2_A_total_base / pe_read2_base_float) * 100;
		pe_read1_G_total_base_rate = (pe_read1_G_total_base / pe_read1_base_float) * 100;
		pe_read2_G_total_base_rate = (pe_read2_G_total_base / pe_read2_base_float) * 100;
		pe_read1_C_total_base_rate = (pe_read1_C_total_base / pe_read1_base_float) * 100;
		pe_read2_C_total_base_rate = (pe_read2_C_total_base / pe_read2_base_float) * 100;
		pe_read1_T_total_base_rate = (pe_read1_T_total_base / pe_read1_base_float) * 100;
		pe_read2_T_total_base_rate = (pe_read2_T_total_base / pe_read2_base_float) * 100;
		pe_read1_N_total_base_rate = (pe_read1_N_total_base / pe_read1_base_float) * 100;
		pe_read2_N_total_base_rate = (pe_read2_N_total_base / pe_read2_base_float) * 100;
		pe_read1_q10_more_80_rate = (pe_read1_q10_more_80 / pe_read1_read_float) * 100;
		pe_read2_q10_more_80_rate = (pe_read2_q10_more_80 / pe_read2_read_float) * 100;
		pe_read1_q20_more_80_rate = (pe_read1_q20_more_80 / pe_read1_read_float) * 100;
		pe_read2_q20_more_80_rate = (pe_read2_q20_more_80 / pe_read2_read_float) * 100;
		pe_read1_q30_more_80_rate = (pe_read1_q30_more_80 / pe_read1_read_float) * 100;
		pe_read2_q30_more_80_rate = (pe_read2_q30_more_80 / pe_read2_read_float) * 100;
		printf("Sequence Type: PE\n");
		printf("Item\tfq1\tfq2\n");
		printf("Read length\t%.0d\t%.0d\n",pe_read1_length, pe_read2_length);
		printf("Total number of reads\t%ld\t%ld\n",pe_read1_reads, pe_read2_reads);
		printf("Total number of bases\t%ld\t%ld\n",pe_read1_bases, pe_read2_bases);
		printf("Q10 number of bases\t%ld ( %.3lf )\t%ld ( %.3lf )\n",pe_read1_Q20_base, pe_read1_q10_rate, pe_read2_Q20_base, pe_read2_q10_rate);
		printf("Q20 number of bases\t%ld ( %.3lf )\t%ld ( %.3lf )\n",pe_read1_Q20_base, pe_read1_q20_rate, pe_read2_Q20_base, pe_read2_q20_rate);
		printf("Q30 number of bases\t%ld ( %.3lf )\t%ld ( %.3lf )\n",pe_read1_Q30_base, pe_read1_q30_rate, pe_read2_Q30_base, pe_read2_q30_rate);
		printf("Q10 than 80%% of per read\t%ld ( %.3lf )\t%ld ( %.3lf )\n",pe_read1_q10_more_80, pe_read1_q10_more_80_rate, pe_read2_q10_more_80, pe_read2_q10_more_80_rate);
		printf("Q20 than 80%% of per read\t%ld ( %.3lf )\t%ld ( %.3lf )\n",pe_read1_q20_more_80, pe_read1_q20_more_80_rate, pe_read2_q20_more_80, pe_read2_q20_more_80_rate);
		printf("Q30 than 80%% of per read\t%ld ( %.3lf )\t%ld ( %.3lf )\n",pe_read1_q30_more_80, pe_read1_q30_more_80_rate, pe_read2_q30_more_80, pe_read2_q30_more_80_rate);
		printf("Number of base A\t%ld ( %.3lf )\t%ld ( %.3lf )\n",pe_read1_A_total_base, pe_read1_A_total_base_rate, pe_read2_A_total_base, pe_read2_A_total_base_rate);
		printf("Number of base G\t%ld ( %.3lf )\t%ld ( %.3lf )\n",pe_read1_G_total_base, pe_read1_G_total_base_rate, pe_read2_G_total_base, pe_read2_G_total_base_rate);
		printf("Number of base C\t%ld ( %.3lf )\t%ld ( %.3lf )\n",pe_read1_C_total_base, pe_read1_C_total_base_rate, pe_read2_C_total_base, pe_read2_C_total_base_rate);
		printf("Number of base T\t%ld ( %.3lf )\t%ld ( %.3lf )\n",pe_read1_T_total_base, pe_read1_T_total_base_rate, pe_read2_T_total_base, pe_read2_T_total_base_rate);
		printf("Number of base N\t%ld ( %.3lf )\t%ld ( %.3lf )\n",pe_read1_N_total_base, pe_read1_N_total_base_rate, pe_read2_N_total_base, pe_read2_N_total_base_rate);
		printf("A/T Difference ratio %% : %.3lf\t%.3lf\n",pe_read1_A_T_diference_rate, pe_read2_A_T_diference_rate);
		printf("G/C Difference ratio %% : %.3lf\t%.3lf\n",pe_read1_G_C_diference_rate, pe_read2_G_C_diference_rate);
		printf("GC rate : %.3lf\t%.3lf\n",pe_read1_GC_rate, pe_read2_GC_rate);
		kseq_destroy(seq);
		gzclose(fp); 		
		return 0;
	}
}
