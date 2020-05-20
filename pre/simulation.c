#include <stdio.h>
#include <curses.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>
#define MATHLIB_STANDALONE 1
#include <Rmath.h>

#include "simulation.h"
#include "cmdline.h"
#include "io.h"
#include "array.h"
#include "myfun.h"

#define PLOIDY 4
int main(int argc, const char *argv[]) {
	int err = NO_ERROR;
	simu_options opt;
	make_simu_opt(&opt);
	if (parse_opt(&opt, argc, argv))
		exit(mmessage(ERROR_MSG, INVALID_CMDLINE, ""));
	
	fastq_data *fds = NULL;
	fastq_options fop = {.read_encoding = XY_ENCODING, .read_names = 1};
	// if use samtools
	//	if (opt.fsa_file) {
	//		char *region = NULL;
	//		unsigned long least, most;
	//		most = rand() % (opt.length + 1); // within 1-opt.length
	//		least = most - opt.len_N;
	//		char temp[strlen(opt.ref_name) + 1];
	//		strcpy(temp, opt.ref_name);
	//		char *ptr = strtok(temp, opt.delim_ref);
	//		size_t length = strlen(ptr) + strlen("-") + strlen(":") + (int)(floor(log10((least))) + 1) + (int)(floor(log10((most))) + 1) + 1;
	//		region = malloc(length);
	//		sprintf(region, "%s%s%zu%s%zu", ptr, ":", least, "-", most);
	//		mmessage(INFO_MSG, NO_ERROR, "sampled region: %s\n", region);
	//		extract_ref(opt.samtools_command, region, opt.fsa_file, opt.extracted_rf);
	//	}
	
	/* read in sampled reference genome */
	if ((err = read_fastq(opt.extracted_rf, &fds, &fop)))
		exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "Reading '%s' "
			      "failed with error '%s' (%d).\n",
			      opt.extracted_rf, fastq_error_message(err),
			      err));
	opt.len_N = fds->n_max_length;
	/* sample another genome with homologous SNPs */
	simu_dat *dat = NULL;
	if (make_data(&dat))
		exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "Make data error \n"));
	
	if (load_data(dat, &opt, fds))
		exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "Load data error \n"));
	
	return(err);
}

void fprint_fsa(FILE *fp, char_t **data, size_t n, size_t p, char const * const prefix) {
	for (size_t i = 0; i < n; ++i) {
		fprintf(fp, ">%s%lu\n", prefix, i);
		for (size_t j = 0; j < p; ++j)
			fprintf(fp, "%c", xy_to_char[(int)data[i][j]]);
		fprintf(fp, "\n");
	}
}

void write_sam(FILE *fp, char_t *B, size_t p, char *name_A, char *name_B) {
	fprintf(fp, "@SQ	SN:");
	fprintf(fp, "%s	LN:%lu\n", name_A, p);
	fprintf(fp, "%s	0	%s	1	255	%luM	*	0	0	\n", name_B, name_A, p);
	for (size_t j = 0; j < p; ++j)
		fprintf(fp, "%c", xy_to_char[(int)B[j]]);
	fprintf(fp, " * \n");
}

void fprint_seq(FILE *fp, char_t *data, size_t p, char const * const prefix) {
	fprintf(fp, ">%s\n", prefix);
	for (size_t j = 0; j < p; ++j)
		fprintf(fp, "%c", xy_to_char[(int)data[j]]);
	fprintf(fp, "\n");
}

void make_simu_opt(simu_options *opt)
{
	//	opt->delim_ref = ":";
	//	opt->delim_len = "-";
	opt->out_sam = NULL;
	opt->extracted_rf = NULL;
	opt->fsa_file = NULL;
	opt->ref_name = NULL;
	opt->len_N = 0;
	opt->out_file = "sample";
	opt->seed = 0;
	opt->alpha = 0.05;
	opt->beta = 0.1;
	opt->heter_rate = 0.001;
	opt->homo_rate = 0.001;
	opt->substitution_rate = 0.3333333;
	opt->num_ind = 1;
} /* make_options */

int make_data(simu_dat **dat) {
	simu_dat *dp;
	*dat = malloc(sizeof(**dat));
	if (*dat == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data object");
	dp = *dat;
	dp->heter_loci = NULL;
	dp->homo_loci = NULL;
	dp->seq_A = NULL;
	dp->seq_B = NULL;
	dp->seq_A2 = NULL;
	dp->seq_B2 = NULL;
	dp->ind = NULL;
	return NO_ERROR;
}

int parse_opt(simu_options *opt, int argc, const char **argv)
{
	int i, j;
	int err = NO_ERROR;
	char a;
	
	for (i = 1; i < argc; ++i) {
		if (strlen(argv[i]) < 2)
			usage_error(argv, i, (void *)opt);
		j = 1;
		a = argv[i][j];
		while (a == '-' && ++j < (int) strlen(argv[i]))
			a = argv[i][j];
		switch(a) {
			case 'f':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				mmessage(INFO_MSG, NO_ERROR, "Fasta file:");
				opt->fsa_file = argv[++i];
				fprintf(stderr, " %s", opt->fsa_file);
				fprintf(stderr, "\n");
				break;
			case 'h':
				fprint_usage(stderr, argv[0], opt);
				exit(EXIT_SUCCESS);
			case 'r':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				mmessage(INFO_MSG, NO_ERROR, "Reference names:");
				opt->out_sam = argv[++i];
				fprintf(stderr, " %s", opt->out_sam);
				fprintf(stderr, "\n");
				break;
			case 's':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->seed = read_uint(argc, argv, ++i, (void *)opt);
				srand(opt->seed);
				mmessage(INFO_MSG, NO_ERROR, "Seed: %lu\n", opt->seed);
				break;
			case 'o':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->out_file = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Out fsa file:");
				fprintf(stderr, " %s", opt->out_file);
				fprintf(stderr, "\n");
				break;
			case 'e':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->extracted_rf = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Fasta file:");
				fprintf(stderr, " %s", opt->extracted_rf);
				fprintf(stderr, "\n");
				break;
			case 'j':
				opt->homo_rate = read_cmdline_double(argc, argv, ++i, (void *)opt);
				break;
			case 'g':
				opt->heter_rate = read_cmdline_double(argc, argv, ++i, (void *)opt);
				break;
			case 'a':
				opt->alpha = read_cmdline_double(argc, argv, ++i, (void *)opt);
				break;
			case 'b':
				opt->beta = read_cmdline_double(argc, argv, ++i, (void *)opt);
				break;
			case 'n':
				opt->num_ind = read_uint(argc, argv, ++i, (void *)opt);
				mmessage(INFO_MSG, NO_ERROR, "No. of individuals:");
				fprintf(stderr, " %d", opt->num_ind);
				fprintf(stderr, "\n");
				break;
			default:
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
		}
	}
	
	return err;
	
CMDLINE_ERROR:
	if (err == NO_ERROR) {
		err = INVALID_CMD_ARGUMENT;
		i--;
	}
	usage_error(argv, i, (void *)opt);
	return err;
} /* parse_options */

int load_data(simu_dat *dat, simu_options *opt, fastq_data *fds) {
	unsigned int i, j, m;
	dat->seq_A = malloc(opt->len_N * sizeof(*dat->seq_A));
	dat->seq_B = malloc(opt->len_N * sizeof(*dat->seq_B));
	dat->seq_A2 = malloc(opt->len_N * sizeof(*dat->seq_A2));
	dat->seq_B2 = malloc(opt->len_N * sizeof(*dat->seq_B2));
	dat->homo_loci = malloc(opt->len_N * sizeof(*dat->homo_loci));
	dat->heter_loci = malloc(opt->len_N * sizeof(*dat->heter_loci));
	MAKE_2ARRAY(dat->ind, 4, opt->len_N);
	
	dat->seq_A = fds->reads;
	// make seq_B
	fprintf(stderr, "start simulation\n");
	for (j = 0; j < opt->len_N; ++j) {
		char_t complement[3];
		unsigned int count = 0;
		double r = rand() / (RAND_MAX + 1.);
		if (r <= opt->homo_rate) {
			dat->homo_loci[j] = 1;
			for(i = 0; i < 4; ++i)
				if (i != dat->seq_A[j])
					complement[count++] = i;
			// sample substitution nuc
			double rn = rand() / (RAND_MAX + 1.);
			//			fprintf(stderr, "%f ",rn);
			if (rn <= opt->substitution_rate)
				dat->seq_B[j] = complement[0];
			else if (rn <= (opt->substitution_rate + opt->substitution_rate))
				dat->seq_B[j] = complement[1];
			else
				dat->seq_B[j] = complement[2];
		}
		else {
			dat->homo_loci[j] = 0;
			dat->seq_B[j] = dat->seq_A[j];
		}
	}
	fprintf(stderr, "\ngenome B simulated, homo rate %lf, mutated:\n", opt->homo_rate);
	for (j = 0; j < opt->len_N; ++j) {
		if(dat->homo_loci[j] == 1)
			fprintf(stderr, "%d: %d|%d \n", j, dat->seq_A[j],dat->seq_B[j]);
	}
	// ref A nad B in the same file for HMM method to use
	FILE *fp = NULL;
	size_t len = strlen(opt->fsa_file) + 5 + 1;
	char *fsa_file = malloc(len);
	sprintf(fsa_file, "%s.fsa", opt->fsa_file);
	if (fsa_file) {
		fp = fopen(fsa_file, "w");
		if (!fp)
			exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, fsa_file));
	}
	size_t name_len = strlen("Genome_A:0-") + (int)log10(opt->len_N) + 1 + 1;
	char *name_A = malloc(name_len);
	char *name_B = malloc(name_len);
	sprintf(name_A, "Genome_A:0-%u", opt->len_N);
	sprintf(name_B, "Genome_B:0-%u", opt->len_N);
	fprint_seq(fp, dat->seq_A, opt->len_N, name_A);
	fprint_seq(fp, dat->seq_B, opt->len_N, name_B);
	fclose(fp);
	fp = NULL;
	if (opt->out_sam) {
		fp = fopen(opt->out_sam, "w");
		if (!fp)
			exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->out_sam));
	}
	// sam file: alignment of A and B
	write_sam(fp, dat->seq_B, opt->len_N, name_A, name_B);
	fclose(fp);
	fp = NULL;
	// ref A nad B in the seperate files for alignment of reads
	sprintf(fsa_file, "%sA.fsa", opt->fsa_file);
	fp = fopen(fsa_file, "w");
	fprint_seq(fp, dat->seq_A, opt->len_N, "Genome_A");
	fclose(fp);
	fp = NULL;
	sprintf(fsa_file, "%sB.fsa", opt->fsa_file);
	fp = fopen(fsa_file, "w");
	fprint_seq(fp, dat->seq_B, opt->len_N, "Genome_B");
	fclose(fp);
	fp = NULL;
	
	// make heter loci
	for (j = 0; j < opt->len_N; ++j) {
		if (dat->homo_loci[j]) {
			dat->heter_loci[j] = 0;
			continue;
		}
		double r = rand() / (RAND_MAX + 1.);
		if (r <= (opt->heter_rate + opt->heter_rate)) {
			double rn = rand() / (RAND_MAX + 1.);
			//			fprintf(stderr, "%f ",rn);
			if (rn <= 0.5)
				dat->heter_loci[j] = 1;
			else
				dat->heter_loci[j] = 2;
		} else {
			dat->heter_loci[j] = 0;
		}
	}
	// make the rest genome A, B
	for (j = 0; j < opt->len_N; ++j) {
		char_t complement[3];
		unsigned int count = 0;
		if(dat->heter_loci[j] == 0) {
			dat->seq_A2[j] = dat->seq_A[j];
			dat->seq_B2[j] = dat->seq_B[j];
		} else if (dat->heter_loci[j] == 1) {
			dat->seq_B2[j] = dat->seq_B[j];
			double rn = rand() / (RAND_MAX + 1.);
			for(i = 0; i < PLOIDY; ++i)
				if (i != dat->seq_A[j])
					complement[count++] = i;
			if (rn <= opt->substitution_rate)
				dat->seq_A2[j] = complement[0];
			else if (rn <= (opt->substitution_rate + opt->substitution_rate))
				dat->seq_A2[j] = complement[1];
			else
				dat->seq_A2[j] = complement[2];
		} else {
			dat->seq_A2[j] = dat->seq_A[j];
			double rn = rand() / (RAND_MAX + 1.);
			for(i = 0; i < PLOIDY; ++i)
				if (i != dat->seq_B[j])
					complement[count++] = i;
			if (rn <= opt->substitution_rate)
				dat->seq_B2[j] = complement[0];
			else if (rn <= (opt->substitution_rate + opt->substitution_rate))
				dat->seq_B2[j] = complement[1];
			else
				dat->seq_B2[j] = complement[2];
		}
	}
	fprintf(stderr, "\nheterzgous sites simulated, heter rate %lf, mutated:\n", opt->heter_rate);
	for (j = 0; j < opt->len_N; ++j) {
		if(dat->heter_loci[j] == 1)
			fprintf(stderr, "A: %d: %d|%d \n", j, dat->seq_A[j], dat->seq_A2[j]);
		else if(dat->heter_loci[j] == 2)
			fprintf(stderr, "B: %d: %d|%d \n", j, dat->seq_B[j], dat->seq_B2[j]);
	}
	set_seed(opt->seed, opt->seed);
	opt->prop_allele = rbeta(opt->alpha, opt->beta);
	double p = opt->prop_allele * opt->prop_allele;
	double pq = 2 * opt->prop_allele * (1 - opt->prop_allele);
	fprintf(stderr, "\nhwe prob:%lf, p^2: %lf 2pq: %lf\n", opt->prop_allele, p, pq);
	
	for (m = 0; m < opt->num_ind; ++m) {
		fprintf(stderr, "Individual %d\n", m);
		for (j = 0; j < opt->len_N; ++j) {
			int hwe = 0; // indicate which genotype
			if (dat->homo_loci[j] == 0 && dat->heter_loci[j] == 0)
				for (i = 0; i < PLOIDY; ++i)
					dat->ind[i][j] = dat->seq_A[j];
			else if (dat->homo_loci[j] == 1 && dat->heter_loci[j] == 0) {
				for (i = 0; i < 2; ++i)
					dat->ind[i][j] = dat->seq_A[j];
				for (i = 2; i < 4; ++i)
					dat->ind[i][j] = dat->seq_B[j];
				fprintf(stderr, "++ %d: %c%c|%c%c\n", j,xy_to_char[dat->ind[0][j]],xy_to_char[dat->ind[1][j]], xy_to_char[dat->ind[2][j]], xy_to_char[dat->ind[3][j]]);
			}
			else if (dat->homo_loci[j] == 0 && dat->heter_loci[j] != 0){
				double rn = rand() / (RAND_MAX + 1.);
				if (rn <= p)
					hwe = 0;
				else if (rn <= p + pq)
					hwe = 1;
				else
					hwe = 2;
				if (dat->heter_loci[j] == 1) {
					for (i = 2; i < 4; ++i)
						dat->ind[i][j] = dat->seq_B[j];
					if (hwe == 0) {
						for (i = 0; i < 2; ++i)
							dat->ind[i][j] = dat->seq_A[j];
					} else if (hwe == 1) {
						fprintf(stderr, "** ");
						dat->ind[0][j] = dat->seq_A[j];
						dat->ind[1][j] = dat->seq_A2[j];
					} else {
						for (i = 0; i < 2; ++i)
							dat->ind[i][j] = dat->seq_A2[j];
					}
					fprintf(stderr, "%d: %c%c|%c%c\n", j,xy_to_char[dat->ind[0][j]],xy_to_char[dat->ind[1][j]], xy_to_char[dat->ind[2][j]], xy_to_char[dat->ind[3][j]]);
				} else {
					for (i = 0; i < 2; ++i)
						dat->ind[i][j] = dat->seq_A[j];
					if (hwe == 0) {
						for (i = 2; i < 4; ++i)
							dat->ind[i][j] = dat->seq_B[j];
					} else if (hwe == 1) {
						fprintf(stderr, "** ");
						dat->ind[2][j] = dat->seq_B[j];
						dat->ind[3][j] = dat->seq_B2[j];
					} else {
						for (i = 2; i < 4; ++i)
							dat->ind[i][j] = dat->seq_B2[j];
					}
					fprintf(stderr, "%d: %c%c|%c%c\n", j,xy_to_char[dat->ind[0][j]],xy_to_char[dat->ind[1][j]], xy_to_char[dat->ind[2][j]], xy_to_char[dat->ind[3][j]]);
				}
			}
		}
		int le = 1;
		if (m != 0)
			le = m;
		
		size_t length = strlen(opt->out_file) + (int)log10(le) + 5 + 1;
		char *file = malloc(length);
		sprintf(file, "%s%u.fsa", opt->out_file, m);
		if (file) {
			fp = fopen(file, "w");
			if (!fp)
				exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, file));
		}
		fprint_fsa(fp, dat->ind, PLOIDY, opt->len_N, "H");
		fclose(fp);
		fp = NULL;
	}
	return NO_ERROR;
}


void fprint_usage(FILE *fp, const char *cmdname, void *obj) {
	simu_options *opt = (simu_options *) obj;
	size_t start = strlen(cmdname) - 1;
	
	while (cmdname[start] != '/' && start) start--;
	if (cmdname[start] == '/') start++;
	
	for (size_t i = start; i < strlen(cmdname); ++i)
		fputc(toupper(cmdname[i]), fp);
	fprintf(fp, "(%d)\n", 1);
	fprintf(fp, "\nNAME\n\t%s - Generate simulated individuals\n", &cmdname[start]);
	fprintf(fp, "\nSYNOPSIS\n\t%s -e <fsa> -f <rf_fsa> -j <homo_rate> -g <heter_rate> -a <alpha> \n -b <beta> -s <seed> -o <output> -n <sample> -r <sam>\t\n", &cmdname[start]);
	fprintf(fp, "\nOPTIONS\n");
	fprintf(fp, "\t-e <fsa> \n\t\tSpecify input fasta file for simulation\n");
	fprintf(fp, "\t-r <sam> \n\t\tSAM file of aligning B genome to A\n");
	fprintf(fp, "\t-f <rf_fsa> \n\t\tSpecify fasta file name for simulated two genomes\n");
	fprintf(fp, "\t-o <outfile> \n\t\tSpecify out file names(Default: %s)\n", opt->out_file);
	fprintf(fp, "\t-j <homo_rate>\n\t\tSpecify homologous rate\n");
	fprintf(fp, "\t-g <heter_rate>\n\t\tSpecify heterzgous rate\n");
	fprintf(fp, "\t-a <alpha>\n\t\tSpecify the shape parameter in beta distribution\n");
	fprintf(fp, "\t-b <beta>\n\t\tSpecify the rate parameter in beta distribution\n");
	fprintf(fp, "\t-s <seed>\n\t\tRandom number generator seed\n");
	fprintf(fp, "\t-n <sample>\n\t\tNo. of individuals simulated\n");
	fprintf(fp, "\t-h \n\t\tThis help\n");
} /* fprint_usage */
