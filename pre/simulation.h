//
//  simulation.h
//  Test
//
//  Created by Yudi Zhang on 5/14/20.
//  Copyright © 2020 Yudi Zhang. All rights reserved.
//

#ifndef simulation_h
#define simulation_h

#include <stdio.h>
#include "fastq.h"
#include "nuc.h"
#include "qual.h"

typedef struct _simu_options simu_options;
typedef struct _simu_data simu_dat;

struct _simu_options {
	unsigned int len_N;		/*<! length of sequence sampled from fsa */
	const char *out_file;		/*<! out_file */
	const char *fsa_file;		/*<! fsa files */
	const char *ref_name;		/*<! random sample ref name */
	const char *out_sam;		/*<! alignment of A and B in a sam file */
	const char *extracted_rf;
	char *delim_ref;
	char *delim_len;
	unsigned long seed;
	double homo_rate;		/*<! rate of homologous SNPs rate */
	double heter_rate;		/*<! rate of homologous SNPs rate */
	double alpha;			/*<! beta distribution(proportion of the reference allele at the jth allelic SNP) */
	double beta;			/*<! beta distribution(proportion of the reference allele at the jth allelic SNP) */
	double substitution_rate;	/*<! snps substition rate */
	double prop_allele;		/*<! HWE */
	unsigned int num_ind;
};

struct _simu_data {
	char_t *seq_A;
	char_t *seq_B;
	char_t *seq_A2;
	char_t *seq_B2;
	char_t **ind;
	unsigned int *homo_loci;
	unsigned int *heter_loci;
};

void make_simu_opt(simu_options *opt);
int parse_opt(simu_options *opt, int argc, const char **argv);
int make_data(simu_dat **dat);
int load_data(simu_dat *dat, simu_options *opt, fastq_data *fds);
void fprint_fsa(FILE *fp, char_t **data, size_t n, size_t p, char const * const prefix);
void fprint_seq(FILE *fp, char_t *data, size_t p, char const * const prefix);
void fprint_usage(FILE *fp, const char *cmdname, void *obj);
void write_sam(FILE *fp, char_t *B, size_t p, char *name_A, char *name_B);
#endif /* simulation_h */