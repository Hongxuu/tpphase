#include <stdio.h>
#include <curses.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>

#include "sam.h"
#include "fastq.h"
#include "nuc.h"
#include "qual.h"
#include "uthash.h"
#include "io.h"
#include "array.h"
#include "pick_reads.h"
#include "options.h"
#include "myfun.h"

int main(int argc, const char *argv[])
{
	options opt;
	int debug_level = QUIET;//ABSOLUTE_SILENCE;//MINIMAL;//DEBUG_I;//
	int err = NO_ERROR;
	sam *sds[N_FILES] = {NULL, NULL};
	fastq_options fop = {.read_encoding = IUPAC_ENCODING, .read_names = 1};
	sam_hash *by_name[N_FILES] = {NULL, NULL};//, NULL, NULL};
	size_t my_refs[N_FILES] = {0, 0};//, 0, 0};
	FILE *fp = NULL;
	unsigned int j, i;
	default_options(&opt);
	if (parse_options(&opt, argc, argv))
		exit(mmessage(ERROR_MSG, INVALID_CMDLINE, ""));
	options_rf opt_rf;
	make_options(&opt_rf);
	opt_rf.sam_file = opt.sam_file;
	ref_info *rf_info = NULL;
	fastq_data *fdr = NULL;
	
	if ((err = read_fastq(opt.uni_genome, &fdr, &fop)))
		exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "Reading '%s' "
			      "failed with error '%s' (%d).\n",
			      opt.uni_genome, fastq_error_message(err),
			      err));
	// get the selected reference (match the name)
	unsigned int A_id = 0, B_id = 0;
	unsigned int rptr = 0, rptr_b = 0;
	char *A_name = NULL;
	char *names = NULL;
//	printf("%lu\n", strlen(opt.ref_names[0]));
//	strlen has to be the same when comparing
	fdr->n_lengths = malloc(fdr->n_reads * sizeof(*fdr->n_lengths));
	if (fdr->n_max_length == fdr->n_min_length)
		for (j = 0; j < fdr->n_reads; ++j)
			fdr->n_lengths[j] = fdr->n_min_length;
		
	for (j = 0; j < fdr->n_reads; ++j) {
		names = fdr->names;
		for (i = 0; i < j; ++i)
			names += fdr->name_lengths[i];
		A_name = malloc((fdr->name_lengths[j] + 1) * sizeof(*A_name));
//		for (i = 0; i < fdr->name_lengths[j]; ++i)
//			A_name[i] = names[i];
		strncpy (A_name, names, fdr->name_lengths[j]);
		A_name[fdr->name_lengths[j]] = '\0';
//		printf("%lu %s\n", strlen(A_name), A_name);
		if (!strcmp(opt.ref_names[0], A_name)) {
			A_id = j;
			B_id = j+1;
			rptr_b = rptr + fdr->n_lengths[j];
			break;
		}
		rptr += fdr->n_lengths[j];
	}
	
	// make universal genome, gap from A as I(4), gap from B as J(5), mismatch mas M (6)
	data_t to_xy[NUM_IUPAC_SYMBOLS] = {
		0, XY_A, XY_C, 0, XY_G, 0, 0, 0, XY_T, 0, 0, 0, 0, 0, 0, 0
	};
	data_t *uni_genome = malloc(fdr->n_lengths[A_id] * sizeof(*uni_genome));
	
	for (j = 0; j < fdr->n_lengths[A_id]; ++j) {
		if (fdr->reads[rptr + j] == fdr->reads[rptr_b + j]) {
			uni_genome[j] = to_xy[fdr->reads[rptr + j]];
		} else if (fdr->reads[rptr + j] != 0 && fdr->reads[rptr_b + j] != 0 &&
			   fdr->reads[rptr + j] != fdr->reads[rptr_b + j]) {
			uni_genome[j] = 6;
		} else if (fdr->reads[rptr + j] != 0 && fdr->reads[rptr_b + j] == 0) {
			uni_genome[j] = 5;
		} else {
			uni_genome[j] = 4;
		}
	}
//	PRINT_VECTOR(uni_genome, fdr->n_lengths[A_id]);
	// store the index after alignment, if gap use -1, this alignment does not contain softclips, when map the reads back, the start and position should consider its length
	unsigned int gap_a = 0;
	unsigned int gap_b = 0;
	int *id_A = malloc(fdr->n_lengths[A_id] * sizeof(*id_A));
	int *id_B = malloc(fdr->n_lengths[A_id] * sizeof(*id_B));
	for (j = 0; j < fdr->n_lengths[A_id]; ++j) {
		id_A[j] = j - gap_a;
		id_B[j] = j - gap_b;
		if (uni_genome[j] == 4) {
			gap_a++;
			id_A[j] = -1;
		} else if (uni_genome[j] == 5) {
			gap_b++;
			id_B[j] = -1;
		}
	}
//	PRINT_VECTOR(id_A, fdr->n_lengths[A_id]);
//	printf("B\n");
//	PRINT_VECTOR(id_B, fdr->n_lengths[A_id]);
	/* store information to the reference targeted sam file */
	make_targets_info(opt_rf, &rf_info);
	
	// read in the alignments to A B
	for (j = 0; j < N_FILES; ++j) {
		fp = fopen(opt.sbam_files[j], "r");
		if (!fp)
			exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR,
				      opt.sbam_files[j]));
		read_sam(fp, &sds[j]);
		fclose(fp);
	}
	fp = NULL;
	// pick the reads, more reads could be picked
	pickreads(rf_info, &opt_rf, sds);
	
	char *strand;
	for (j = 0; j < N_FILES; ++j) {
		printf("Genome %d\n", j);
		/* find selected references index in sam files */
		size_t rchar = 0;
		unsigned char found = 0;
		for (size_t i = 0; i < sds[j]->n_se; ++i) {
			sam_entry *se = &sds[j]->se[i];
			
			/* skip unmapped */
			if ((se->flag & (1 << 2)))
				continue;
			if (se->which_ref == -1)
				continue;
			
			if (!strcmp(opt.ref_names[j], se->ref_name)) {
				se->name_s = NULL;
				/* strand for hashing on strand and name */
				if ((se->flag & 16) == 0) {
					strand = "+";
					// flip the strand if A is aligned to reverse complement of B
					if (j == 1 && rf_info->info[se->which_ref].strand_B == 1)
						strand = "-";
				} else {
					strand = "-";
					if (j == 1 && rf_info->info[se->which_ref].strand_B == 1)
						strand = "+";
				}
				
				size_t length = strlen(se->name) + strlen(strand) + 1;
				se->name_s = malloc(length);
				sprintf(se->name_s, "%s%s", se->name, strand);
				my_refs[j] = se->which_ref; // this my_refs index should be adjusted
				found = 1;
				printf("%s\t", se->name_s);
			}
			rchar += strlen(se->ref_name) + 1;
		}
		printf("\n");
		if (!found)
			exit(mmessage(ERROR_MSG, INVALID_USER_INPUT, "no "
				      "reference '%s' in fasta file '%s'",
				      opt.ref_names[j], opt.sbam_files[j]));
		
		/* hash sam file to reference (use n_se since some references are repeated in the targted sam file) */
		hash_sam(sds[j], &by_name[j], HASH_REFERENCE, my_refs[j], rf_info->ref_sam->n_se,
			 opt.drop_unmapped, opt.drop_secondary,
			 opt.drop_soft_clipped, opt.drop_indel,
			 opt.min_length, opt.max_length, opt.max_eerr);
		
		mmessage(INFO_MSG, NO_ERROR, "Number of %u alignments: %zu\n",
			 j, sds[j]->n_per_ref[my_refs[j]]);
	}
	
	//merge read pair, if read only align to one genome, keep the alignment information, but do not pick the best alignment
	merge_hash *mh = NULL;
	size_t total_read = hash_merge(&mh, N_FILES, sds, my_refs);
	printf("total picked reads %lu\n", total_read);

	// oh well, this is 0-based, so plus 1 to match with what the rest related code designed for
	sam_entry *fse = &rf_info->ref_sam->se[my_refs[0]];
	ref_entry *re = &rf_info->info[my_refs[0]];
	re->start_B++;
//	re->end_B;
	re->start_A++;
//	re->end_A;
	// store aligned index used in the real genome
	
	long *real_id_A = malloc(fdr->n_lengths[A_id] * sizeof(*real_id_A));
	long *real_id_B = malloc(fdr->n_lengths[A_id] * sizeof(*real_id_B));
	
	real_id_A[0] = re->start_A + fse->pos - 1; // A start index, 1 based
	if (re->strand_B) { // if B reversed
		printf("Genome B is reverse complemented\n");
		real_id_B[0] = re->end_B;
		if (fse->cig->ashes[fse->cig->n_ashes - 1].type == CIGAR_SOFT_CLIP || fse->cig->ashes[fse->cig->n_ashes - 1].type == CIGAR_HARD_CLIP)
			real_id_B[0] -= fse->cig->ashes[fse->cig->n_ashes - 1].len;
	} else {
		real_id_B[0] = re->start_B;
		if (fse->cig->ashes[0].type == CIGAR_SOFT_CLIP || fse->cig->ashes[0].type == CIGAR_HARD_CLIP)
			real_id_B[0] += fse->cig->ashes[0].len; // length of unaligned in B
	}
	
	for (j = 1; j < fdr->n_lengths[A_id]; ++j) {
		if (id_A[j] == -1) {
			real_id_A[j] = -1;
			if (re->strand_B) // if B is reversed
				real_id_B[j] = real_id_B[0] - id_B[j];
			else
				real_id_B[j] = real_id_B[0] + id_B[j];
			
		} else if (id_B[j] == -1) {
			real_id_B[j] = -1;
			real_id_A[j] = real_id_A[0] + id_A[j];
		} else {
			real_id_A[j] = real_id_A[0] + id_A[j];
			if (re->strand_B) // if B is reversed
				real_id_B[j] = real_id_B[0] - id_B[j];
			else
				real_id_B[j] = real_id_B[0] + id_B[j];
		}
	}
	for (j = 1; j < fdr->n_lengths[A_id]; ++j)
		printf("%ld: %c\t\t\t", real_id_A[j], iupac_to_char[fdr->reads[rptr + j]]);
	for (j = 1; j < fdr->n_lengths[A_id]; ++j)
		printf("%ld: %c\t", real_id_B[j], iupac_to_char[fdr->reads[rptr_b + j]]);
	printf("\nreal start-end in genome A %ld-%ld B %ld-%ld\n", real_id_A[0], real_id_A[fdr->n_lengths[A_id] - 1], real_id_B[0], real_id_B[fdr->n_lengths[A_id] - 1]);
	PRINT_VECTOR(real_id_A, fdr->n_lengths[A_id]);
	printf("B\n");
	PRINT_VECTOR(real_id_B, fdr->n_lengths[A_id]);
	
	unsigned int strand_genome;
	int *id_uni = NULL;
//	size_t uni_len;
	long *real_id = NULL;
	unsigned int rf_id = 0;
	
//	for (i = 0; i < fdr->n_lengths[A_id]; ++i) {
//		printf("%d:%d:%c:%ld\t", i, id_A[i], iupac_to_char[fdr->reads[i + rptr]], real_id_A[i]);
//	}
//	printf("\n");
//	for (i = 0; i < fdr->n_lengths[B_id]; ++i) {
//		printf("%d:%d:%c:%ld\t", i, id_B[i], iupac_to_char[fdr->reads[i + rptr_b]], real_id_B[i]);
//	}
//	printf("\n");
	printf("start!\n");
	if (opt.out_file) {
		fp = fopen(opt.out_file, "w");
		if (!fp)
			exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt.out_file));
	}
	unsigned int n_read = 1;
	for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
		sam_entry *se;
		
		// only mapped to one reference, adjust the alignment to universal
		if (me->nfiles != N_FILES) {
			me->exclude = 1;
			if (me->indices[0]) {
				se = &sds[0]->se[me->indices[0][0]]; // indices[0] represents A
				strand_genome = re->strand_A;
				id_uni = id_A;
				rf_id = rptr;
//				uni_len = fdr->n_lengths[A_id] - gap_a;
				real_id = real_id_A;
			} else {
				se = &sds[1]->se[me->indices[1][0]];
				strand_genome = re->strand_B;
				id_uni = id_B;
				rf_id = rptr_b;
//				uni_len = fdr->n_lengths[A_id] - gap_b;
				real_id = real_id_B;
			}
			printf("%s is ajusted\ndata\n", se->name_s);
			adjust_alignment(se, &fdr->reads[rf_id], strand_genome, id_uni, fdr->n_lengths[A_id], real_id);
			// output the final alignment
			output_data(fp, se, n_read);
			n_read++;
			continue;
		}
		
		/* force one alignment per sub-genome */
		double max_ll = -INFINITY;
		unsigned int max_id = 0;
		for (j = 0; j < N_FILES; ++j) {
//			printf("genome %d\n", j);
			if (me->count[j] > 1)
				exit(mmessage(ERROR_MSG, INTERNAL_ERROR,
					      "Read %u aligns twice in genome %s.\n",
					      j, sds[j]->se[me->indices[j][0]].name_s));
			
			se = &sds[j]->se[me->indices[j][0]];
			if (j == 0) {
				strand_genome = re->strand_A;
				id_uni = id_A;
				rf_id = rptr;
//				uni_len = fdr->n_lengths[A_id] - gap_a;
				real_id = real_id_A;
			} else {
				strand_genome = re->strand_B;
				id_uni = id_B;
				rf_id = rptr_b;
//				uni_len = fdr->n_lengths[A_id] - gap_b;
				real_id = real_id_B;
			}
			printf("%s is ajusted\n", se->name_s);
			adjust_alignment(se, &fdr->reads[rf_id], strand_genome, id_uni, fdr->n_lengths[A_id], real_id);
			if (max_ll < se->ll_aln) {
				max_ll = se->ll_aln;
				max_id = j;
			}
		}
		printf("choose %d\n", max_id);
		se = &sds[max_id]->se[me->indices[max_id][0]];
		// output the needed format
		output_data(fp, se, n_read);
		n_read++;
	}
	printf("read: %d\n", n_read);
	fclose(fp);
}
