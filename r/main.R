
datafile = "../../data/hmm/WGS/simu/ref3/low_cov/out.txt"
alignment = "../../data/hmm/WGS/simu/ref3/ref.fsa"
call_aln(ref_nameA = "Genome_A:0-2000", ref_nameB = "Genome_B:0-2000",
         ref_fsa = alignment,
         ref_sam = "../../data/hmm/WGS/simu/ref3/ref.sam",
         alnA = "../../data/hmm/WGS/simu/ref3/low_cov/aln0A.sam",
         alnB = "../../data/hmm/WGS/simu/ref3/low_cov/aln0B.sam",
         out_file = datafile)
ref1.low <- altragenotype(datafile = "../../data/hmm/WGS/simu/ref1/L_SNP/low_cov/out.txt", 
              alignment = "../../data/hmm/WGS/simu/ref1/L_SNP/ref.fsa")
ref1.low.snp <- confusion_matric(res = ref1.low, truth_file = "../../data/hmm/WGS/simu/ref1/L_SNP/sim0.fsa")

ref1.high <- altragenotype(datafile = "../../data/hmm/WGS/simu/ref1/L_SNP/high_cov/out.txt", 
                          alignment = "../../data/hmm/WGS/simu/ref1/L_SNP/ref.fsa")
ref1.high.snp <- confusion_matric(res = ref1.high, truth_file = "../../data/hmm/WGS/simu/ref1/L_SNP/sim0.fsa")
ref2.low <- altragenotype(datafile = "../../data/hmm/WGS/simu/ref2/LOW_SNP/low_cov/out.txt", 
                          alignment = "../../data/hmm/WGS/simu/ref2/LOW_SNP/ref.fsa")
ref2.low.snp <- confusion_matric(res = ref2.low, truth_file = "../../data/hmm/WGS/simu/ref2/LOW_SNP/indiv0.fsa")
ref2.high <- altragenotype(datafile = "../../data/hmm/WGS/simu/ref2/LOW_SNP/mid_cov/out.txt", 
                          alignment = "../../data/hmm/WGS/simu/ref2/LOW_SNP/ref.fsa")
ref2.high.snp <- confusion_matric(res = ref2.high, truth_file = "../../data/hmm/WGS/simu/ref2/LOW_SNP/indiv0.fsa")
## mnlogit took a relative same time when data are more, otherwise, it takes shorter (1/2 time as bw)

ref3.low <- altragenotype(datafile = "../../data/hmm/WGS/simu/ref3/low_cov/out.txt", 
                           alignment = "../../data/hmm/WGS/simu/ref3/ref.fsa")
ref3.low.snp <- confusion_matric(res = ref3.low, truth_file = "../../data/hmm/WGS/simu/ref3/indiv0.fsa")





