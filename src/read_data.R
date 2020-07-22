
#####TODO: MAKE IT READ MULTIPLE LENGTH READS
read_fasta <- function(datafile = NULL)
{
  checkmate::expect_file_exists(datafile, access = "r")
  if ((utils::tail(unlist(strsplit(datafile, "[.]")), 1) != "fasta") && 
      (utils::tail(unlist(strsplit(datafile, "[.]")), 1) != "fsa") &&
      (utils::tail(unlist(strsplit(datafile, "[.]")), 1) != "fa"))
    stop("The input datafile has to be fasta file!")
  
  if (!is.loaded("r_read_fasta", PACKAGE = "sync_data_r")) {
    dyn.load("./src/sync_data_r.so")
  }
  res <- .Call("r_read_fasta", datafile)
  
  names(res) <- c("reads", "dim")
  if (length(res$reads) == res$dim[1] * res$dim[2])
    res$reads <- matrix(res$reads, ncol = res$dim[2], byrow = TRUE)
  
  return(res)
} #read_fasta(datafile = FastaFile)

read_sam <- function(samfile = NULL, ref_name = NULL, 
                     fastq_file = NULL, datafile = NULL) {
  checkmate::expect_file_exists(samfile, access = "r")
  if (utils::tail(unlist(strsplit(samfile, "[.]")), 1) != "sam")
    stop("The input datafile has to be sam file!")
  
  if (!is.loaded("./src/r_read_sam", PACKAGE = "sync_data_r")) {
    dyn.load("./src/sync_data_r.so")
  }
  
  if(is.null(datafile) == TRUE | is.null(samfile) == TRUE)
    stop("Must input a data file and a samfile!")
  
  .Call("r_read_sam", samfile, ref_name, fastq_file, datafile)
}
#read_sam(samfile = samfile, ref_name = ref_name, fastq_file = "sub.fastq", 
#datafile= "sub.txt")

read_fastq <- function(datafile = NULL)
{
  checkmate::expect_file_exists(datafile, access = "r")
  if (utils::tail(unlist(strsplit(datafile, "[.]")), 1) != "fastq")
    stop("The input datafile has to be fastq file!")
  
  if (!is.loaded("r_read_fastq", PACKAGE = "sync_data_r")) {
    dyn.load("./src/sync_data_r.so")
  }
  
  res <- .Call("r_read_fastq", datafile)
  
  names(res) <- c("reads", "quality", "dim")
  res$reads <- matrix(res$reads, ncol = res$dim[2], byrow = TRUE)
  res$quality <- matrix(res$quality, ncol = res$dim[2], byrow = TRUE)

  return(res)
} #read_fastq
#read_fastq(datafile = fastq_file)

call_ampliclust <- function(ampliclust_command = NULL, fastq_file = NULL, ac_outfile = NULL) {
  checkmate::expect_file_exists(fastq_file, access = "r")
  checkmate::expect_file_exists(ampliclust_command, access = "r")
  if (utils::tail(unlist(strsplit(fastq_file, "[.]")), 1) != "fastq" |
      utils::tail(unlist(strsplit(fastq_file, "[.]")), 1) != "fq")
    stop("The input datafile has to be fastq file!")
  if (tail(unlist(strsplit(ampliclust_command, "/")), 1) != "run_ampliclust")
    stop("Ampliclust usage is wrong!")
  if (!dir.exists(sub('/[^/]*$', '', ac_outfile)))
    stop("ac_outfile path does not exist!")
  
  if (!is.loaded("r_ampliclust_init", PACKAGE = "sync_data_r")) {
    dyn.load("./src/sync_data_r.so")
  }
  
  res <- .Call("r_ampliclust_init", ampliclust_command, fastq_file, ac_outfile)
  names(res) <- c("reads", "dim")
  reads <- matrix(res$reads, ncol = res$dim[2], byrow = TRUE)
  return(reads)
}

call_aln <- function(ref_nameA = NULL, ref_nameB = NULL, ref_fsa = NULL, ref_sam = NULL,
                     alnA = NULL, alnB = NULL, out_file = NULL, uni_geno_file = NULL) {
  checkmate::expect_file_exists(ref_fsa, access = "r")
  checkmate::expect_file_exists(ref_sam, access = "r")
  checkmate::expect_file_exists(alnA, access = "r")
  checkmate::expect_file_exists(alnB, access = "r")
  if (utils::tail(unlist(strsplit(ref_fsa, "[.]")), 1) != "fsa" &
      utils::tail(unlist(strsplit(ref_fsa, "[.]")), 1) != "fa" &
      utils::tail(unlist(strsplit(ref_fsa, "[.]")), 1) != "fasta")
    stop("The input ref_fsa has to be fasta file!")
  if (utils::tail(unlist(strsplit(ref_sam, "[.]")), 1) != "sam" &
      utils::tail(unlist(strsplit(alnA, "[.]")), 1) != "sam" &
      utils::tail(unlist(strsplit(alnB, "[.]")), 1) != "sam")
    stop("The input alnA, alnB, ref_sam must be sam file!")
  if (!dir.exists(sub('/[^/]*$', '', out_file)))
    stop("out_file path does not exist!")
  if (!dir.exists(sub('/[^/]*$', '', uni_geno_file)))
    stop("uni_geno_file path does not exist!")
  
  if (!is.loaded("r_make_aln", PACKAGE = "sync_data_r"))
    dyn.load("./src/sync_data_r.so")
  
  .Call("r_make_aln", ref_nameA, ref_nameB, ref_fsa, ref_sam, 
        alnA, alnB, out_file, uni_geno_file)
}
# samfile = "../../data/tetraploid/308-TAN-A.sam"
# ref_name = "Adur420_2:199509_P2"
# read_sam(samfile, ref_name, fastq_file = "res.fastq", datafile = "res.txt")
#haps <- call_ampliclust(ampliclust_command, fastq_file, ac_outfile = "../init")
#ampliclust_command = "../../amplici/run_ampliclust"
#fastq_file = "res.fastq"
# call_aln(ref_nameA = "Genome_A:0-1373", ref_nameB = "Genome_B:0-1373",
#          ref_fsa = "../../data/tpphase/WGS/simu/L_SNP/ref.fsa",
#          ref_sam = "../../data/tpphase/WGS/simu/L_SNP/ref.sam",
#          alnA = "../../data/tpphase/WGS/simu/L_SNP/high_cov/alnA1.sam",
#          alnB = "../../data/tpphase/WGS/simu/L_SNP/high_cov/alnB1.sam",
#          out_file = "./out.txt")



