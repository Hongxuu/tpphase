
#####TODO: MAKE IT READ MULTIPLE LENGTH READS
read_fasta <- function(datafile = NULL)
{
  checkmate::expect_file_exists(datafile, access = "r")
  if ((utils::tail(unlist(strsplit(datafile, "[.]")), 1) != "fasta") & 
      (utils::tail(unlist(strsplit(datafile, "[.]")), 1) != "fa"))
    stop("The input datafile has to be fasta file!")
  
  if (!is.loaded("r_read_fasta", PACKAGE = "sync_data_r")) {
    dyn.load("sync_data_r.so")
  }
  
  res <- .Call("r_read_fasta", datafile)
  
  names(res) <- c("reads", "dim")
  reads <- matrix(res$reads, ncol = res$dim[2], byrow = TRUE)
  
  return(reads)
} #read_fasta(datafile = FastaFile)

read_sam <- function(samfile = NULL, ref_name = NULL, 
                     fastq_file = NULL, datafile = NULL) {
  checkmate::expect_file_exists(samfile, access = "r")
  if (utils::tail(unlist(strsplit(samfile, "[.]")), 1) != "sam")
    stop("The input datafile has to be sam file!")
  
  if (!is.loaded("r_read_sam", PACKAGE = "sync_data_r")) {
    dyn.load("sync_data_r.so")
  }
  
  if(is.null(datafile) == TRUE | is.null(samfile) == TRUE | is.null(fastq_file) == TRUE | is.null(ref_name) == TRUE)
    stop("Must input data file and a reference name!")
  
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
    dyn.load("sync_data_r.so")
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
  if (utils::tail(unlist(strsplit(fastq_file, "[.]")), 1) != "fastq")
    stop("The input datafile has to be fastq file!")
  if (tail(unlist(strsplit(ampliclust_command, "/")), 1) != "run_ampliclust")
    stop("Ampliclust usage is wrong!")
  if (!dir.exists(sub('/[^/]*$', '', ac_outfile)))
    stop("ac_outfile path does not exist!")
  
  if (!is.loaded("r_ampliclust_init", PACKAGE = "sync_data_r")) {
    dyn.load("sync_data_r.so")
  }
  
  res <- .Call("r_ampliclust_init", ampliclust_command, fastq_file, ac_outfile)
  names(res) <- c("reads", "dim")
  reads <- matrix(res$reads, ncol = res$dim[2], byrow = TRUE)
  return(reads)
}
# samfile = "../../data/tetraploid/308-TAN-A.sam"
# ref_name = "Adur420_2:199509_P2"
# read_sam(samfile, ref_name, fastq_file = "res.fastq", datafile = "res.txt")
#haps <- call_ampliclust(ampliclust_command, fastq_file, ac_outfile = "../init")
#ampliclust_command = "../../amplici/run_ampliclust"
#fastq_file = "res.fastq"



