
rm(list=ls())
suppressMessages(library("optparse"))
suppressMessages(library("bigreadr"))
suppressMessages(library("stringr"))
suppressMessages(library("caret"))
suppressMessages(library("Rcpp"))
suppressMessages(library("RcppArmadillo"))
suppressMessages(library("inline"))
suppressMessages(library("glmnet"))
suppressMessages(library("MASS"))
suppressMessages(library("doMC"))
suppressMessages(library("foreach"))

# ## progress bar function from https://stackoverflow.com/questions/51213293/is-it-possible-to-get-a-progress-bar-with-foreach-and-a-multicore-kind-of-back
progBar <- function(ii, N, per = 10) {
  #ii is current iteration.
  #N is total number of iterations to perform.
  #per is step (percent) when to update the progress bar. We need this as when multiple iterations are being performed progress bar gets updated too often and output is messed up.
  if (ii %in% seq(1, N, per)) {
    x <- round(ii * 100 / N)
    message("[ ",
            paste(rep("=", x), collapse = ""),
            paste(rep("-", 100 - x), collapse = ""),
            " ] ", x, "%", "\r",
            appendLF = FALSE)
    if (ii == N) cat("\r")
  }
}


option_list = list(
  make_option("--PATH_package", action="store", default=NA, type='character',
              help="Full path to the directory for the downloaded file after decompression [required]"),
  make_option("--PATH_out", action="store", default=NA, type='character',
              help="Output directory of results [required]"),
  make_option("--PATH_plink", action="store", default=NA, type='character',
              help="Path to plink2 executable [required]"),
  make_option("--FILE_sst", action="store", default=NA, type='character',
              help="Full path and the file name of the GWAS summary statistics, separated by comma [required] [must have columns: rsid, chr, beta, beta_se, a1, a0, n_eff]"),
  make_option("--pop", action="store", default=NA, type='character',
              help="Population of the GWAS sample, separated by comma [required]"),
  make_option("--num_samples", action="store", default=NA, type='character',
              help="Population-wise number of samples (comma-separated values), first for target"),
  make_option("--chrom", action="store", default="1-22", type='character',
              help="The chromosome on which the model is fitted, separated by comma or dash for consecutive chromosomes [required]"),
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]"),
  make_option("--NCORES", action="store", default=1, type="integer",
              help="How many cores to use [default: %default]"),
  make_option("--bfile_testing", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus bed/bim/fam) for testing [required]"),
  make_option("--pheno_testing", action="store", default=NA, type='character',
              help="Path to phenotype file for testing (PLINK format) [optional, taken from bfile otherwise]"),
  make_option("--covar_testing", action="store", default=NA, type='character',
              help="Path to quantitative covariates for testing (PLINK format) [optional]"),
  make_option("--testing", action="store", default=F, type='logical',
              help="Whether to perform testing in seperate dataset [required]"),
  make_option("--cleanup", action="store", default=T, type="logical",
              help="Cleanup temporary files or not [default: %default]")
  
)
opt = parse_args(OptionParser(option_list=option_list))
num_samples <- as.numeric(unlist(strsplit(opt$num_samples, ",")))

NCORES <- opt$NCORES

ethnic = str_split(opt$pop,",")[[1]]; M <- length(ethnic)
sumdata_path = str_split(opt$FILE_sst,",")[[1]]


# Perform i/o checks here:
for ( f in sumdata_path ) {
  if ( !file.exists(f) ){
    cat( "ERROR: ", f , " input file does not exist\n" , sep='', file=stderr() )
    q()
  }
}

opt$chrom <- gsub("-",":",opt$chrom)
eval(parse(text=paste0("allchrom = c(",opt$chrom,")")))

suppressWarnings(dir.create(opt$PATH_out))

if(! dir.exists(opt$PATH_out)){
  cat( "ERROR: output path does not exist\n" )
  q()
}

suppressWarnings(dir.create(paste0(opt$PATH_out, "/tmp")))
suppressWarnings(dir.create(paste0(opt$PATH_out, "/tmp/PRS_in_all_settings_bychrom")))
suppressWarnings(dir.create(paste0(opt$PATH_out, "/before_ensemble")))
suppressWarnings(dir.create(paste0(opt$PATH_out,"/tmp/sample_scores_",ethnic[1])))

sourceCpp("/dcs04/nilanjan/data/Anagh/PRS_proj/code_github/Promit_projects/grad_func_PRS.cpp")


########################################################################
########################################################################

if ( opt$verbose >= 1 ) cat("\n** Step 1. Preprocessing data **\n")

# Initialize a data frame to store R^2 values
results <- data.frame(iter = integer(), eta = numeric(), alpha = numeric(), R2 = numeric())

# Define parameter grids (non-ADAM)
iters <- c(10, 100)
etas <- c(0.1, 0.01, 0.001, 0.0001)  # Use same eta for all
alphas <- c(0.1, 0.01, 0.001, 0)  # Use same alpha for all

# Define parameter grids (non-ADAM)
# iters <- c(1:10)
# etas <- c(0.0001)  # Use same eta for all
# alphas <- c(0)  # Use same alpha for all

for (iter in iters) {
  for (eta in etas) {
    for (alpha in alphas) {

############
## Step 1.1. Loading summary data and matching to reference data

ref <- fread2(paste0(opt$PATH_package,"/ref_bim.txt"))

df_beta_list <- vector("list", length = M)
summ_max <- numeric(length = M)
N0 <- integer(length = M)
for (l in 1:M){
  
  # Match summary statistics to reference data
  df_beta <- fread2(sumdata_path[l])
  df_beta <- df_beta[!is.na(df_beta$beta) & !is.na(df_beta$beta_se), ]
  df_beta$n_eff[is.na(df_beta$n_eff)] <- mean(df_beta$n_eff, na.rm=T)
  df_beta <- df_beta[df_beta$rsid %in% ref$V2,]
  ref_tmp <- ref[match(df_beta$rsid, ref$V2),]
  
  # Match effect alleles to that in reference data
  tmp1 <- paste0(df_beta$a1, ":", df_beta$a0); tmp2 <- paste0(df_beta$a0, ":", df_beta$a1)
  tmp0 <- paste0(ref_tmp$V5, ":", ref_tmp$V6)
  flip <-  tmp2 == tmp0
  keep <- (tmp1 == tmp0 | tmp2 == tmp0)
  if ( opt$verbose == 2 ){
    messageflip <- paste0(". ",sum(flip)," SNPs are flipped and corrected \n")
  }else{
    messageflip <- paste0(" \n")
  }
  
  # Adjust direction of effect for flipped variants
  if(sum(flip)>0){
    df_beta$beta[flip] <- -1 * df_beta$beta[flip]
  }
  df_beta <- df_beta[keep,,drop=F]
  
  
  # Get scale factor for standardized effect size
  N0[l] <- max(df_beta$n_eff, na.rm=T)
  df_beta$n_eff[is.na(df_beta$n_eff)] <- N0[l]
  df_beta$snps_scale <- sqrt(df_beta$n_eff * df_beta$beta_se^2 + df_beta$beta^2)
  df_beta$beta_hat <- df_beta$beta / df_beta$snps_scale
  df_beta_list[[l]] <- df_beta
  
  rm(list = c("df_beta","ref_tmp","tmp0","tmp1","tmp2","flip","keep"))
}

########################################################################
########################################################################

if ( opt$verbose >= 1 ) cat("\n** Step 2. Fitting models by chromosome **\n")

registerDoMC(NCORES)

# Run algorithm parallelled by chromosomes

ff <- foreach(j = 1:length(allchrom), ii = icount(), .final = function(x) NULL) %dopar% {
  
  chr <- allchrom[j]
  
  ############
  ## Step 2.1. Extract variants in both provided summary statistics and reference data
  if (opt$verbose == 2) cat("\n** Step 2.1 started for chromosome ", chr, " **\n")
  Nsnps0 <- vector("list", length = M)
  snps_list0 <- vector("list", length = M)
  LD_list0 <- vector("list", length = M)
  summ_list0 <- vector("list", length = M)
  snps_scale0 <- vector("list", length = M)
  for (l in 1:M){
    
    load(paste0(opt$PATH_package,"/",ethnic[l],"/chr",chr,"_LD.RData"))
    df_beta <- df_beta_list[[l]]
    
    # Mark overlapped variants
    m <- lapply(snps_list, FUN=function (x){x %in% df_beta$rsid})
    tmpLD <- LD_list
    tmpSNP <- snps_list
    for(i in 1:length(m)){
      if(length(tmpLD[[i]])>0){
        # Subset reference data by overlapped variants
        tmpLD[[i]] <- tmpLD[[i]][m[[i]],m[[i]],drop=F]; tmpLD[[i]][is.nan(tmpLD[[i]])] <- 1
        tmpSNP[[i]] <- tmpSNP[[i]][m[[i]]]
        
        # Remove variants due to homozygosity or perfect correlations
        if(nrow(tmpLD[[i]])>1){ drop = findCorrelation(tmpLD[[i]],cutoff = 0.99999) }else{ drop <- integer(0) }
        if(length(drop)>0){
          tmpLD[[i]] <- tmpLD[[i]][-drop, -drop, drop=F]
          tmpSNP[[i]] <- tmpSNP[[i]][-drop]
        }
      }
    }
    LD_list0[[l]] <- tmpLD
    snps_list0[[l]] <- tmpSNP
    Nsnps0[[l]] <- unlist(lapply(tmpSNP, length))
    
    # Match standardized effect size and scale factors
    tmp <- lapply(snps_list0[[l]], FUN=function (x){ df_beta[match(x, df_beta$rsid),] } )
    summ_list0[[l]] <- lapply(tmp, FUN=function (x){ x$beta_hat } )
    snps_scale0[[l]] <- lapply(tmp, FUN=function (x){ x$snps_scale } )
    
    if ( opt$verbose == 2 ) cat(paste0(sum(Nsnps0[[l]])," SNPs are included in the analysis of ",ethnic[l]," on CHR",chr,messageflip))
    
    rm(list = c("i","LD_list","Nsnps","snps_list","tmp","tmpLD","tmpSNP","m","df_beta"))
  }
  
  nblock <- length(LD_list0[[l]])
  if (opt$verbose == 2) cat("\n** Step 2.1 ended for chromosome ", chr, " **\n")
  ############
  ## Step 2.2. Transform to standard data format
  if (opt$verbose == 2) cat("\n** Step 2.2 started for chromosome ", chr, " **\n")
  # Organize data in a structure fit in the algorithm
  indx_block1 <- integer(length = nblock)
  snp_list1 <- vector("list", length = nblock)
  Nsnps1 <- integer(length = nblock)
  indx1 <- vector("list", length = nblock)
  summ_list1 <- vector("list", length = nblock)
  snps_scale1 <- vector("list", length = nblock)
  LD_list1 <- vector("list", length = nblock)

  for (bl in 1:nblock){

    ## snp_list1, Nsnps1
    snp_list_tmp <- vector("list", length = M)
    tmp <- character()
    for (l in 1:M){
      if(is.null(snps_list0[[l]][[bl]])) { next }
      snp_list_tmp[[l]] <- snps_list0[[l]][[bl]]
      tmp <- c(tmp, snp_list_tmp[[l]])
    }
    tmp <- unique(tmp)
    Nsnps1[bl] <- length(tmp)
    if(Nsnps1[bl]==0){ indx_block1[bl] <- 0; next }
    snp_list1[[bl]] <- tmp
    indx_block1[bl] <- 1

    ## indx1: the position of ref SNP in original summ_list
    ## summ_list1: summ stat matched to reference snp list (set to 0 for snps not in a certain ethnic group)
    ## LD_list1: LD correlations matched to reference snp list (set to 0 for snps not in a certain ethnic group)
    indx_tmp <- matrix(nrow = Nsnps1[bl], ncol = M)
    summ_list_tmp <- matrix(nrow = Nsnps1[bl], ncol = M)
    snps_scale_tmp <- matrix(nrow = Nsnps1[bl], ncol = M)
    LD_list_tmp <- vector("list", length = M)
    for (l in 1:M){
      m <- match(snp_list1[[bl]], snp_list_tmp[[l]]); m1 <- m; m1[is.na(m1)] <- 0; indx_tmp[,l] <- m1
      m1 <- summ_list0[[l]][[bl]][m]; summ_list_tmp[,l] <- m1
      m1 <- snps_scale0[[l]][[bl]][m]; snps_scale_tmp[,l] <- m1
      m1 <- as.matrix(LD_list0[[l]][[bl]][m,m]); m1[is.na(m1)] <- 0; diag(m1)[is.na(m)] <- 1; LD_list_tmp[[l]] <- m1
    }
    indx1[[bl]] <- indx_tmp
    summ_list1[[bl]] <- summ_list_tmp
    snps_scale1[[bl]] <- snps_scale_tmp
    LD_list1[[bl]] <- LD_list_tmp

    rm(list=c("indx_tmp","summ_list_tmp","snps_scale_tmp","LD_list_tmp","snp_list_tmp","tmp","m1"))
  }

  summ_list <- summ_list1
  snps_scale <- snps_scale1
  LD_list <- LD_list1
  indx <- indx1
  indx_block <- indx_block1
  snp_list <- snp_list1
  Nsnps <- Nsnps1
  N <- N0
  #testing
  #print(N0)
  rm(list=c("summ_list1","snps_scale1","LD_list1","indx1","indx_block1","snp_list1","Nsnps1",
            "Nsnps0","snps_list0","snps_scale0","summ_list0","l","LD_list0","tmp"))
  if (opt$verbose == 2) cat("\n** Step 2.2 ended for chromosome ", chr, " **\n")
  
  ############
  ## Step 2.3. Run algorithm
  
  # Ensure verbose logging is enabled if verbose == 2
  if (opt$verbose == 2) cat("\n** Step 2.3 started for chromosome ", chr, " **\n")
        
  tryCatch({
    res <- gradient_descent_transfer_learning_all_blocks(
    summ_list,
    LD_list,
    M,
    indx,
    indx_block,
    n0 = N0[1],
    nk_list = N0[-1],
    alpha1 = alpha,
    alpha2 = alpha,
    alpha3 = alpha,
    alpha4 = alpha,
    eta_l = eta, #Use ADAM
    eta_m = eta,
    max_iter = iter,
    adaptive = TRUE,
    eta = 0.001,
    beta1 = 0.9,
    beta2 = 0.999,
    epsilon = 1e-8
  )  
    }, error = function(e) {
    cat("Error encountered during gradient descent for chromosome", chr, "\n")
    cat("Error message:", e$message, "\n")
  })
  
  if (opt$verbose == 2) cat("\n** Step 2.3 ended for chromosome ", chr, " **\n")
  
  ############
  ## Step 2.4. Clean PRSs into a vector
  if (opt$verbose == 2) cat("\n** Step 2.4 started for chromosome ", chr, " **\n")

  # Step 2.4. Clean PRSs into a matrix (#variant X #block)
  for (bl in 1:nblock) {
    tmp1 <- res[[bl]]$b  # Extract beta vector for block 'bl'
    tmp2 <- snps_scale[[bl]][,1]  # SNP scaling for block 'bl'
    dim(res[[bl]]$b);dim(snps_scale[[bl]])
    tmp2[is.na(tmp2)] <- 0    # Replace NA values in SNP scale with 0
    
    # Multiply beta vector (tmp1) by SNP scaling (tmp2)
    if (bl == 1) {
      b_tmp <- tmp1 * tmp2  # For the first block, initialize 'b_tmp'
    } else {
      b_tmp <- c(b_tmp, tmp1 * tmp2)  # Concatenate to build PRS vector
    }
  }
  
  # Final PRS matrix
  prs <- b_tmp
  #testing
  
  # Clean the PRS matrix
  prs[is.na(prs)] <- 0  # Set NA values to 0
  prs[prs > 10] <- 0    # Set values greater than 10 to 0
  prs[prs < -10] <- 0   # Set values less than -10 to 0
  #print(summary(prs))
  # Clean up temporary variables
  rm(list = c("tmp1", "tmp2", "b_tmp"))
  
  if (opt$verbose == 2) cat("\n** Step 2.4 ended for chromosome ", chr, " **\n")
  ############
  
  ############
  ## Step 2.6. Save files
  #
  # Note: 1. In the final prs file, the columns are: (rsid, a1: effect allele, a0: reference allele, PRSs...)
  #       2. For the param file, the order of its rows is same as the order of columns for PRSs. The param file indicate the tuning parameters and score source of the PRSs.
  
  snps <- unlist(snp_list)
  ref_tmp <- ref[match(snps, ref$V2),]
  df <- data.frame(rsid = snps, a1= ref_tmp$V5, a0= ref_tmp$V6, prs, stringsAsFactors=F)
  
  fwrite2(df, paste0(opt$PATH_out,"/tmp/PRS_in_all_settings_bychrom/prs_chr",chr,".txt"), col.names = F, sep="\t", nThread=1)
  
  rm(list=c("snps","snp_list","res", "df","prs"))
  
  progBar(ii, length(allchrom), per=5)
}

rm(list=c("df_beta_list"))
if (opt$verbose == 2) cat("\n** Step 2.6 ended **\n")
############
## Step 2.7. Combine all chromosomes
if (opt$verbose == 2) cat("\n** Step 2.7 started **\n")
score <- foreach(j = 1:length(allchrom), .combine='rbind') %dopar% {
  chr <- allchrom[j]
  prs <- fread2(paste0(opt$PATH_out,"/tmp/PRS_in_all_settings_bychrom/prs_chr",chr,".txt"))
  return(prs)
}
registerDoMC(1)
tmp <- score[,4] != 0; m <- !(tmp==0)
score <- score[m,,drop=F]
colnames(score) <- c("rsid","a1","a0",paste0("score",1:(ncol(score)-3)))
fwrite2(score, paste0(opt$PATH_out,"/before_ensemble/score_file.txt"), col.names = T, sep="\t", nThread=NCORES)

if ( opt$verbose >= 1 ) cat(paste0("PRS saved in ", opt$PATH_out,"/before_ensemble/score_file.txt \n"))
if (opt$verbose == 2) cat("\n** Step 2.7 ended **\n")
################
if(opt$testing){
  
  if ( opt$verbose >= 1 ) cat("\n** Step 3. Testing **\n")
  
  ############
  ## Step 3.1. Load data
  if (opt$verbose == 2) cat("\n** Step 3.1 started **\n")
  # Make/fetch the phenotype file
  fam <- read.table(paste(opt$bfile_testing,".fam",sep=''),as.is=T)
  if ( !is.na(opt$pheno_testing) ) {
    pheno <- read.table(opt$pheno_testing, as.is=T)
    # Match up data
    m <- match( paste(fam[,1],fam[,2]) , paste(pheno[,1],pheno[,2]) )
    m.keep <- !is.na(m)
    fam <- fam[m.keep,,drop=F]
    pheno <- pheno[m[m.keep],,drop=F]
  } else {
    pheno <- fam[,c(1,2,6)]
  }
  m <- is.na(pheno[,3]) # Remove samples with missing phenotype
  fam <- fam[!m,,drop=F]
  pheno <- pheno[!m,,drop=F]
  
  # Load in the covariates if needed
  if ( !is.na(opt$covar_testing) ) {
    covar <- ( read.table(opt$covar_testing,as.is=T,head=T) )
    if ( opt$verbose >= 1 ) cat( "Loaded",ncol(covar)-2,"covariates\n")
    # Match up data
    m <- match( paste(fam[,1],fam[,2]) , paste(covar[,1],covar[,2]) )
    m.keep <- !is.na(m)
    fam <- fam[m.keep,]
    pheno <- pheno[m.keep,]
    covar <- covar[m[m.keep],]
    reg <- summary(lm( pheno[,3] ~ as.matrix(covar[,3:ncol(covar)]) ))
    if ( opt$verbose >= 1 ) cat( reg$r.sq , "variance in phenotype explained by covariates in testing samples \n" )
    pheno[,3] <- scale(reg$resid)
  }
  if (opt$verbose == 2) cat("\n** Step 3.1 ended **\n")
  ############
  ## Step 3.2. Calculate scores for all tuning parameter settings on tuning samples
  if (opt$verbose == 2) cat("\n** Step 3.2 started  **\n")
  if (opt$verbose == 2) cat("\n** PLINK step started **\n")
  arg <- paste0(opt$PATH_plink ," --threads ",NCORES,
                " --bfile ",opt$bfile_testing,
                " --score ", opt$PATH_out,"/before_ensemble/score_file.txt header-read",
                " cols=+scoresums,-scoreavgs --score-col-nums 4",
                " --out ",opt$PATH_out,"/tmp/sample_scores_",ethnic[1],"/before_ensemble_testing")
  # system( arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
  system(arg, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (opt$verbose == 2) cat("\n** PLINK step ended **\n")
  if (!file.exists(paste0(opt$PATH_out,"/tmp/sample_scores_",ethnic[1],"/before_ensemble_testing.sscore"))) {
    stop("Error: PLINK output file does not exist. Check if PLINK ran successfully.")
  }
  SCORE <- fread2(paste0(opt$PATH_out,"/tmp/sample_scores_",ethnic[1],"/before_ensemble_testing.sscore"))
  
  m <- match( paste(fam[,1],fam[,2]) , paste(SCORE[,1],SCORE[,2]) )
  m.keep <- !is.na(m)
  fam <- fam[m.keep,]
  pheno <- pheno[m.keep,]
  SCORE <- SCORE[m[m.keep],]
  SCORE_id <- SCORE[,1:4]
  SCORE <- SCORE[,-1:-4]
  if (opt$verbose == 2) cat("\n** score vec derived **\n")
  
  # Get testing R2
  fit <- lm(pheno[,3]~SCORE)
  R2 <- summary(fit)$r.square
  # print(R2)
  
  # Store results
  results <- rbind(results, data.frame(iter = iter, eta = eta, alpha = alpha, R2 = R2))
  print(paste("Iteration:", iter, "| Eta:", eta, "| Alpha:", alpha, "| R²:", R2, "-> Appending results to the dataframe!"))
  # /dcs04/nilanjan/data/Anagh/PRS_proj/PROSPER/PROSPER_example_results/PROSPER/after_ensemble_AFR/R2.txt
}
    }
  }
}
write.csv(results, "R2_results.csv", row.names = FALSE)
if ( opt$verbose >= 1 ) cat(paste0("** !COMPLETED! \n"))

# if(opt$cleanup){
#   arg = paste0("rm -rf " , opt$PATH_out, "/tmp")
#   system(arg)
# }