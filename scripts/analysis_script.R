library(tools)

args <- commandArgs(trailingOnly = TRUE)

# Genotype file (.tsv)
gen_file <- args[1]
gen_file_name <- file_path_sans_ext(gen_file)
gen_corr_file <- paste(gen_file_name, "corr.tsv", sep = "_")
gen_rqtl_file <- paste(gen_file_name, "rqtl.csv", sep = "_")
gen <- read.csv(gen_file,
                sep = "\t",
                header = TRUE,
                check.names = FALSE)
scn <- sapply(colnames(gen), function(x) { paste("snp", x, sep = "_") })
colnames(gen) <- scn
gen_t <- t(gen)
gen_t <- gen_t[rowSums(is.na(gen_t)) != ncol(gen_t), ]
gen <- gen[,colSums(is.na(gen)) < nrow(gen)]
gen <- apply(gen, c(1,2), function(x) {ifelse(x == 1, "B", "A")})
gen_col <- c()
for(i in 1:nrow(gen)) {
    nc <- paste("Sample", i, sep = "_")
    gen_col <- append(gen_col, nc)
}
colnames(gen_t) <- gen_col
write.table(gen_t,
            file = gen_corr_file,
            col.names = NA,
            sep = "\t"
            )
gen_rqtl <- cbind(data.frame(id = gen_col), gen)
gen_chr <- c(NA)
for(c in 1:ncol(gen)) {
    gen_chr <- append(gen_chr, 1)
}
gr <- rbind(gen_chr)
colnames(gr) <- colnames(gen_rqtl)
gg <- rbind(gr, gen_rqtl)
write.table(gg,
            file = gen_rqtl_file,
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            na = "",
            sep = ","
            )

# Expression file (.tsv)
expr_file <- args[2]
expr_file_name <- file_path_sans_ext(expr_file)
expr_corr_file <- paste(expr_file_name, "corr.tsv", sep = "_")
expr_rqtl_file <- paste(expr_file_name, "rqtl.csv", sep = "_")
expr <- read.csv(expr_file,
                 sep = "\t",
                 header = TRUE,
                 check.names = TRUE)
ecn <- sapply(colnames(expr), function(x) { paste("gene", x, sep = "_") })
colnames(expr) <- ecn
expr_t <- t(expr)
expr_t <- expr_t[rowSums(is.na(expr_t)) != ncol(expr_t), ]
expr_col = c()
for(i in 1:nrow(expr)) {
    nc <- paste("Sample", i, sep = "_")
    expr_col <- append(expr_col, nc)
}
colnames(expr_t) <- expr_col
write.table(expr_t,
            file = expr_corr_file,
            col.names = NA,
            sep = "\t"
            )

expr_rqtl <- cbind(data.frame(id = expr_col), expr)
expr_rqtl <- expr_rqtl[,colSums(is.na(expr_rqtl)) < nrow(expr_rqtl)]
write.table(expr_rqtl,
            file = expr_rqtl_file,
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = ","
            )

# Output prefix
eqtl_out_prefix <- args[3]

# tmp folder
tmp_dir <- args[4]

###############################################################
# Matrix eQTL
library("MatrixEQTL")
###############################################################
## Settings

# Covariates file name
# Set to character() for no covariates
covariates_file_name <- character()

# Only associations significant at this level will be saved
pvOutputThreshold <- 1e-2

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance <- numeric()

## Load genotype data

snps <- SlicedData$new()
snps$fileDelimiter <- "\t"      # the TAB character
snps$fileOmitCharacters <- "NA" # denote missing values;
snps$fileSkipRows <- 1          # one row of column labels
snps$fileSkipColumns <- 1       # one column of row labels
snps$fileSliceSize <- 2000      # read file in slices of 2,000 rows
snps$LoadFile(gen_corr_file)

## Load gene expression data

gene <- SlicedData$new()
gene$fileDelimiter <- "\t"      # the TAB character
gene$fileOmitCharacters <- "NA" # denote missing values;
gene$fileSkipRows <- 1          # one row of column labels
gene$fileSkipColumns <- 1       # one column of row labels
gene$fileSliceSize <- 2000      # read file in slices of 2,000 rows
gene$LoadFile(expr_corr_file)

## Load covariates

cvrt <- SlicedData$new()
cvrt$fileDelimiter <- "\t"      # the TAB character
cvrt$fileOmitCharacters <- "NA" # denote missing values;
cvrt$fileSkipRows <- 1          # one row of column labels
cvrt$fileSkipColumns <- 1       # one column of row labels
if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name)
}

## Run the analysis LINEAR

# Output file name
output_file_name <- paste(eqtl_out_prefix,
                          "MatrixEQTL_LINEAR.tsv",
                          sep = "_")

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel <- modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

me <- Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
)


## Run the analysis ANOVA

# Output file name
output_file_name <- paste(eqtl_out_prefix,
                          "MatrixEQTL_ANOVA.tsv",
                          sep = "_")

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel <- modelANOVA; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

me <- Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
)

###############################################################
# mRMR program (windows version, mrmr_win32.exe; Linux version, mrmr) is 
# downloaded from http://penglab.janelia.org/proj/mRMR/
###############################################################

# The genotype data
snp.dat <- gen_t
# The gene expression data
gene.dat <- expr_t

colnames(expr_t)

# The output eQTL file name
output_file_name <- paste(eqtl_out_prefix,
                          "mRMR.tsv",
                          sep = "_")

List <- 'MaxRel'
#List<- 'mRMR'

# The maximal number of features to be analyzed
feature_count <- 10

# The cutoff of accuracy
Accuracy.cutoff <- 0.6 

table(colnames(gene.dat) == colnames(snp.dat))

dis <- function(x,y){
    1 - crossprod(x, y)/sqrt(crossprod(x) * crossprod(y))
}

loocv.NNA <- function(dat,cls) {
    predict.cls <- rep(NA, length(cls))
    for(j in 1:nrow(dat)) {
        dat.test <- dat[j,]
        dat.train <- dat[-j,]
        cls.tarin <- cls[-j]
        dis.train <- apply(dat.train,1,function(x){dis(x,dat.test)})
        cls.test <- table(cls.tarin[dis.train == min(dis.train)])
        predict.cls[j] <- names(cls.test)[cls.test == max(cls.test)][1]
    }
    predict.cls <- as.numeric(predict.cls)
    acc <- sum(predict.cls == cls) / length(cls)
    acc.list <- acc
    for(x in names(table(cls))) {
        acc.x <- sum((predict.cls == x) & (cls == x)) / sum(cls == x)
        acc.list <- c(acc.list, acc.x)
    }
    acc.list
}

mrmr_prog <- paste(args[5], "mrmr", sep="/")
dir.create(tmp_dir, showWarnings= FALSE)

for(i in 1:nrow(snp.dat)) {
    snp.id <- rownames(snp.dat)[i]
    snp.cls <- as.numeric(snp.dat[snp.id,])
    cat(paste("Computing", "-", snp.id, sep = " "))
    cat("\n")
    in_file <- paste(eqtl_out_prefix,#args[3],
		     snp.id,
                     'input.csv',
                     sep='.')
    write.table(cbind(snp.cls,t(gene.dat)),
                file = in_file,
                row.names = FALSE,
                quote = FALSE,
                sep = ',')
    cmd <- paste(mrmr_prog,
                 ' -i ',
                 in_file,
                 ' -t 1 -n ',
                 feature_count,
                 ' -s ',
                 ncol(gene.dat),
                 ' -v ',
                 nrow(gene.dat),
                 sep = '')
    tmp.output <- system(cmd,
                         intern = TRUE,
                         ignore.stderr = TRUE)
    unlink(in_file)
    if(List == 'mRMR') {
        strat.line <- grep("*** mRMR features *** ",
                           tmp.output,
                           fixed = TRUE)+1
        end.line <- grep(" *** This program and the respective minimum Redundancy Maximum Relevance (mRMR) ",
                         tmp.output,
                         fixed = TRUE)-1
    }
    if(List == 'MaxRel') {
        strat.line <- grep("*** MaxRel features ***",
                           tmp.output,
                           fixed = TRUE)+1
        end.line <- grep("*** mRMR features *** ",
                         tmp.output,
                         fixed = TRUE)-1
    }
    txt_file <- paste(tmp_dir,
                      snp.id,
                      '.mRMR.txt',
                      sep = '')
    csv_file <- paste(tmp_dir,
                      snp.id,
                      '.mRMR.IFS.csv',
                      sep = '')
    sink(txt_file)
    cat(tmp.output[strat.line:end.line],
        file = txt_file,
        sep = '\n')
    sink()
    fea.list <- as.character(read.table(txt_file,
                                        header = TRUE)$Name)
    sink(csv_file)
    cat('Number.of.feature',
        'Accuracy',
        paste('Accuracy', names(table(snp.cls)), sep = '.'),
        sep = ',')
    cat('\n', sep = '')
    for(j in 2:length(fea.list)) {
        fea <- fea.list[1:j]
        performance <- loocv.NNA(t(gene.dat[fea,]), snp.cls)
        cat(j, performance, sep = ',')
        cat('\n', sep = '')
    }
    sink()
}

snp.list <- rownames(snp.dat)
print(length(snp.list))

eQTL.table <- NULL
for(i in 1:length(snp.list)) {
    snp <- snp.list[i]
    txt_file <- paste(tmp_dir,
                      snp,
                      '.mRMR.txt',
                      sep = '')
    csv_file <- paste(tmp_dir,
                      snp,
                      '.mRMR.IFS.csv',
                      sep = '')
    mRMR <- read.delim(txt_file,
                       header = TRUE)
    Accuracy <- read.csv(csv_file,
                         header = TRUE)
    Accuracy <- rbind(c(1,
                        rep(NA, ncol(Accuracy)-1)),
                      Accuracy)
    tab <- cbind(SNP = snp,
                 Order = mRMR[,'Order'],
                 Gene = gsub(' ',
                             '',
                             mRMR[,'Name']),
                 Score = mRMR[,'Score'],
                 Accuracy = Accuracy[,'Accuracy'])
    index.i <- which(as.numeric(tab[,'Accuracy']) >= Accuracy.cutoff)
    if(length(index.i) != 0) {
        eQTL.table <- rbind(eQTL.table,tab[1:min(index.i),])
    }
}

write.table(eQTL.table,
            file = output_file_name,
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

unlink(tmp_dir, recursive = TRUE)

###############################################################
# r/QTL analysis
library(qtl)
###############################################################

rqtl_data <- read.cross("csvs",
                        ".",
                        gen_rqtl_file,
                        expr_rqtl_file,
                        estimate.map = FALSE)

## rqtl_data <- est.rf(rqtl_data)
## pdf("RecombFactorPair.pdf")
## plotRF(rqtl_data)
## dev.off()

rqtl_data <- calc.genoprob(rqtl_data,
                           step = 0)


# EM algorithm
output_file_name <- paste(eqtl_out_prefix,
                          "rqtl_em.csv",
                          sep = "_"
                          )
out.em <- scanone(rqtl_data,
                  pheno.col = 2:nphe(rqtl_data)
                  )
write.table(out.em,
            file = output_file_name,
            col.names = NA,
            row.names = TRUE,
            quote = FALSE,
            sep = ","
            )

output_file_name <- paste(eqtl_out_prefix,
                          "rqtl_hk.csv",
                          sep = "_"
                          )
out.hk <- scanone(rqtl_data,
                  pheno.col = 2:nphe(rqtl_data),
                  method = "hk")
write.table(out.hk,
            file = output_file_name,
            col.names = NA,
            row.names = TRUE,
            quote = FALSE,
            sep = ","
            )


## pdf("Test3.pdf")
## plot(out.em, out.hk, out.imp, col=c("blue", "red", "green"))
## plot(out.imp - out.em, out.hk - out.em, col=c("green", "red"), ylim=c(-1,1))
## dev.off()
