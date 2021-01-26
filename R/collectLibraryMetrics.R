install.packages("BiocManager",repos = "http://cran.us.r-project.org")
library(BiocManager)
BiocManager::install("bamsignals")
library(IRanges)
library(bamsignals)
library(GenomeInfoDb)
library(GenomicRanges)
BiocManager::install("Rsamtools")
library(Rsamtools)
#collectLibraryStats("../../Downloads/aneufinder-master/backupBAM/")
args = commandArgs(trailingOnly=TRUE)


qc.spikiness <- function(counts) {
	if (is.null(counts)) {
		return(NA)
	}
	counts <- as.vector(counts)
	sum.counts <- sum(counts)
	spikiness <- sum(abs(diff(counts))) / sum.counts
	return(spikiness)
}

qc.entropy <- function(counts) {
	if (is.null(counts)) {
		return(NA)
	}
	counts <- as.vector(counts)
	total.counts <- sum(counts)
	n <- counts/total.counts
	entropy <- -sum( n * log(n) , na.rm=TRUE)
	return(entropy)
}

bamToGRanges <- function(bamfile, bamindex=bamfile,chromosomes=NULL,pairedEndReads=FALSE,remove.duplicate.reads=FALSE,min.mapq=10,max.fragment.width=1000,blacklist=NULL,what='mapq') {

	## Input checks
	if (!is.null(blacklist)) {
		if ( !(is.character(blacklist) | class(blacklist)=='GRanges') ) {
			stop("'blacklist' has to be either a bed(.gz) file or a GRanges object")
		}
	}

	## Check if bamindex exists
	bamindex.raw <- sub('\\.bai$', '', bamindex)
	bamindex <- paste0(bamindex.raw,'.bai')
	if (!file.exists(bamindex)) {
		ptm <- startTimedMessage("Making bam-index file ...")
		bamindex.own <- Rsamtools::indexBam(bamfile)
		warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
		bamindex <- bamindex.own
		stopTimedMessage(ptm)
	}
	chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(bamfile))
	chroms.in.data <- names(chrom.lengths)
	if (is.null(chromosomes)) {
		chromosomes <- chroms.in.data
	}
	chroms2use <- intersect(chromosomes, chroms.in.data)
	## Stop if non of the specified chromosomes exist
	if (length(chroms2use)==0) {
		chrstring <- paste0(chromosomes, collapse=', ')
		stop('The specified chromosomes ', chrstring, ' do not exist in the data. Pay attention to the naming convention in your data, e.g. "chr1" or "1".')
	}
	## Issue warning for non-existent chromosomes
	diff <- setdiff(chromosomes, chroms.in.data)
	if (length(diff)>0) {
		diffs <- paste0(diff, collapse=', ')
		warning(paste0('Not using chromosomes ', diffs, ' because they are not in the data.'))
	}

	## Import the file into GRanges
	ptm <- startTimedMessage("Reading file ",basename(bamfile)," ...")
	gr <- GenomicRanges::GRanges(seqnames=chroms2use, ranges=IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))
	if (!remove.duplicate.reads) {
		if (pairedEndReads) {
			data.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=what, mapqFilter = min.mapq))
		} else {
			data.raw <- GenomicAlignments::readGAlignments(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=what, mapqFilter = min.mapq))
		}
	} else {
		if (pairedEndReads) {
			data.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=what, mapqFilter = min.mapq, flag=scanBamFlag(isDuplicate=FALSE)))
		} else {
			data.raw <- GenomicAlignments::readGAlignments(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=what, mapqFilter = min.mapq, flag=scanBamFlag(isDuplicate=FALSE)))
		}
	}
	stopTimedMessage(ptm)

	if (length(data.raw) == 0) {
		if (pairedEndReads) {
			stop(paste0("No reads imported. Does your file really contain paired end reads? Try with 'pairedEndReads=FALSE'"))
		}
		stop(paste0('No reads imported! Check your BAM-file ', bamfile))
	}

	## Convert to GRanges and filter
	if (pairedEndReads) {
		ptm <- startTimedMessage("Converting to GRanges ...")
		data <- GenomicAlignments::granges(data.raw, use.mcols = TRUE, on.discordant.seqnames='drop') # treat as one fragment
		stopTimedMessage(ptm)

		ptm <- startTimedMessage("Filtering reads ...")
		# if (!is.na(min.mapq)) {
		# 	mapq.first <- mcols(GenomicAlignments::first(data.raw))$mapq
		# 	mapq.last <- mcols(GenomicAlignments::last(data.raw))$mapq
		# 	mapq.mask <- mapq.first >= min.mapq & mapq.last >= min.mapq
		# 	if (any(is.na(mapq.mask))) {
		# 		warning(paste0(bamfile,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NA' to keep all reads."))
		# 	}
		# 	data <- data[which(mapq.mask)]
		# }
		# Filter out too long fragments
		data <- data[width(data)<=max.fragment.width]
		stopTimedMessage(ptm)
	} else {
		ptm <- startTimedMessage("Converting to GRanges ...")
		data <- GenomicAlignments::granges(data.raw, use.mcols = TRUE) # treat as one fragment
		stopTimedMessage(ptm)

		ptm <- startTimedMessage("Filtering reads ...")
		# if (!is.na(min.mapq)) {
		# 	if (any(is.na(mcols(data)$mapq))) {
		# 		warning(paste0(bamfile,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NA' to keep all reads."))
		# 		mcols(data)$mapq[is.na(mcols(data)$mapq)] <- -1
		# 	}
		# 	data <- data[mcols(data)$mapq >= min.mapq]
		# }
		# Filter out too long fragments
		data <- data[width(data)<=max.fragment.width]
		stopTimedMessage(ptm)
	}

	## Exclude reads falling into blacklisted regions
	if (!is.null(blacklist)) {
		ptm <- startTimedMessage("Filtering blacklisted regions ...")
		if (is.character(blacklist)) {
			if (grepl('^chr', seqlevels(data)[1])) {
				chromosome.format <- 'UCSC'
			} else {
				chromosome.format <- 'NCBI'
			}
			black <- importBed(blacklist, skip=0, chromosome.format=chromosome.format)
		} else if (class(blacklist)=='GRanges') {
			black <- blacklist
		} else {
			stop("'blacklist' has to be either a bed(.gz) file or a GRanges object")
		}
		overlaps <- findOverlaps(data, black)
		idx <- setdiff(1:length(data), S4Vectors::queryHits(overlaps))
		data <- data[idx]
		stopTimedMessage(ptm)
	}

	return(data)

}

fixedWidthBins <- function(bamfile=NULL, assembly=NULL, chrom.lengths=NULL, chromosome.format, binsizes=1e6, stepsizes=NULL, chromosomes=NULL) {

	### Check user input ###
	if (length(binsizes) == 0) {
		return(list())
	}
	if (is.null(bamfile) & is.null(assembly) & is.null(chrom.lengths)) {
		stop("Please specify either a 'bamfile', 'assembly' or 'chrom.lengths'")
	}
	if (is.null(bamfile) & is.null(chrom.lengths)) {
		trigger.error <- chromosome.format
	}
	if (!is.null(stepsizes)) {
		if (length(stepsizes) != length(binsizes)) {
			stop("Need one element in 'stepsizes' for each element in 'binsizes'.")
		}
		if (any(binsizes < stepsizes)) {
			stop("'stepsizes' must be smaller/equal than 'binsizes'")
		}
	}

	### Get chromosome lengths ###
	if (!is.null(bamfile)) {
		ptm <- startTimedMessage(paste0("Reading header from ", bamfile, " ..."))
		chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(bamfile))
		stopTimedMessage(ptm)
	} else if (!is.null(assembly)) {
		if (is.character(assembly)) {
			ptm <- startTimedMessage("Fetching chromosome lengths from UCSC ...")
			df <- GenomeInfoDb::getChromInfoFromUCSC(assembly)
			stopTimedMessage(ptm)
		} else if (is.data.frame(assembly)) {
			df <- assembly
		} else {
			stop("Unknown assembly")
		}
		chrom.lengths <- df$size
		if (chromosome.format=='UCSC') {
		} else if (chromosome.format=='NCBI') {
			df$chrom = sub('^chr', '', df$chrom)
		}
		names(chrom.lengths) <- df$chrom
		chrom.lengths <- chrom.lengths[!is.na(names(chrom.lengths))]
		chrom.lengths <- chrom.lengths[!is.na(chrom.lengths)]
	} else if (!is.null(chrom.lengths)) {
		chrom.lengths <- chrom.lengths[!is.na(names(chrom.lengths))]
		chrom.lengths <- chrom.lengths[!is.na(chrom.lengths)]
	}
	chroms.in.data <- names(chrom.lengths)
	if (is.null(chromosomes)) {
		chromosomes <- chroms.in.data
	}
	chroms2use <- intersect(chromosomes, chroms.in.data)
	## Stop if none of the specified chromosomes exist
	if (length(chroms2use)==0) {
		chrstring <- paste0(chromosomes, collapse=', ')
		stop('Could not find length information for any of the specified chromosomes: ', chrstring, '. Pay attention to the naming convention in your data, e.g. "chr1" or "1".')
	}
	## Issue warning for non-existent chromosomes
	diff <- setdiff(chromosomes, chroms.in.data)
	if (length(diff)>0) {
		diffs <- paste0(diff, collapse=', ')
		warning('Could not find length information for the following chromosomes: ', diffs)
	}

	### Making fixed-width bins ###
	bins.list <- list()
	for (ibinsize in 1:length(binsizes)) {
		binsize <- binsizes[ibinsize]
		ptm <- startTimedMessage("Making fixed-width bins for bin size ", binsize, " ...")
		chrom.lengths.floor <- floor(chrom.lengths / binsize) * binsize
		bins <- unlist(GenomicRanges::tileGenome(chrom.lengths.floor[chroms2use], tilewidth=binsize), use.names=FALSE)
		bins <- bins[end(bins) > 0] # no chromosomes that are smaller than binsize
		if (any(width(bins)!=binsize)) {
			stop("tileGenome failed")
		}
		# seqlengths(bins) <- as.integer(chrom.lengths[names(seqlengths(bins))])
		seqlengths(bins) <- chrom.lengths[chroms2use]
		if (!is.null(stepsizes)) {
			shift.bp <- 0
			stepsize <- stepsizes[ibinsize]
			bins.list.step <- GRangesList()
			while (shift.bp < binsize) {
				bins.list.step[[as.character(shift.bp)]] <- suppressWarnings( trim(shift(bins, shift.bp)) )
				shift.bp <- stepsize + shift.bp
			}
			bins.list[[paste0('binsize_', format(binsize, scientific=TRUE, trim=TRUE), '_stepsize_', format(stepsize, scientific=TRUE, trim=TRUE))]] <- bins.list.step
		} else {
			bins.list[[paste0('binsize_', format(binsize, scientific=TRUE, trim=TRUE))]] <- bins
		}

		skipped.chroms <- setdiff(seqlevels(bins), as.character(unique(seqnames(bins))))
		if (length(skipped.chroms)>0) {
			warning("The following chromosomes were skipped because they are smaller than binsize ", binsize, ": ", paste0(skipped.chroms, collapse=', '))
		}
		stopTimedMessage(ptm)

	}

	return(bins.list)

}

startTimedMessage <- function(...) {

	x <- paste0(..., collapse='')
	message(x, appendLF=FALSE)
	ptm <- proc.time()
	return(ptm)

}

stopTimedMessage <- function(ptm) {

	time <- proc.time() - ptm
	message(" ", round(time[3],2), "s")

}

estimateComplexity <- function(reads) {
	message("Calculating complexity")
	downsample.sequence <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1) # downsampling sequence for MM approach
	vm <- vector()
	k <- vector()
	multiplicity <- list()
	total.reads.sans.dup <- vector()
	total.reads.unique <- vector()
	total.reads <- vector()
	for (p in downsample.sequence) {
		message("    p = ",p, appendLF=FALSE)
		if (p != 1) {
			down.reads <- reads[sort(sample(1:length(reads), size=p*length(reads), replace=FALSE))]
		} else {
			down.reads <- reads
		}
		sp <- start(down.reads)[as.logical(strand(down.reads)=='+')]
		sp1 <- c(sp[length(sp)], sp[-length(sp)])
		sm <- start(down.reads)[as.logical(strand(down.reads)=='-')]
		sm1 <- c(sm[length(sm)], sm[-length(sm)])

		if (length(sp)==1) {
			rlep <- rle(FALSE)
		} else {
			rlep <- rle(sp==sp1)
		}
		if (length(sm)==1) {
			rlem <- rle(FALSE)
		} else {
			rlem <- rle(sm==sm1)
		}
		tab.p <- table(rlep$lengths[rlep$values])	# table of number of duplicates
		names(tab.p) <- as.numeric(names(tab.p)) + 1
		tab.m <- table(rlem$lengths[rlem$values])
		names(tab.m) <- as.numeric(names(tab.m)) + 1
		multiplicities <- sort(as.numeric(c(1, union(names(tab.p), names(tab.m)))))
		m <- matrix(0, nrow=length(multiplicities), ncol=2)
		rownames(m) <- multiplicities
		m[names(tab.p),1] <- tab.p
		m[names(tab.m),2] <- tab.m
		dups <- apply(m, 1, sum)
		dups['1'] <- length(down.reads) - sum(dups*as.numeric(names(dups)))
		dups <- data.frame(multiplicity=as.numeric(names(dups)), frequency=dups)
		multiplicity[[as.character(p)]] <- dups

		total.reads.sans.dup[as.character(p)] <- dups[1,2]
		if (nrow(dups)>1) {
			total.reads.unique[as.character(p)] <- total.reads.sans.dup[as.character(p)] + sum(dups[2:nrow(dups),2])
		} else {
			total.reads.unique[as.character(p)] <- total.reads.sans.dup[as.character(p)]
		}
		total.reads[as.character(p)] <- length(down.reads)
	}
	message("")
	df <- data.frame(x=total.reads, y=total.reads.unique)

	## Complexity estimation with Michaelis-Menten
	complexity.MM <- NA
	complexity.MM.ggplt <- NA
	vm.init <- quantile(df$y, 1)
	k.init <- quantile(df$x, .25)
	tC <- tryCatch({
		complexity.MM.fit <- stats::nls(y ~ vm * x/(k+x), data=df, start=list(vm=vm.init, k=k.init))
		complexity.MM <- as.numeric(stats::coefficients(complexity.MM.fit)[1])
		max.x <- max(0.9 * stats::coefficients(complexity.MM.fit)['k'] / (1-0.9), max(total.reads))
		x <- seq(from=0, to=max.x, length.out=1000)
		df.fit <- data.frame(x=x, y=stats::predict(complexity.MM.fit, data.frame(x)))
		complexity.MM.ggplt <- ggplot(df) + geom_point(aes_string(x='x', y='y'), size=5) + geom_line(data=df.fit, mapping=aes_string(x='x', y='y')) + xlab('total number of reads') + ylab('unique reads') + theme_bw() + ggtitle('Michaelis-Menten complexity estimation')
	}, error = function(err) {
		# warning("Complexity estimation with Michaelis-Menten failed.")
	})

	rl <- list(complexity=c(MM=complexity.MM), MM.ggplt=complexity.MM.ggplt)
	return(rl)
}

collectLibraryStats <- function(folder){
	met=data.frame(file=c(),complexity=c(),coverage=c(),spikiness=c(),entropy=c())
	#file=list.files(folder,full.names = T,pattern="\\.bam$")[1]
	for (file in list.files(folder,full.names = T,pattern="\\.bam$")){
		bamindex=file
		bamindex.raw <- sub('\\.bai$', '', bamindex)
		bamindex <- paste0(bamindex.raw,'.bai')
		if (!file.exists(bamindex)) {
			ptm <- startTimedMessage("Making bam-index file ...")
			bamindex.own <- Rsamtools::indexBam(file)
			warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
			bamindex <- bamindex.own
			stopTimedMessage(ptm)
		}
		chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(file))
		chrom.lengths.df <- data.frame(chromosome=names(chrom.lengths), length=chrom.lengths)
		bins <- fixedWidthBins(chrom.lengths=chrom.lengths,  binsizes=1e+06, stepsizes=1e+06)
		counts.ss <- bamsignals::bamCount(file, bins[[1]][[1]], mapqual=10, paired.end="filter", tlenFilter=c(0,1000), verbose=FALSE, ss=TRUE, filteredFlag=1024)
		pcounts <- counts.ss[1,]
		mcounts <- counts.ss[2,]
		counts <- pcounts + mcounts
		readsperbin <- round(sum(as.numeric(counts)) / length(counts), 2)
		countmatrix <- matrix(c(counts,mcounts,pcounts), ncol=3)
		colnames(countmatrix) <- c('counts','mcounts','pcounts')
		mcols(bins[[1]][[1]]) <- as(countmatrix, 'DataFrame')

		data=bamToGRanges(file)
		complexity=estimateComplexity(data)[[1]]

		genome.length <- sum(as.numeric(seqlengths(data)))
		data.strand <- data
		strand(data.strand) <- '*'
		coverage <- sum(as.numeric(width(data.strand))) / genome.length
		#genome.covered <- sum(as.numeric(width(reduce(data.strand)))) / genome.length
		spikiness=qc.spikiness(bins[[1]][[1]]$counts)
		entropy=qc.entropy(bins[[1]][[1]]$counts)

		print(paste0("Complexity:",round(complexity,digits = -4)))
		print(paste0("Coverage:",coverage))
		#print(paste0("Genome.covered:",genome.covered))
		print(paste0("Spikiness",spikiness))
		print(paste0("Entropy:",entropy))

		row=data.frame(file=basename(file),complexity=complexity,coverage=coverage,spikiness=spikiness,entropy=entropy)
		met=rbind(met,row)
	}
	write.table(met,paste0("Input/NewMetrics/",basename(folder),".txt"),sep="\t",quote=F,row.names = F,col.names = T)
}

collectLibraryStats(args[1])



