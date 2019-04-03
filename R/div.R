div <- function
(file=NULL, sequence=NULL, sqs=NULL, method="LWL85", pairwise=TRUE, max.lim=3)
{
	stopifnot(
		!(is.null(file) && is.null(sequence)),
		method[1] %in% c("LWL85","JC69","K80","F81","K81","F84","BH87","T92","TN93","GG95"),
		requireNamespace("gtools", quietly=TRUE))

	if (!is.null(file)) {
		if (!requireNamespace("seqinr", quietly=TRUE)) {
			stop("Reading FASTA files require package seqinr")
		}		
	  if (!requireNamespace("ape", quietly=TRUE)) {
	    stop("Reading FASTA files require package ape")
	  }
	  if (method[1]=="LWL85") {
	    sequence <- seqinr::read.fasta(file)
	  }
	  else {
	  sequence=read.FASTA(file, type = "DNA")
		}
	}
  if (method[1]=="LWL85") {
    sequence <- .checkseq(sequence, gene.name=if (is.null(file)) "" else file)
  }
  else {
    sequence <- .checkseq2(sequence, gene.name=if (is.null(file)) "" else file)
  }
	if (is.null(sqs)) {
		sqs <- names(sequence)
	}
	combn <- gtools::combinations(n=length(sqs), r=2, v=sqs)
	if (method[1]=="LWL85") {
		return(data.frame(div=.LWL85(sequence, combn[,1], combn[,2], pairwise=pairwise, max.lim=max.lim), sp1=combn[,1], sp2=combn[,2]))
	} 
	else {
	  return(data.frame(div=.cal_div(sequence, combn[,1], combn[,2], pairwise=pairwise, method=method), sp1=combn[,1], sp2=combn[,2]))
	} 
	stop("Method ", method, " unknown.") # This should never happen
}
