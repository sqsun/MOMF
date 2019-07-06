#' Deconvoluation Analysis via Nonnegative Matrix Factorization
#'
#' The data set input \code{DataX} and \code{DataW} are necessary. \code{DataU} is needed if you do the multiple omics analysis
#' \code{DataH} is needed if you do the single omic analysis
#' Model fitting using the alternating direction method of multipliers (ADMM). 


#' Note: All data input should be nonnegative value, and the class should be list structure. All data sets have the same samples
#'
#' @param DataX nonnegative data matrix (list): bulk RNA-seq (#individuals x #genes) and scRNA-seq (#cells x #genes)
#' @param DataW initial nonnegative matrix (list)
#' @param DataH initial nonnegative matrix (list) when fitting single omic model
#' @param DataU initial nonnegative matrix (list) when fitting multiple omics model.
#' @param method two options: Kullback-Leibler (KL) divergence; Itakura-Saito (IS) divergence.
#' Default is "KL".
#' @param rho ADMM parameter, default is 2.
#'
#' @return a list.
#' @export
momf.fit <- function(DataX, DataW=NULL, DataH=NULL, DataU=NULL, DataPriorU=NULL, method="KL", rho=2, num_iter=5000){
	stopifnot(!missing(DataX))
	
	if(method=="KL"){beta_val<-1}
	if(method=="IS"){beta_val<-0}
	num_data <- length(DataX)
	# two divergences are provided currently
	if(method %in% c("KL", "IS") ){
		if(length(DataX)==1 & !is.null(DataH) & is.list(DataX) & is.null(DataPriorU)){# single data fitting
			stopifnot(!missing(DataW))
			cat(paste("## the details of inputs data (SOMF)\n") )
			cat(paste("## number of data sets: ", num_data,"\n") )
			cat(paste("## number of features: ", dim(DataX[[1]])[1],"\n") )
			cat(paste("## sample size: ", dim(DataX[[1]])[2],"\n") )
			cat(paste("## divergence: ", method,"\n") )
			
			cat(paste("## single omic data is fitting ...\n") )
			res_mf <- SOMF_cpp(DataX[[1]], DataW[[1]], DataH[[1]], beta_val, rho, num_iter)
			cat(paste("## finished! \n") )
			
			return( list("W"=res_mf$W, "H"=res_mf$H) )
		}else if(is.null(DataH) & !is.list(DataX) & is.null(DataPriorU)){
			stop("DataH should be provided/Data class should be list!")
		}# end fi
		
		# multiple omics data fitting
		if(length(DataX)>1 & !is.null(DataU) & is.list(DataX) & is.null(DataPriorU) ){
			stopifnot(!missing(DataW))
			cat(paste("## the details of inputs data (MOMF)\n") )
			cat(paste("## number of data sets: ", num_data,"\n") )
			cat(paste("## number of samples: ", ncol(DataX[[1]]),"\n") )
			for(idata in 1:num_data){
				stopifnot(ncol(DataX[[1]]) == ncol(DataX[[idata]]))
				cat(paste("## number of features for data set",idata,": ", nrow(DataX[[idata]]),"\n") )
			}
			cat(paste("## divergence: ", method,"\n") )
			cat(paste("## multiple omics data is fitting ...\n") )
			res_mf <- MOMF_cpp(DataX, DataW, DataU, beta_val, rho, num_data, num_iter)
			cat(paste("## finished! \n") )
	
			return(list("W"=res_mf$multiW, "U"=res_mf$U) )
		}else if(is.null(DataU) & !is.list(DataX) & is.null(DataPriorU)){
			stop("DataU should be provided/Data class should be 'list' structure!")
		}# end fi
		
		# with fix U
		if(length(DataX)>1 & is.list(DataX) & !is.null(DataPriorU) ){# deconvolution fitting
			if(is.null(DataW)){
				DataW <- list(W1 = DataX[[1]]%*%DataPriorU%*%inv(t(DataPriorU)%*%DataPriorU), 
					W2 = DataX[[2]]%*%DataPriorU%*%inv(t(DataPriorU)%*%DataPriorU))
			}
			if(is.null(DataU)){DataU <- t(DataPriorU)}
			num_data <- length(DataX)
			cat(paste("## the details of inputs data (MOMF)\n") )
			cat(paste("## number of data sets: ", num_data,"\n") )
			cat(paste("## number of genes: ", ncol(DataX[[1]]),"\n") )
			stopifnot(ncol(DataX[[1]]) == ncol(DataX[[2]]))
			cat(paste("## number of cells in scRNA-seq data: ", nrow(DataX[[1]]),"\n") )
			cat(paste("## number of individuals in bulk RNA-seq data: ", nrow(DataX[[2]]),"\n") )
			cat(paste("## divergence: ", method,"\n") )
			cat(paste("## deconvolution model is fitting ...\n") )
			res_mf <- MOMF_fixU_cpp(DataX, DataW, DataU, t(DataPriorU), beta_val, rho, num_data, num_iter)
			cat(paste("## finished! \n") )
			
			cell_prop <- t(apply(res_mf$multiW$W2, 1,  function(x) x/sum(x)))
			colnames(cell_prop) <- colnames(DataPriorU)
			rownames(cell_prop) <- rownames(DataX[[2]])
			cell_specific <- res_mf$U
			colnames(cell_specific) <- colnames(DataPriorU)
			rownames(cell_specific) <- rownames(DataPriorU)
			return(list("cell.prop"=cell_prop, "cell.specifc"=cell_specific) )
		}else if(is.null(DataU) & !is.list(DataX) & is.null(DataPriorU)){
			stop("DataU should be provided/Data class should be 'list' structure!")
		}# end fi
	}#end if
	
}# end funcs

#' Compute cell type specific mean in single cell RNA-seq data
#' @param sc_counts single-cell RNA-seq data 
#' @param cell_type the corresponding single cell labels (character)
#' @export
momf.computeRef <- function(sc_counts, cell_type){
	
	# cell type weight
    weight <- sapply(unique(cell_type), function(ct){
            y = sc_counts[ ,cell_type %in% ct, drop = FALSE]
            sum(y)/sum(sc_counts)
        })
    U <- sapply(unique(cell_type), function(ct){
            y = sc_counts[ ,cell_type %in% ct, drop = FALSE]
            rowSums(y)/sum(y)*weight[names(weight) == ct]
        })
	U <- 1000000 * U
	colnames(U) <- unique(cell_type)
	rownames(U) <- rownames(sc_counts)
	# return estimated mean for each cell type
	return(U)
}# end func



#' Compute cell type specific mean in single cell RNAseq data
#' @export
momf.ct_mean = function(DataX, cell_type, ct_select = NULL, filtering = TRUE, min_cell = 2, min_read = 3, verbose = TRUE){
	
	## filter out lowly expressed genes
	if(filtering){  
		selected_gene = rownames(DataX)[rowSums(DataX > min_read) > min_cell]
		DataX <- DataX[selected_gene, , drop = FALSE]
	}# end fi

	# relative abundance of gene g for cell type k
	ct_mean <- sapply(unique(cell_type), function(ct){
      y = DataX[ ,cell_type %in% ct, drop = FALSE]
      rowSums(y)/sum(y)
	})# end computing
	
	# average number of total mRNA molecules for cells of cell type k
	cell_size <- sapply(unique(cell_type), function(ct){
      y = DataX[ ,cell_type %in% ct, drop = FALSE]
      sum(y)/sum(cell_type %in% ct)
	})# end computing
	
	res <- sweep(ct_mean, 2, cell_size, "*")
	colnames(res) <- unique(cell_type)
	
	# select the cell type we are interested in
	if(!is.null(ct_select)){
		res <- res[, match(ct_select, colnames(res))]
	}
	# return estimated mean for each cell type
	return(res)
}# end func


