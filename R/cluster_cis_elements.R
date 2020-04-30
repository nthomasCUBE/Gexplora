function(){
	M=matrix(rep(0,length(ROWS)*length(COLS)),nrow=length(ROWS))
	rownames(M)=ROWS
	colnames(M)=COLS
	for(x in 1:dim(data)[1]){
		M[data[x,1],data[x,3]]=1
	}
	library(pheatmap)
	pheatmap(M)
}
