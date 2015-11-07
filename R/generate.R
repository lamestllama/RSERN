generate = function(s, q, N)
{
	if(!(is.numeric(s) && 
	     is.numeric(q) && 
	     is.numeric(N)))
	{
		cat('s is',class(s),', beta is',class(q),
		', N is',class(N),'\n')
		stop("all must be numeric")
	}
  
	if(s < 0)
	{
		cat('s is',s,'\n')
		stop('must be non-negative')
	}
  
	if(q < 0 || q > 1)
	{
		cat('q is',q,'\n')
		stop('must be in (0,1]')
	}
  
	if(N <=0)
	{
		cat('N is',N,'\n')
		stop('N must be greater than zero')
	}
	#Safety net
	if(s > 20)
	{
	 	cat("s too large\n")
	 	tmp <- list()
	 	return(tmp)
	}
  
  tmp <- .Call("generate", c(s, q, N))
  
	return(graph.data.frame(tmp[1], FALSE, tmp[2]))
}
