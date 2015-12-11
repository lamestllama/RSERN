generate = function(s, q, N, probFunc = 0, distFunc = 0, shape = 0)
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

    if(probFunc < 0 || probFunc > 8)
    {
        cat('probFunc is',probFunc,'\n')
        stop('Valid values are:
              waxman = 0,
              clipped_waxman = 1,
              waxman_transition_threshold =2,
              threshold = 3,
              constant = 4,
              powerlaw = 5,
              cauchy = 6,
              exponential = 7,
              maxentropy = 8')
    }

    if(distFunc < 0 || distFunc > 3)
    {
        cat('distFunc is',distFunc,'\n')
        stop('Valid values are:
              euclidean = 0,
              manhattan = 1,
              discrete = 2,
              maxdist = 3.')
    }

    if(shape < 0 || shape > 1)
    {
        cat('shape is',shape,'\n')
        stop('Valid values are:
              square = 0,
              circle = 1.')
    }



  tmp <- .Call("generate", c(s, q, N, probFunc, distFunc, shape))

	return(graph.data.frame(tmp[1], FALSE, tmp[2]))
}
