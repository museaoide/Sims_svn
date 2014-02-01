plotir <- function(resp, var, shock, type=c("multiple","single"),main="Impulse response plot", yax.flip=TRUE, vfirst=TRUE) {
    ## resp: impulse responses, as produced by impulsdtrf
    ## var: indexes of variables whose responses are to be plotted
    ## shock:  indexes of shocks whose effects are to be plotted
    ## yax.flip: Unless very few responses are plotted, vertical axis
    ##           labels start to overlap, unless yax.flip=TRUE
    ## vfirst: the variable index increases down each column.  Otherwise
    ##         (vfirst=FALSE) the shock index does so.
    ## As opposed to more direct plotting commands, this one has the
    ## advantage that (1) it scales all responses of the same variable
    ## to the same height on a given graph, (2) that it allows printing
    ## a matrix of graphs with rows corresponding to variables and columns
    ## to shocks, or vice versa, and (3) that it prompts for the legend
    ## (in the "single" case).
    ##-------------------------------
    ## Program seems to be unfinished or truncated somehow.  Incomplete repair 13.10.26.
    pr <- resp[var,shock, , drop=FALSE]
    ybound <- matrix(0,2,length(var))
    for (iv in 1:length(var)) ybound[ , iv] <- range(pr[iv, , ])
    if(vfirst) pr <- aperm(pr, c(2,1,3))
    nr <- dim(pr)[1]
    nc <- dim(pr)[2]
    np <- dim(pr)[3]
    layout(matrix(1:(nr * nc), nr, nc))
    for (ic in 1:nc) {
        if (nr > 1) {
            for (ir in 1:(nr-1))
                plot(0:(np-1), pr[ir, ic, ], type="l", ylim=ybound[ , if(vfirst) ir else ic], frame=TRUE)    
        }
    }
}
