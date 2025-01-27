#Process alleles for the regression model
make_allele_data <- function(lin_alleles, allele_covar_mat, susc_types=NULL)
{
    n_alleles <- ncol(lin_alleles)
    n_covar <- nrow(allele_covar_mat)
    stopifnot(ncol(allele_covar_mat) == n_alleles)
    allele_mats <- apply(allele_covar_mat, 2, function(x) {
        n <- sum(x)
        m <- matrix(0L, nrow = n, ncol = n_covar)
        m[matrix(c(1:n, which(x > 0)), n , 2)] <- 1L       
        list(m) #necessary so that R doesn't attempt to reinterpret results as a vector if ncol=1
    })
    allele_mats <- unlist(allele_mats, recursive = FALSE) #undo additional list nesting
    names(allele_mats) <- paste0("allele_",1:n_alleles,"_interactions")
    n_interactions <- apply(allele_covar_mat, 2, sum)

    allele_levels <- NA
    allele_idx <- NA
    if (is.null(susc_types))
    {
        allele_idx <- apply(lin_alleles, 2, function(x) levels(factor(x)))
        lin_alleles <- as.matrix(apply(lin_alleles, 2, function(x) as.integer(factor(x))))
        allele_levels <- apply(lin_alleles, 2, max)
    } else 
    {
        stopifnot("Number of susceptible types must match number of alleles"=
            n_alleles==length(susc_types))
        allele_idx <- sapply(1:ncol(lin_alleles), function(x) levels(factor(lin_alleles[,x], exclude=susc_types[x])))
        lin_alleles <- as.matrix(sapply(1:ncol(lin_alleles), function(x) as.integer(factor(lin_alleles[,x], exclude=susc_types[x]))))
        lin_alleles <- ifelse(is.na(lin_alleles), 0L, lin_alleles)
        allele_levels <- apply(lin_alleles, 2, max)
    }

    allele_data <- data.frame(allele=1:n_alleles, n=n_interactions, m=allele_levels, k=n_covar)
    
    return(list(
        allele_mats=allele_mats,
        allele_data=allele_data,
        lin_alleles=lin_alleles,
        allele_idx=allele_idx
    )) 
}

make_grid <- function(t_max, covars, n_points)
{
    n_covar <- nrow(covars)-1
    grid_width <- t_max/n_points
    grid <- seq(from=grid_width, to=t_max, length.out=n_points)
    covars_interp <- matrix(data=0.0, nrow=n_covar, ncol=n_points)
    j <- 1
    for(i in seq_along(grid))
    {
        t_grid <- grid[i]
        while(covars[1, j] < (t_grid-1e-6)) 
        {
            j <- j+1        

        }
        covars_interp[, i] <- covars[-1, j]
    }
    return(list(grid_endpts=grid, 
        grid_width=grid_width,
        covars=covars_interp,
        K_covar=nrow(covars_interp)))
}

#Process trees into a poissonian representation
make_tree_data <- function(trees, mrsts, grid_endpts)
{
    N_tot_obs <- 0L
    N_tot_coal <- 0L
    lin_sizes <- c()
    lin_coal_counts <- c()

    wt <- c()
    aik <- c()
    ce <- c()
    grid_idx <- c()

    for (i in seq_along(trees))
    {
        evts <- extract_events(trees[[i]], mrsts[[i]], grid_endpts)
        wt <- c(wt, evts$wt)
        aik <- c(aik, evts$aik)
        ce <- c(ce, evts$ce)
        grid_idx <- c(grid_idx, evts$idx)
        
        lin_sizes <- c(lin_sizes, length(evts$wt))
        lin_coal_counts <- c(lin_coal_counts, length(evts$ce))
    }

    N_tot_obs <- sum(lin_sizes)
    N_tot_coal <- sum(lin_coal_counts)

    return(list(
        N_grid=length(grid_endpts),
        N_lin=length(trees),
        N_tot_obs=N_tot_obs,
        N_tot_coal=N_tot_coal,
        lin_sizes=lin_sizes,
        lin_coal_counts=lin_coal_counts,
        wt=wt, 
        aik=aik,
        ce=ce,
        grid_idx=grid_idx
    ))
}

#Project tree onto a grid
extract_events <- function(tree, mrst, grid_endpts) 
{
    #Extract tree events

    nsamp <- length(tree$tip.label)
    ncoal <- tree$Nnode
    nevt <- ncoal + nsamp
    ndepths <- node.depth.edgelength(tree)
    ndepths <- max(ndepths) - ndepths + mrst
    ntypes <- c(rep(1L, nsamp), rep(-1L, ncoal))

    node_ord <- order(ndepths)

    ndepths_inord <- ndepths[node_ord]
    ntypes_inord <- ntypes[node_ord]
    
    #Compute At
    At <- cumsum(ntypes_inord)
    At <- c(0, At[1:(nevt-1)])
    #Truncate tree to fit grid 
    trunc_pt <- max(grid_endpts)
    max_inbounds <- nevt
    
    exceeds_grid <- F
    At_oob <- NA
    if(ndepths_inord[nevt] >= trunc_pt)
    {
        #Include first event out of bounds so that grid cells get correctly interpolated
        max_inbounds <- max(which(ndepths_inord < trunc_pt))
        exceeds_grid <- T
        At_oob <- At[max_inbounds + 1L]
    } 

    ndepths_trunc <- ndepths_inord[1:max_inbounds]
    ntypes_trunc <- ntypes_inord[1:max_inbounds]
    At_trunc <- At[1:max_inbounds]

    #Combine grid events and tree events

    ngrid <- length(grid_endpts)

    grid_idx <- c(rep(NA, max_inbounds), 1:ngrid)
    evts_all <- c(ndepths_trunc, grid_endpts) 
    At_all <- c(At_trunc, rep(NA, ngrid))
    types_all <- c(ntypes_trunc, rep(1L, ngrid))
    
    all_ord <- order(evts_all)

    grid_idx <- grid_idx[all_ord]
    evts_all <- evts_all[all_ord]
    At_all <- At_all[all_ord]
    types_all <- types_all[all_ord]

    #Fill in grid cell index for tree events
    last_cell <- NA
    
    for(i in length(grid_idx):1) 
    {
        if(!is.na(grid_idx[i]))
        {
            last_cell <- grid_idx[i]
        } else 
        {
            stopifnot("Invalid grid."=!is.na(last_cell))
            grid_idx[i] <- last_cell
        }
    }

    #Fill in At for grid events
    last_At <- At_oob

    for (i in length(At_all):1)
    {
        if(!is.na(At_all[i]))
        {
            last_At <- At_all[i]
        } else 
        {
            At_all[i] <- last_At
        }
    }

    #Trim combined events to those that contribute to the likelihood
    #Last event is always a grid cell by construction
    last_evt <- max(which(ifelse(!is.na(At_all), At_all >= 2, FALSE)))

    grid_idx <- grid_idx[1:last_evt]
    evts_all <- evts_all[1:last_evt]
    At_all <- At_all[1:last_evt]
    types_all <- types_all[1:last_evt]

    #Trim first event
    first_evt <- min(which(ifelse(!is.na(At_all), At_all >= 2, FALSE))) - 1

    grid_idx <- grid_idx[first_evt:last_evt]
    evts_all <- evts_all[first_evt:last_evt]
    At_all <- At_all[first_evt:last_evt]
    types_all <- types_all[first_evt:last_evt]

    N <- length(grid_idx)
    
    wt <- diff(evts_all)
    aik <- choose(At_all[2:N], 2)
    idx <- grid_idx[2:N]
    ce <- which(types_all[2:N] == -1L)
    ts <- evts_all[2:N]
    
    return(list(
        ts=ts,
        wt=wt,
        aik=aik,
        idx=idx,
        ce=ce
    ))
}