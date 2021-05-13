library(blockTools)
library(tibble)

optgreed1 <- function(x, block.vars, id.vars, vcov, n.tr, l2, l1names, valid, validvar, validlb, validub, verbose, ismahal, dist){
  vcd <- as.matrix(x[, block.vars])
  if(nrow(vcd) == 1){
    result <- data.frame(matrix(c(x[, id.vars], rep(NA, n.tr)), ncol = n.tr + 1))
  }
  else{
    vcovi <- ginv(vcov)
    if(is.numeric(validvar) && length(validvar)==1){
      p = (length(unique(l1names)) %/% n.tr) + as.integer((length(unique(l1names)) %% n.tr) > 0)
    }
    else{
      p = nrow(x)
    }
  
    out <- .C("optgreed",
              data = as.double(vcd),
              vec = as.double(dist[row(dist) < col(dist)]),
              nrow = as.integer(nrow(x)),
              ncol = as.integer(ncol(vcd)),
              vcovi = as.double(vcovi),
              ntr = as.integer(n.tr),
              l2 = as.integer(l2),
              l1names = as.integer(as.factor(l1names)),
              valid = as.integer(valid),
              validvar = as.double(validvar),
              validlb = as.double(validlb),
              validub = as.double(validub),
              verbose = as.integer(verbose),
              pairdist = numeric(p),
              ismahal = as.integer(ismahal),
              result = integer(p * n.tr),
              p = as.integer(p))
    result <- data.frame(matrix(out$result, ncol=(n.tr), byrow=TRUE), out$pairdist)
  }
  return(result)
}

block1 <- function (data, vcov.data = NULL, groups = NULL, n.tr = 2, id.vars, 
    block.vars = NULL, algorithm = "optGreedy", distance = "mahalanobis", 
    weight = NULL, optfactor = 10^7, row.sort = NULL, level.two = FALSE, 
    valid.var = NULL, valid.range = NULL, seed.dist, namesCol = NULL, 
    verbose = FALSE, ...) 
{
    if (is.null(algorithm)) {
        stop("Blocking algorithm is unspecified.  See documentation for\noptions.")
    }
    if (!(algorithm %in% c("optGreedy", "optimal", "naiveGreedy", 
        "randGreedy", "sortGreedy"))) {
        stop("Blocking algorithm must be one of the following: optGreedy, optimal, naiveGreedy, randGreedy, sortGreedy")
    }
    if (algorithm == "randGreedy") {
        row.sort <- sample(seq(1:nrow(data)))
        data <- data[row.sort, ]
    }
    if (algorithm == "sortGreedy") {
        if (is.null(row.sort)) {
            stop("Blocking algorithm is 'sortGreedy', but vector of new row\npositions unspecified.  See documentation for details.")
        }
        if (length(row.sort) != nrow(data)) {
            stop("Length of vector 'row.sort' does not equal number of rows in\nthe data.  Respecify 'row.sort'.")
        }
        data <- data[row.sort, ]
    }
    if (is.matrix(data)) {
        data <- as.data.frame(data)
    }
    if (is.tibble(data)) {
        data <- as.data.frame(data)
    }
    if (is.null(block.vars)) {
        block.vars <- names(data)[!(names(data) %in% c(id.vars, 
            groups))]
    }
    if (!is.null(namesCol)) {
        if (level.two == FALSE && length(namesCol) != n.tr) {
            stop("The length of namesCol should equal n.tr")
        }
        else if (level.two == TRUE && length(namesCol) != 2 * 
            n.tr) {
            stop("When level.two == TRUE, the length of namesCol should equal 2*n.tr")
        }
    }
    if (is.null(vcov.data)) {
        vcov.data <- data[, block.vars]
    }
    else {
        vcov.data <- vcov.data[, block.vars]
    }
    vcov.data <- as.data.frame(vcov.data)
    if (ncol(vcov.data) == 1) {
        names(vcov.data) <- block.vars
    }
    if (is.character(distance)) {
        if (distance %in% c("mve", "mcd")) {
            iqr.idx <- 0
            while (iqr.idx < length(block.vars)) {
                iqr.idx <- iqr.idx + 1
                iqr.tmp <- unname(quantile(vcov.data[, iqr.idx], 
                  c(0.25, 0.75)))
                if (isTRUE(all.equal(iqr.tmp[1], iqr.tmp[2]))) {
                  warning(paste("Variable ", colnames(vcov.data)[iqr.idx], 
                    " has IQR 0; blocking will proceed using nonresistant Mahalanobis distance scaling matrix.", 
                    sep = ""))
                  distance <- "mahalanobis"
                  iqr.idx <- length(block.vars)
                }
            }
        }
        if (distance == "mahalanobis") {
            blvar.cut <- NULL
            for (blvar.idx in 1:length(block.vars)) {
                if (isTRUE(all.equal(var(vcov.data[, block.vars[blvar.idx]]), 
                  0))) {
                  blvar.cut <- append(blvar.cut, blvar.idx)
                }
            }
            if (length(blvar.cut) != 0) {
                warning(paste("The following blocking variables have zero variance and are dropped: ", 
                  paste(block.vars[blvar.cut], collapse = ", "), 
                  sep = ""))
                block.vars <- block.vars[-(blvar.cut)]
            }
            vcov.data <- vcov.data[, block.vars]
            vc.all <- var(vcov.data)
        }
        if (distance == "mcd") {
            vc.all <- cov.rob(vcov.data, method = "mcd", seed = seed.dist, 
                ...)$cov
        }
        if (distance == "mve") {
            vc.all <- cov.rob(vcov.data, method = "mve", seed = seed.dist, 
                ...)$cov
        }
        if (distance == "euclidean") {
            vc.all <- diag(ncol(vcov.data))
        }
    }
    if (!is.null(weight)) {
        if (is.vector(weight)) {
            if (length(weight) != ncol(vc.all)) {
                stop("Weight vector length must equal number of blocking variables.  Respecify 'weight'.")
            }
            weight <- diag(weight)
        }
        if (is.matrix(weight)) {
            if (sum(dim(weight) == dim(vc.all)) != 2) {
                stop("Weight matrix dimensions must equal number of blocking variables.  Respecify 'weight'.")
            }
        }
        vc.all <- solve(t(solve(t(chol(vc.all)))) %*% weight %*% 
            solve(t(chol(vc.all))))
    }
    gp <- 0
    out <- list()
    if (is.null(groups)) {
        data[, "groups"] <- 1
        groups <- "groups"
    }
    if (is.factor(data[, groups])) {
        data[, groups] <- as.character(data[, groups])
    }
    gp.names <- unique(data[, groups])
    for (i in unique(data[, groups])) {
        gp <- gp + 1
        if (verbose == TRUE) {
            cat("Blocking group ", i, "\n")
        }
        data.gp <- data[data[, groups] == i, c(id.vars, block.vars)]
        level.one.names <- data.gp[, id.vars[1]]
        if (is.factor(level.one.names)) {
            level.one.names <- as.character(level.one.names)
        }
        if (level.two == TRUE) {
            if (length(id.vars) < 2) {
                stop("Blocking requested at second level, but second level not\nidentified.  Specify a second ID variable and re-block.")
            }
            row.names(data.gp) <- data.gp[, id.vars[2]]
        }
        else {
            if (length(unique(data.gp[, id.vars[1]])) != length(data.gp[, 
                id.vars[1]])) {
                stop("Blocking requested at top level, but some units have\nidentical values of the identification variable.  Respecify first\nidentification variable and re-block.")
            }
            row.names(data.gp) <- data.gp[, id.vars[1]]
        }
        data.block <- data.frame(data.gp[, !(names(data.gp) %in% 
            id.vars)])
        if (is.null(valid.var)) {
            valid <- 0
            validvar <- numeric(1)
            validlb <- numeric(1)
            validub <- numeric(1)
        }
        else {
            if (is.null(valid.range)) {
                stop("A valid.var has been specified, but the valid.range has not.  Specify both or neither and re-block.")
            }
            valid <- 1
            validvar <- data.gp[, valid.var]
            validlb <- valid.range[1]
            validub <- valid.range[2]
        }
        if (!is.character(distance)) {
            dist.mat <- distance[data[, groups] == i, data[, 
                groups] == i]
        }
        else {
            nnn <- sum(as.integer(data[, groups] == i))
            dist.mat <- matrix(0, nrow = nnn, ncol = nnn)
        }
        if (algorithm != "optimal") {
            if (algorithm == "optGreedy") {
                out1 <- optgreed1(x = data.gp, block.vars = block.vars, 
                  id.vars = id.vars, dist = dist.mat, vcov = vc.all, 
                  n.tr = n.tr, l2 = level.two, l1names = level.one.names, 
                  valid = as.integer(valid), validvar = as.double(validvar), 
                  validlb = as.double(validlb), validub = as.double(validub), 
                  verbose = as.integer(verbose), ismahal = is.character(distance))
            }
            else if (algorithm %in% c("naiveGreedy", "randGreedy", 
                "sortGreedy")) {
                out1 <- naive(x = data.gp, block.vars = block.vars, 
                  id.vars = id.vars, vcov = vc.all, n.tr = n.tr, 
                  l2 = level.two, l1names = level.one.names, 
                  valid = as.integer(valid), validvar = as.double(validvar), 
                  validlb = as.double(validlb), validub = as.double(validub), 
                  verbose = as.integer(verbose), dist = dist.mat, 
                  ismahal = is.character(distance))
            }
        }
        if (algorithm == "optimal") {
            if (n.tr > 2) {
                warning("You specified algorithm = optimal and n.tr > 2.  However, optimal blocking only implemented for exactly two treatment conditions.  If no other error is encountered, optimal blocks for n.tr = 2 are returned here.")
            }
            if (is.character(distance)) {
                dist.mat <- mahal(data.block, vc.all)
            }
            else {
                dist.mat <- distance[data[, groups] == i, data[, 
                  groups] == i]
            }
            if (!is.null(valid.var)) {
                d.mat <- expand.grid(data.block[, valid.var], 
                  data.block[, valid.var])
                diffs <- abs(d.mat[, 1] - d.mat[, 2])
                valid.vec <- (valid.range[1] <= diffs) & (diffs <= 
                  valid.range[2])
            }
            dist.mat <- matrix(as.integer(optfactor * dist.mat), 
                nrow = nrow(dist.mat), ncol = ncol(dist.mat))
            if (!is.null(valid.var)) {
                warning("You specified algorithm = optimal and valid.var.  However, valid.var and valid.range are only implemented for other algorithms.  If no other error is encountered, optimal blocks ignoring the restriction are returned here.")
            }
            optimalDistanceMatrix <- nbpMatching::distancematrix(dist.mat)
            optimalOutput <- nbpMatching::nonbimatch(optimalDistanceMatrix, 
                precision = 9)
            optimalOutput$halves$Distance <- as.double(optimalOutput$halves$Distance)/optfactor
            out1 <- optimalOutput$halves
            out1 <- data.frame(`Unit 1` = out1$Group1.Row, `Unit 2` = out1$Group2.Row, 
                Distance = out1$Distance)
        }
        storage1 <- out1
        storage1[storage1 == 0 & col(storage1) < ncol(storage1)] <- NA
        count <- 1
        if (algorithm != "optimal") {
            for (col.idx in 1:(ncol(out1) - 1)) {
                storage1$temp <- as.character(data.gp[storage1[, 
                  col.idx], id.vars[1]])
                storage1$temp2 <- as.character(data.gp[storage1[, 
                  col.idx], id.vars[length(id.vars)]])
                names(storage1)[ncol(out1) + count] <- paste("Unit", 
                  col.idx)
                count <- count + 1
                names(storage1)[ncol(out1) + count] <- paste("Subunit", 
                  col.idx)
                count <- count + 1
            }
        }
        else if (algorithm == "optimal") {
            for (col.idx in 1:(ncol(out1) - 1)) {
                if (is.null(namesCol)) {
                  names(storage1)[col.idx] <- paste("Unit", col.idx)
                }
                else {
                  names(storage1)[col.idx] <- namesCol[col.idx]
                }
                storage1[, col.idx] <- as.character(data.gp[storage1[, 
                  col.idx], id.vars[1]])
            }
        }
        if (algorithm != "optimal") {
            storage1$Distance <- storage1[, ncol(out1)]
            storage <- storage1[, (ncol(out1) + 1):ncol(storage1)]
        }
        else if (algorithm == "optimal") {
            storage <- storage1
        }
        if (n.tr > 2) {
            names(storage)[ncol(storage)] <- "Max Distance"
        }
        odd.col <- seq(1:ncol(storage))[seq(1:ncol(storage))%%2 == 
            1]
        even.col <- (odd.col + 1)[1:(length(odd.col) - 1)]
        if (level.two == FALSE && algorithm != "optimal") {
            storage <- storage[, odd.col]
        }
        if (nrow(storage) > 1) {
            o <- order(as.numeric(as.character(storage[, ncol(storage)])), 
                storage[, 1])
            storage <- data.frame(storage[o, ], check.names = FALSE)
        }
        tmp.names <- unique(unlist(storage[, 1:(ncol(storage) - 
            1)]))
        left.out.names <- level.one.names[!(level.one.names %in% 
            tmp.names)]
        len.lon <- length(left.out.names)
        if (len.lon > 0) {
            left.out.df <- as.data.frame(matrix(c(left.out.names, 
                rep(NA, len.lon * (ncol(storage) - 1))), len.lon, 
                ncol(storage)))
            names(left.out.df) <- names(storage)
            storage.tmp <- data.frame(rbind(storage, left.out.df))
            names(storage.tmp) <- names(storage)
            storage <- storage.tmp
            rm(storage.tmp)
        }
        if (!is.null(namesCol)) {
            names(storage)[1:(ncol(storage) - 1)] <- namesCol
        }
        sum.na <- function(sum.na.vector) {
            return(sum(is.na(sum.na.vector)))
        }
        storage <- storage[apply(storage[, 1:(ncol(storage) - 
            1)], 1, sum.na) != (ncol(storage) - 1), ]
        rownames(storage) <- 1:(nrow(storage))
        out[[gp]] <- storage
    }
    names(out) <- gp.names
    o <- order(names(out))
    out <- out[o]
    output <- list(blocks = out, level.two = level.two)
    output$call <- match.call()
    class(output) <- "block"
    return(output)
}
