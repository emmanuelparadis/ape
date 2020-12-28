## by KS
replace.single.quotes <- function(x, start = 1L)
{
    z <- unlist(gregexpr("'", x))
    if (length(z) %% 2) {
        #warning("wrong number of single quotes around labels")
        warning(paste0("odd number of single quotes (", length(z), "): label(s) unchanged"))
        return(x)
    }
    i <- 1
    while (i < length(z)) {
        tmp <- substr(x, z[i], z[i + 1])
        substr(x, z[i], z[i + 1]) <- gsub("\\s+", "_", tmp)
        i <- i + 2
    }
    gsub("'", "", x)
}

"read.nexus.data" <- function (file)
{
    # Simplified NEXUS data parser.
    #
    # Version: 09/13/2006 01:01:59 PM CEST
    #          (modified by EP 2011-06-01 and by TG 2019-06-25)
    #
    # By:      Johan Nylander, nylander @ scs.fsu.edu
    #
    # WARNING: This is parser reads a restricted nexus format,
    #          see README for details.
    #
    # Argument (x) is a nexus formatted data file.
    #
    # Returns  (Value) a list of data sequences each made of a single
    #          vector of mode character where each element is a character.
    #
    # TODO:    Error checking, gap/missing, find.datatype, etc.
    #------------------------------------------------------------------

    "find.ntax" <- function (x)
    {
        for (i in 1:NROW(x)) {
            if(any(f <- grep("\\bntax", x[i], ignore.case = TRUE))) {
                ntax <- as.numeric(sub("(.+?)(ntax\\s*\\=\\s*)(\\d+)(.+)",
                                       "\\3", x[i], perl = TRUE, ignore.case = TRUE))
                break
            }
        }
        ntax
    }

    "find.nchar" <- function (x)
    {
        for (i in 1:NROW(x)) {
            if(any(f <- grep("\\bnchar", x[i], ignore.case = TRUE))) {
                nchar <- as.numeric(sub("(.+?)(nchar\\s*\\=\\s*)(\\d+)(.+)",
                                        "\\3", x[i], perl = TRUE, ignore.case = TRUE))
                break
            }
        }
        nchar
    }

    "find.matrix.line" <- function (x)
    {
        for (i in 1:NROW(x)) {
            if(any(f <- grep("\\bmatrix\\b", x[i], ignore.case = TRUE))) {
                matrix.line <- as.numeric(i)
                break
            }
        }
        matrix.line
    }

    "trim.whitespace" <- function (x)
    {
        gsub("\\s+", "", x)
    }

    "trim.semicolon" <- function (x)
    {
        gsub(";", "", x)
    }

    #TG: Added get polymorphism function
    "get.polymorphism" <- function (x)
    {
        ## Detect polymorphism function
        is.poly.start <- function(x) {return("(" == x || "{" == x)}
        is.poly.end   <- function(x) {return(")" == x || "}" == x)}
        ## Position increment
        position <- 1
        ## Check which position contains a polymorphism
        while(position <= length(x)) {
            ## Check whether the position is polymorphic
            if(is.poly.start(x[position])){
                ## Find the polymorphism end
                poly_end <- position + 1
                while(!is.poly.end(x[poly_end])) {
                    poly_end <- poly_end + 1
                    if(is.poly.start(x[poly_end]) || poly_end > length(x)) {
                        stop("missing closing bracket for a polymorphism at position ", position)
                    }
                }
                ## Replace the position by what's in the middle of the polymorphism
                x[position] <- paste0(x[(position+1):(poly_end-1)], collapse = "/")
                ## Remove the polymorphism
                x <- x[-c((position+1):poly_end)]
            }
            ## Increment the position
            position <- position + 1
        }
        return(x)
    }


    X <- scan(file = file, what = character(), sep = "\n",
              quiet = TRUE, comment.char = "[", strip.white = TRUE)
    ntax <- find.ntax(X)
    nchar <- find.nchar(X)
    matrix.line <- find.matrix.line(X)
    start.reading <- matrix.line + 1
    Obj <- list()
    length(Obj) <- ntax
    i <- 1
    pos <- 0
    tot.nchar <- 0
    tot.ntax <- 0

    ## by KS
    single.quotes <- grepl("'", X)
    if (any(single.quotes)) {
        to.replace <- which(single.quotes)
        for (j in to.replace) {
            X[[j]] <- replace.single.quotes(X[[j]])
        }
    }

    for (j in start.reading:NROW(X)) {
        Xj <- trim.semicolon(X[j])
        if(Xj == "") {
            break
        }
        if(any(jtmp <- grep("\\bend\\b", X[j], perl = TRUE, ignore.case = TRUE))) {
            break
        }
        ts <- unlist(strsplit(Xj, "(?<=\\S)(\\s+)(?=\\S)", perl = TRUE))
        if (length(ts) > 2) {
            stop("nexus parser does not handle spaces in sequences or taxon names (ts>2)")
        }
        if (length(ts) !=2) {
            stop("nexus parser failed to read the sequences (ts!=2)")
        }
        Seq <- trim.whitespace(ts[2])
        Name <- trim.whitespace(ts[1])
        nAME <- paste(c("\\b", Name, "\\b"), collapse = "")
        if (any(l <- grep(nAME, names(Obj)))) {
            tsp <- strsplit(Seq, NULL)[[1]]
            #TG: Convert polymorphisms
            if(any("(" %in% tsp || "{" %in% tsp)) {
                tsp <- get.polymorphism(tsp)
            }
            for (k in 1:length(tsp)) {
                p <- k + pos
                Obj[[l]][p] <- tsp[k]
                chars.done <- k
            }
        }
        else {
            names(Obj)[i] <- Name
            tsp <- strsplit(Seq, NULL)[[1]]
            #TG: Convert polymorphisms
            if(any("(" %in% tsp || "{" %in% tsp)) {
                tsp <- get.polymorphism(tsp)
            }
            for (k in 1:length(tsp)) {
                p <- k + pos
                Obj[[i]][p] <- tsp[k]
                chars.done <- k
            }
        }
        tot.ntax <- tot.ntax + 1
        if (tot.ntax == ntax) {
            i <- 1
            tot.ntax <- 0
            tot.nchar <- tot.nchar + chars.done
            if (tot.nchar == nchar*ntax) {
                print("ntot was more than nchar*ntax")
                break
            }
            pos <- tot.nchar
        }
        else {
            i <- i + 1
        }
    }
    if (tot.ntax != 0) {
        cat("ntax:",ntax,"differ from actual number of taxa in file?\n")
        stop("nexus parser did not read names correctly (tot.ntax!=0)")
    }
    for (i in 1:length(Obj)) {
        if (length(Obj[[i]]) != nchar) {
            cat(names(Obj[i]),"has",length(Obj[[i]]),"characters\n")
            stop("nchar differ from sequence length (length(Obj[[i]])!=nchar)")
        }
    }
    Obj <- lapply(Obj, tolower)
    Obj
}

## by KS:
nexus2DNAbin <- function(x)
{
    bs <- as.raw(._bs_)
    cs <- ._cs_
    fun <- function(x) {
        res <- as.raw(0)
        for (i in x) res <- res | bs[match(i, cs)]
        res <- res & bs[15] # bs[15] == n
        if (!(res %in% bs)) return(cs[17]) # return(?)
        cs[match(res, bs)]
    }
    y <- unique(unlist(x))
    y <- tolower(y)
    y <- y[is.na(match(y, cs))]
    if (length(y)) {
        tmp <- strsplit(y, "/")
        repl <- sapply(tmp, fun)
        for (i in seq_along(x)) {
            for (j in seq_along(y)) {
                x[[i]] <- gsub(y[j], repl[j], x[[i]])
            }
        }
    }
    as.DNAbin(x)
}
