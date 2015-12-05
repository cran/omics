mlm <- function(formula, data, vars, save.residuals=FALSE) {
    Y <- get(response.name(formula))

    formula <- update(formula, "NULL ~ .")
    mf <- model.frame(formula, data, na.action=na.omit, drop.unused.levels=TRUE)
    mm <- model.matrix(formula, mf)

    if (!is.null(attr(mf, "na.action"))) {
        Y <- Y[-attr(mf, "na.action"),]
    }

    labs <- labels(terms(formula))
    if (missing(vars)) {
        vars <- labs
    }

    idx <- which(attr(mm, "assign") %in% match(vars, labs))
    if (length(idx) == 0 & !save.residuals) {
        stop("No variables selected")
    }
    save.coefs <- !(length(idx) == 0 & save.residuals)
    vars <- colnames(mm)[idx]
    colnames(mm) <- sprintf("V%d", 1:ncol(mm))
    new.vars <- colnames(mm)[idx]
    mm <- as.data.frame(mm)
    formula <- as.formula(sprintf("y ~ %s - 1",
        paste0(colnames(mm), collapse=" + ")
    ))

    if (save.coefs) {
        coefs <- array(NA, c(ncol(Y), length(vars), 3),
                       dimnames=list(colnames(Y), vars,
                                     c("coef", "coef.se", "pval")))
    }

    if (save.residuals) {
        residuals <- matrix(NA, nrow(Y), ncol(Y), dimnames=dimnames(Y))
    }

    opts <- options(
        show.error.messages=FALSE,
        warn=2
    )
    on.exit(options(opts))

    for (i in 1:ncol(Y)) {
        mm$y <- Y[,i]

        model <- try(lm(formula, data=mm, na.action=na.exclude), silent=TRUE)
        if (inherits(model, "try-error")) {
            next
        }

        if (save.coefs) {
            tmp <- try(coef(summary(model)), silent=TRUE)
            if (inherits(tmp, "try-error")) {
                next
            }

            for (j in 1:length(vars)) if (new.vars[j] %in% rownames(tmp)) {
                coefs[i,vars[j],"coef"] <- tmp[new.vars[j],"Estimate"]
                coefs[i,vars[j],"coef.se"] <- tmp[new.vars[j],"Std. Error"]
                coefs[i,vars[j],"pval"] <- tmp[new.vars[j],"Pr(>|t|)"]
            }
        }

        if (save.residuals) {
            residuals[,i] <- resid(model, type="response")
        }
    }

    if (save.coefs & length(vars) == 1) {
        coefs <- as.data.frame(coefs[,1,])
    }

    tmp <- list()
    if (save.coefs) {
        tmp$coefficients <- coefs
    }
    if (save.residuals) {
        tmp$residuals <- residuals
    }
    tmp
}
