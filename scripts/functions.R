# logging
log <- function(..., func = "cli_alert_info") {
  if (func == "cli_bullets") {
    msg2 <- list(...)
    msg1 <- msg2[[1]]
  } else {
    msg2 <- msg1 <- paste(...)
  }
  f <- get(func, envir = asNamespace("cli"))
  f(msg1)
  sink(log_file, append = TRUE)
  on.exit(sink(), add = TRUE)
  cat(cli::cli_fmt(do.call(func, as.list(msg2))), "\n")
}


# run a model
suppressPackageStartupMessages(library(lmerTest))
suppressPackageStartupMessages(library(lme4))
library(partR2)
library(rsq)
library(broom.mixed)
library(broom)
run_model <- function(dependent, fixed, model_dat, ns_func, main_study, data_type, platform, model_name, random=NULL, part_r2=NULL, extract=NULL, p=NULL) {

  # progress
  if(!is.null(p)) p()

  # checks
  stopifnot("`dependent` should be length 1" = length(dependent)==1)
  stopifnot("All of `dependent`, `fixed` and `random` must be columns in `model_dat`" = all(c(dependent, fixed, random) %in% names(model_dat)))
  ns_func <- match.arg(ns_func, choices = c("lmerTest::lmer", "stats::lm", "stats::glm", "lme4::glmer"))
  stopifnot("Mixed models must have a `random` effects variable" = (ns_func %in% c("lme4::glmer", "lmerTest::lmer") && !is.null(random)) || (ns_func %in% c("stats::lm", "stats::glm") && is.null(random)))


  # get the function
  nsfunc <- strsplit(ns_func, "::", fixed=TRUE)[[1]]
  ns     <- nsfunc[1]
  func   <- nsfunc[2]


  # extract data
  dat <- model_dat[, .SD, .SDcols = names(model_dat)[grepl(paste0("^", c(dependent, fixed, random), "$", collapse="|"), names(model_dat))]]
  dat <- dat[!is.na(get(dependent)), ]


  # create formula
  dep  <- paste0(dependent)
  fix  <- paste0(fixed, collapse=" + ")
  f    <- paste0(dep, " ~ ", fix)
  rand <- ifelse(!is.null(random), paste0("(1|", random, ")", collapse = " + "), NA)
  f    <- paste0(stats::na.omit(c(f, rand)), collapse = " + ")


  # get args
  args <- list(formula = f,
               REML    = ifelse(grepl("^lmer", func), TRUE, NA),
               family  = ifelse(grepl("^glm(er)?", func), "binomial", NA),
               data    = model_dat)
  args <- args[!is.na(args)]
  if(grepl("^glmer", func)) args$control <- glmerControl(optimizer="bobyqa")


  # log and try
  res <- tryCatch({

    fun     <- getFromNamespace(func, ns)
    fit     <- do.call(fun, args)
    rsq     <- tryCatch(rsq::rsq(fit, adj=FALSE, type="v"), error=function(e) { e; list(model=NA_real_, fixed=NA_real_, random=NA_real_) })
    adj_rsq <- tryCatch(rsq::rsq(fit, adj=TRUE, type="v"), error=function(e) { e; list(model=NA_real_, fixed=NA_real_, random=NA_real_) })
    r       <- if(func %in% c("glmer", "lmer")) broom.mixed::tidy(fit) else broom::tidy(fit)
    r       <- as.data.table(r)
    r[, `:=`(rsq_tot            = if(func %in% c("lm", "glm")) rsq else rsq$model,
             rsq_fixed          = if(func %in% c("lm", "glm")) rsq else rsq$fixed,
             rsq_random         = if(func %in% c("lm", "glm")) NA_real_ else rsq$random,
             adj_rsq_tot        = if(func %in% c("lm", "glm")) adj_rsq else adj_rsq$model,
             adj_rsq_fixed      = if(func %in% c("lm", "glm")) adj_rsq else adj_rsq$fixed,
             adj_rsq_random     = if(func %in% c("lm", "glm")) NA_real_ else adj_rsq$random,
             #rsq_timepoint      = if(is.null(part_r2)) NA_real_ else partR2(fit, data = dat, partvars = part_r2, R2_type = "marginal", nboot=10),
             #rsq_timepoint_cond = if(is.null(part_r2)) NA_real_ else partR2(fit, data = dat, partvars = part_r2, R2_type = "conditional", nboot=10),
             n                  = nrow(dat),
             main_study         = main_study,
             data_type          = data_type,
             platform           = platform,
             model_name         = model_name,
             formula            = f,
             error              = if(func %in% c("lmer", "glmer")) toString(fit@optinfo$conv$lme4$messages) else NA_character_)]

  }, error=function(e) {

    r <- data.table(term       = NA_character_,
                    estimate   = NA_real_,
                    p.value    = NA_real_,
                    n          = nrow(dat),
                    main_study = main_study,
                    data_type  = data_type,
                    platform   = platform,
                    model_name = model_name,
                    formula    = f,
                    error      = paste0(unlist(e), collapse = "; "))

  })

  # adjust p values
  res[p.value==0, p.value := .Machine$double.xmin]


  # pull out estimates
  if (!is.null(extract)) {
    regex <- paste0(extract, collapse="|")
    res <- res[grepl(regex, term), ]
  }

  return(res)

}
