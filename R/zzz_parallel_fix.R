.prepare_adaptivejump_cluster <- function(workers) {
  cl <- parallel::makeCluster(workers)
  library_path <- dirname(find.package("adaptivejump"))
  parallel::clusterCall(cl, function(path) .libPaths(c(path, .libPaths())), library_path)
  parallel::clusterEvalQ(cl, library(adaptivejump))
  cl
}

# Replaces the initial definition so Windows socket workers explicitly inherit
# the library containing the installed package.
run_jump_simulation <- function(model = c("null", "dense", "sparse"),
                                test = c("frictionless", "noise"),
                                n_rep = 1000L,
                                simulation_args = list(),
                                test_args = list(),
                                alpha = c(0.10, 0.05),
                                workers = 1L,
                                seed = NULL) {
  model <- match.arg(model)
  test <- match.arg(test)
  n_rep <- as.integer(n_rep)
  workers <- as.integer(workers)
  if (n_rep < 1L || workers < 1L) {
    stop("`n_rep` and `workers` must be positive integers.", call. = FALSE)
  }
  if (workers == 1L) {
    if (!is.null(seed)) set.seed(as.integer(abs(seed) %% (.Machine$integer.max - 1L)))
    values <- lapply(seq_len(n_rep), function(i) {
      .simulation_one_rep(i, model, test, simulation_args, test_args)
    })
  } else {
    cl <- .prepare_adaptivejump_cluster(min(workers, n_rep))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if (!is.null(seed)) parallel::clusterSetRNGStream(cl, iseed = seed)
    values <- parallel::parLapplyLB(
      cl, seq_len(n_rep),
      function(i, model, test, simulation_args, test_args) {
        adaptivejump:::.simulation_one_rep(i, model, test, simulation_args, test_args)
      },
      model = model, test = test, simulation_args = simulation_args, test_args = test_args
    )
  }
  p_values <- as.data.frame(do.call(rbind, values))
  p_values$replication <- seq_len(n_rep)
  p_values <- p_values[c("replication", setdiff(names(p_values), "replication"))]
  methods <- setdiff(names(p_values), "replication")
  rates <- do.call(rbind, lapply(alpha, function(a) {
    data.frame(
      alpha = a,
      method = methods,
      rejection_rate = vapply(methods, function(m) mean(p_values[[m]] <= a), numeric(1L)),
      stringsAsFactors = FALSE
    )
  }))
  list(
    p_values = p_values,
    rejection_rates = rates,
    settings = list(model = model, test = test, n_rep = n_rep,
                    simulation_args = simulation_args, test_args = test_args,
                    alpha = alpha, workers = workers, seed = seed)
  )
}

daily_jump_tests <- function(data,
                             method = c("noise", "frictionless"),
                             test_args = list(),
                             alpha = c(0.01, 0.05, 0.10),
                             workers = 1L,
                             clean = TRUE) {
  method <- match.arg(method)
  workers <- as.integer(workers)
  if (workers < 1L) stop("`workers` must be positive.", call. = FALSE)
  cleaned <- if (clean) clean_intraday_data(data) else list(data = data, excluded = data.frame())
  daily <- split(cleaned$data, cleaned$data$date)
  daily <- daily[sort(names(daily))]
  if (!length(daily)) {
    return(list(p_values = data.frame(), excluded = cleaned$excluded, rejection_rates = data.frame()))
  }
  if (workers == 1L) {
    values <- lapply(daily, function(day) .daily_test_one(day, method, test_args))
  } else {
    cl <- .prepare_adaptivejump_cluster(min(workers, length(daily)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    values <- parallel::parLapplyLB(
      cl, daily,
      function(day, method, test_args) adaptivejump:::.daily_test_one(day, method, test_args),
      method = method, test_args = test_args
    )
  }
  p_values <- do.call(rbind, values)
  methods <- c("aj", "lm", "cc")
  rejection <- do.call(rbind, lapply(alpha, function(a) {
    data.frame(
      alpha = a,
      method = toupper(methods),
      rejection_rate = vapply(methods, function(m) mean(p_values[[paste0(m, "_p")]] <= a,
                                                          na.rm = TRUE), numeric(1L)),
      stringsAsFactors = FALSE
    )
  }))
  list(p_values = p_values, excluded = cleaned$excluded, rejection_rates = rejection)
}
