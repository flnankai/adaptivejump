test_that("frictionless functions return finite p-values", {
  set.seed(42)
  sim <- simulate_jump_path(
    model = "dense", days = 1, seconds_per_day = 300,
    sampling_seconds = 1, theta = 1, mark_sd = 1
  )
  price <- exp(sim$latent_log_price)
  aj <- aj_test(price, k = 2L)
  lm <- lm_test(price, calibration = "bootstrap", bootstrap_rep = 20L, seed = 1L)
  combined <- adaptive_jump_test(
    price, k = 2L, lm_calibration = "bootstrap", bootstrap_rep = 20L, seed = 1L
  )
  expect_s3_class(aj, "jump_test")
  expect_s3_class(lm, "jump_test")
  expect_true(is.finite(aj$p_value))
  expect_true(is.finite(lm$p_value))
  expect_true(is.finite(combined$cc$p_value))
})

test_that("noisy functions estimate q and produce p-values", {
  set.seed(43)
  sim <- simulate_jump_path(
    model = "sparse", days = 1, seconds_per_day = 300,
    sampling_seconds = 1, lambda = 1, mark_sd = 0.2, noise_sd = 0.005
  )
  price <- exp(sim$observed_log_price)
  expect_true(is.finite(estimate_noise_sd(price)))
  ajj <- ajj_test(price, kn = 15L, rn = 75L)
  lm <- lm_noise_test(price, seconds_per_day = 300, bootstrap_rep = 20L, seed = 2L)
  expect_s3_class(ajj, "jump_test")
  expect_s3_class(lm, "jump_test")
  expect_true(is.finite(ajj$p_value))
  expect_true(is.finite(lm$p_value))
})

test_that("Cauchy combination and BH selection have valid outputs", {
  cc <- cauchy_combine(c(0.01, 0.20))
  selection <- bh_select(c(0.001, 0.03, 0.4), fdr = 0.2)
  expect_s3_class(cc, "jump_test")
  expect_true(cc$p_value > 0 && cc$p_value < 1)
  expect_length(selection$selected, 3L)
})
