test_that("intraday cleaning and last-tick resampling follow documented rules", {
  data <- data.frame(
    date = c(rep("2026-01-05", 4), rep("2026-01-06", 4)),
    time = c("09:30:00", "09:30:02", "09:30:04", "09:30:06",
             "09:30:00", "09:30:02", "09:30:04", "09:30:06"),
    price = c(100, 101, 102, 103, 200, 200, 200, 200),
    stringsAsFactors = FALSE
  )
  cleaned <- clean_intraday_data(data, min_obs = 3L)
  expect_equal(unique(cleaned$data$date), "2026-01-05")
  expect_equal(cleaned$excluded$reason, "constant_price")
  coarse <- resample_last_tick(cleaned$data, interval_seconds = 5L)
  expect_equal(nrow(coarse), 2L)
  expect_equal(coarse$price, c(102, 103))
})
