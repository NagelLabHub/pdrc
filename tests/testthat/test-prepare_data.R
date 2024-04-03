
test_that("prepare_data function works correctly", {
  # Generate sample dataset
  dataset <- data.frame(
    Sample = c("001","002","003"),
    c_0 = c(40,50,60),
    c_15=c(25,27,30),
    c_30=c(20,21,22),
    c_60=c(15,17,18),
    c_120=c(12,11,10) )

  # Call the prepare_data function
  prepared_data <- prepare_data(dataset, "001")

  # Check if the output has the expected structure and values
  expect(nrow(prepared_data) >= 5, "Number of rows should be >= 5.")
  expect_equal(ncol(prepared_data), 3)  # Check the number of columns

  # Add more specific checks based on your function's expected output

  # Example checks for specific values
  expect_equal(prepared_data$time[1], 0)
})
