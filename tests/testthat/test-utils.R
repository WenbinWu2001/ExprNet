test_that("pirwin.hall works", {
  expect_equal(pirwin.hall(1, 2), 0.5)
  expect_equal(round(pirwin.hall(sqrt(2)/2, 2), 2), 0.25)
})

test_that("pirwin.hall invalid inputs", {
  expect_error(pirwin.hall("asd", 2))
  expect_error(pirwin.hall(1, c(1, 2)))
  expect_error(pirwin.hall(-0.1, 2))
  expect_error(pirwin.hall(2.5, 2))
  expect_error(pirwin.hall(2, 2.7))
})
