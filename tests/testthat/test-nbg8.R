test_that("multiplication works", {
  m <- matrix(1:12, nrow = 3, ncol = 4, byrow = TRUE)

  expect_equal(ngb8(3,4)[[1]], c(m[1,2], m[2,1], m[2,2]))
  expect_equal(ngb8(3,4)[[8]], c(m[1,3], m[1,4], m[2,3], m[3,3], m[3,4]))
  expect_equal(ngb8(3,4)[[6]], c(m[1,1], m[1,2], m[1,3], m[2,1], m[2,3], m[3,1], m[3,2], m[3,3]))
})
