test_that("dchyper errors on n length", {
  expect_error(dchyper(3, 5, 6, 5, verbose = F))
})

test_that("dchyper errors on s length", {
  expect_error(dchyper(3, c(5,10), c(6, 10), c(4,8), verbose = F))
})

test_that("dchyper errors on m length not equal m length", {
  expect_error(dchyper(3, 5, c(6, 10, 8), c(4,8), verbose = F))
})

test_that("dchyper errors on negative values", {
  expect_error(dchyper(3, 5, c(6, -2), c(4,8), verbose = F))
})

test_that("dchyper errors on sample bigger than population", {
  expect_error(dchyper(3, 5, c(6, 2), c(4,8), verbose = F))
})

test_that("dchyper errors on non-integers", {
  expect_error(dchyper(3, 5, c(6, 4.5), c(4,8), verbose = F))
})

test_that("pchyper works", {
  expect_equal(pchyper(3, 5, c(6, 10), c(4,8), verbose = F), sum(dchyper(0, 5, c(6, 10), c(4,8), verbose = F),
                                                    dchyper(1, 5, c(6, 10), c(4,8), verbose = F),
                                                    dchyper(2, 5, c(6, 10), c(4,8), verbose = F),
                                                    dchyper(3, 5, c(6, 10), c(4,8), verbose = F)))
})

test_that("qchyper at 1 works", {
  expect_equal(qchyper(1, 5, c(6, 10), c(5,8), verbose = F), 5)
})

test_that("qchyper at 0 works", {
  expect_equal(qchyper(0, 5, c(6, 10), c(5,8), verbose = F), 0)
})

test_that("meanchyper works", {
  expect_equal(meanchyper(5, c(6, 10), c(5,8)), 120/99)
})

test_that("pvalchyper works", {
  expect_equal(1-pvalchyper(3, 5, c(6, 10), c(5,8), "upper", verbose = F),
               pchyper(2, 5, c(6, 10), c(5,8), verbose = F))
})

test_that("pvalchyper tail is upper or lower", {
  expect_error(pvalchyper(3, 5, c(6, 6), c(4,8), verbose = F, tail = "a"))
})

test_that("mleS works", {
  set.seed(1)
  expect_equal(mleS(rchyper(10^5, 5, c(6, 10), c(5,8), verbose = F),
                    c(6, 10), c(5,8), verbose = F), 5)
})

test_that("mleS input intersection not greater than sample size", {
  expect_error(mleS(6, c(6, 10), c(5,8), verbose = F))
})

test_that("mleN works", {
  set.seed(1)
  expect_equal(mleN(1, rchyper(10^5, 5, c(6, 10), c(5,8), verbose = F),
                    5, c(0, 10), c(5,8), verbose = F), 6)
})

test_that("mleM works", {
  set.seed(1)
  expect_equal(mleM(1, rchyper(10^5, 5, c(6, 10), c(5,8), verbose = F),
                    5, c(6, 10), c(0,8), verbose = F), 5)
})

test_that("momN works", {
  set.seed(1)
  expect_equal(momN(1, rchyper(10^5, 5, c(6, 10), c(5,8), verbose = F),
                    5, c(0, 10), c(5,8)), 6.01027533946055)
})

test_that("momM works", {
  set.seed(1)
  expect_equal(momM(1, rchyper(10^5, 5, c(6, 10), c(5,8), verbose = F),
                    5, c(6, 10), c(0,8)), 4.99533375)
})





