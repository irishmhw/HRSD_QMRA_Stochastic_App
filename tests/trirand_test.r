require(testthat)
source('./TriRand.r')

test_that('creating TriRand floats in a loop and sapply are the same', {
    # this was a sanity check for converting everything to vector operations
    iter <- 100000
    a <- rtri(iter, 0.5, 1, 2.5)

    b <- matrix(ncol=1, nrow=iter)

    for (i in seq(1,iter)){
        b[i] <- TriRand(0.5, 1, 2.5)
    }
    
    c <- matrix(ncol=1, nrow=iter)
    
    for (i in seq(1,iter)){
        c[i] <- TriRand(0.5, 1, 2.5)
    }
    
    expect_equal(mean(b), mean(c), tolerance=1e-2)
    expect_equal(mean(a), mean(b), tolerance=1e-2)
})
