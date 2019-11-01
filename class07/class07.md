class07
================
Jelani Lyda
10/23/2019

## Revisit our functions from last day

``` r
source("http://tinyurl.com/rescale-R")
```

Let’s try our rescale() function from last day

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
rescale( c(3, 10, 7))
```

    ## [1] 0.0000000 1.0000000 0.5714286

``` r
rescale2(c(3, 10, 4, 7))
```

    ## [1] 0.0000000 1.0000000 0.1428571 0.5714286

## Write a function both\_NA()

We want to write a function, called both\_na(), that counts how many
positions in two input vectors, x and y, both have a missing value

``` r
x <- c(1, 2, NA, 3, NA)
y <- c(NA, 2, NA, 3, 4)

is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
which( is.na(x))
```

    ## [1] 3 5

``` r
#Working snippit of code
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

``` r
both_na <- function(x,y){
  sum(is.na(x) & is.na(y))
}

both_na(x,y)
```

    ## [1] 1

``` r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA, NA)

both_na(x,y2)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 4

``` r
both_na2 <- function(x,y){
  if (length(x) != length(y)){
    stop("Vectors are different lenghts")
  }
  sum(is.na(x) & is.na(y))
}

both_na2(x, y1)
```

    ## [1] 2

\#\#Make a function to grade students

``` r
# student 1
a <- c(100, 100, 100, 100, 100, 100, 100, 90)
# student 2
b <- c(100, NA, 90, 90, 90, 90, 97, 80)

c <- b[-(which.min(b))]
mean(c, na.rm = TRUE)
```

    ## [1] 92.83333

``` r
grade <- function(data_set){
  data_set[is.na(data_set)] <- 0
  data_set2 <- data_set[-(which.min(data_set))]
  mean(data_set2, na.rm = TRUE)
}

grade(a)
```

    ## [1] 100

``` r
f <- c(1, 1, 1, 1)
f[is.na(f)] <- 0
f
```

    ## [1] 1 1 1 1
