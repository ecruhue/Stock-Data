price <- read.csv("price_data.csv")

# calcuate the difference 
price.change <- apply(price, 2, diff.default)
price.change <- as.data.frame(price.change)

# get the 30% and 70% qunatile
fun <- function(x){
  quantile(x, c(0.3, 0.7))
}

quan.price <- apply(price.change, 2, fun)
quan.price <- t(quan.price)
quan.price

# transform the data 
# type get 30% -1, 40% 0 and 30% 1
trans.fun <- function(x){
  sort.price.change <- sort(x)
  low.bound <- quantile(sort.price.change, 0.3)
  up.bound <- quantile(sort.price.change, 0.7)
  ifelse(x > up.bound, 1, 
         ifelse(x > low.bound, 0 , -1))
}

trans.price <- apply(price.change, 2, trans.fun)
trans.price <- as.data.frame(trans.price)


