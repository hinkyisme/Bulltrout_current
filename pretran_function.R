f <- function(x) {
  if (x[2] < x[1]) {
    Temp = x[1] - x[2]
  } else if (x[2] > x[1]) {
    Temp = (x[2] - x[1])*10
  } else {
    Temp = 1
  }
  (Temp)
}
# need a function that takes raster (temp)
# looks at first cell, compares to second cell
# if first cell is less than second cell, assign a value of 1
# if first cell is greater than second cell, assign a value of negative 1
# continue through whole raster until all values assigned, then run traditional transition 
## matrix stuff.

f <- function(x) {
  if (x[2] < x[1]) {
    R = x[1] - x[2]
  } else if (x[2] > x[1]) {
    R = (x[2] - x[1])*10
  } else {
    R = 1
  }
  (R)
}