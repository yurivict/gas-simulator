buckets <- read.delim("~/gas-simulator/buckets.txt", sep=" ", header=F)

fn <- buckets

for(vv in buckets[1]) {
  r <- 1
  for (v in vv) {
    fn[r,2] <- fn[r,2]/(4*pi*fn[r,1]^2)*2900000 # 10 is for scale
    r <- r+1
  }
}

plot(buckets, type="o",
  xlab="Velocity",
  ylab="Number of particles"
)
lines(fn, type="o", col="red")
