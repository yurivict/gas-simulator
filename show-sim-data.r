buckets <- read.delim("~/gas-simulator/buckets.txt", sep=" ", header=F)

fn <- buckets

#for(i in buckets){for(j in i){print("--");print(j)}}
#for(i in buckets){print("--");print(i)}

for(vv in buckets[1]) {
  r <- 1
  for (v in vv) {
    fn[r,2] <- fn[r,2]/(4*3.141592*fn[r,1]*fn[r,1])*10 # 10 is for scale
    r <- r+1
  }
}

plot(buckets, type="o",
  xlab="Velocity",
  ylab="Number of particles"
)
lines(fn, type="o", col="red")
