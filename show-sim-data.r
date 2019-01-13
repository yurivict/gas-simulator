buckets <- read.delim("~/gas-simulator/buckets.txt", sep=" ", header=F)

fn <- buckets

for(vv in buckets[1]) {
  r <- 1
  for (v in vv) {
    fn[r,2] <- fn[r,2]/(4*pi*fn[r,1]^2)*2900000 # 10 is for scale
    r <- r+1
  }
}

#
# Maxwell distribution
#

# f(v)dv = √(m/2pikT)^3*4*pi*v^2*exp(-mv^2/(2kT))*dv

# assume homogeneous buckets distribution in velocity
dv <- buckets[3,1] - buckets[2,1]
maxwell <- buckets

# physical constants
Na = 6.022e+23
R = 8.3144598
k = R/Na

# params
N = 1000000
m = 6.646476410e-27
T = 293.321

# compute maxwell distribution
for(vv in buckets[1]) {
  r <- 1
  for (v in vv) {
    maxwell[r,2] <- N*((m/(2*pi*k*T))^(3/2))*(4*pi*v^2)*exp(-m*(v^2)/(2*k*T))*dv
    r <- r+1
  }
}

sum(buckets[2])
sum(maxwell[2])

plot(buckets, type="o",
  xlab="Velocity",
  ylab="Number of particles"
)
lines(fn, type="o", col="red")
lines(maxwell, type="o", col="green")
