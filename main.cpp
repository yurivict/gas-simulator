//
// Copyright (C) 2019 by Yuri Victorovich. All rights reserved.
//

#include <stdlib.h>

#include <string>
#include <array>
#include <list>
#include <vector>
#include <set>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include <ctime>
#include <limits>

#include <cxxopts.hpp>
#include <nlohmann/json.hpp>
#include <boost/format.hpp>

#if USE_PARALLELISM
#include <thread>
#include "ThreadPool.h"
#endif

#include "Vec3.h"

//
// types
//

typedef double Float;

//
// physics constants and laws
//

static constexpr Float Na = 6.022e+23;  // Avogadro constant
static constexpr Float R = 8.3144598;   // 8.3144598(48) in kg m² s⁻² K−1 mol−1 , see https://en.wikipedia.org/wiki/Gas_constant
static constexpr Float k = R/Na;        // Boltzmann constant, see https://en.wikipedia.org/wiki/Boltzmann_constant
static constexpr Float Troom = 293.15;  // room temperature in Kelvins: 20°C=293.15K
static constexpr Float Patm = 101325;   // 101.325 kPa

// quick gas computations:
// Vmol(0°C)  = 1mol*R*T=273.15K/101325 ≈ 0.0224m³ = 22.4dm³
// Vmol(20°C) = 1mol*R*T=293.15K/101325 ≈ 0.0241m³ = 24.1dm³
// V(1mil, 20°C) = (n=1mil/Na)RT/Patm ≈ 0.342e-6 (~0.3μm)

#if USE_PARALLELISM
static ThreadPool *threadPool = nullptr; // static thread pool user by all code here
#endif

#if DBG_TRACK_PARTICLES
#define DBG_TRACK_PARTICLE_MSG(msg) std::cout << "...[track]: " << msg << std::endl;
#define DBG_TRACK_PARTICLE_CODE(cmds ...) cmds
#else
#define DBG_TRACK_PARTICLE_MSG(msg) {/*do nothing*/}
#define DBG_TRACK_PARTICLE_CODE(cmds ...) {/*do nothing*/}
#endif

#if USE_PARALLELISM
#define TID "[" << ThreadPool::tid() << "]"
#else
#define TID ""
#endif

//
// helper classes
//

class AOut : public std::ostringstream { // atomic message printer
public:
  ~AOut() {
    std::cout << str();
  }
};

//
// helper functions
//

namespace constexpr_funcs {

// borrowed from the answer in https://stackoverflow.com/questions/23781506/compile-time-computing-of-number-of-bits-needed-to-encode-n-different-states/23795363
constexpr unsigned floorlog2(unsigned x) {
  return x == 1 ? 0 : 1+floorlog2(x >> 1);
}

constexpr unsigned ceillog2(unsigned x) {
  return x == 1 ? 0 : floorlog2(x - 1) + 1;
}

namespace Detail {
  static double constexpr sqrtNewtonRaphson(double x, double curr, double prev) {
    return curr == prev
      ? curr
      : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
  }

  static double constexpr cbrtNewtonRaphson(double x, double curr, double prev) {
    return curr == prev
      ? curr
      : cbrtNewtonRaphson(x, curr - (curr*curr*curr - x)/(3*curr*curr), curr);
  }
}

static double constexpr sqrt(double x) {
  return x >= 0 && x < std::numeric_limits<double>::infinity()
    ? Detail::sqrtNewtonRaphson(x, x, 0)
    : std::numeric_limits<double>::quiet_NaN();
}

static double constexpr cbrt(double x) {
  return x >= 0 && x < std::numeric_limits<double>::infinity()
    ? Detail::cbrtNewtonRaphson(x, x, 0)
    : std::numeric_limits<double>::quiet_NaN();
}

unsigned constexpr lim1(unsigned i) {
  return i >= 1 ? i : 1;
}

static Float constexpr temperatureToEnergy(Float T) {return k*T;}
static Float constexpr energyToTemperature(Float E) {return E/k;}

}; // constexpr_funcs

namespace xasm {

  static inline uint64_t getCpuCycles(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
  }

}; // xasm

static void formatUInt64(uint64_t x, std::ostream &os) {
  if (x >= 1000) {
    formatUInt64(x/uint64_t(1000), os);
    os << ',';
    os << std::setw(3) << std::setfill('0') << x%1000;
  } else {
    os << x;
  }
}

static std::string formatUInt64(uint64_t x) {
  std::ostringstream ss;
  formatUInt64(x, ss);
  return ss.str();
}

//
// relevant physics laws
//

// Ideal gas law: https://en.wikipedia.org/wiki/Ideal_gas_law
// PV = nRT    // n is number of moles
// E = kT      // relation between the average kinetic energy and temperature

//
// params
//
static constexpr unsigned N = 1000000;
static constexpr unsigned Ncreate = N; // (DEBUG) should be =N, but allows to override the number of particles for debugging
static constexpr Float particlesPerBucket = 1.; // average number of particles per bucket. Performance drops when the number is >1 (like 2)
static constexpr Float m = 6.646476410e-27; // He atomic mass in kg
static constexpr unsigned Pairwise_RandomizeRounds = 1500*N;
static constexpr Float SZa[3] = {4e-6,1e-7,1e-7};    // size of the box, when the particles are outside this box they reflect to the other side
static constexpr Float SZ(int i) {return SZa[i-1];}
static constexpr Float particleRadius = 140e-12; // He atomic radius
static constexpr Float particleRadius2 = (2*particleRadius)*(2*particleRadius);
static unsigned numCycles;
static unsigned cyclePrintPeriod; // how often to print the cycle
static constexpr Float Teff = Troom;  // effective temperature (same as default temperature)
static constexpr Float Vthermal = constexpr_funcs::sqrt((k*Teff)*2/m);  // thermal velocity
static constexpr Float Vcutoff = Vthermal*4.; // velocity over which there is a negligent number of particles XXX TODO 5. should be percentage estimate based
static constexpr Float penetrationCoefficient = 0.3; // fraction of the particle radius that we allow to be penetrated at worst, considering Vcutoff
static constexpr Float dt = particleRadius*penetrationCoefficient/Vcutoff;
static constexpr Float initE1 = constexpr_funcs::temperatureToEnergy(Troom)*0.95;
static constexpr Float initE2 = constexpr_funcs::temperatureToEnergy(Troom)*1.05;
static constexpr std::array<Float,2> energyLimits = {{0., 2.*initE2}}; // consider energies up to (..energyUpTo)
static constexpr unsigned NumOutputBuckets = 200; // how many buckets we build
static constexpr Float InteractionPct = 0.001; // how much energy is transferred
#if DBG_SAVE_IMAGES
static constexpr std::array<Float,2> imageAreaX = {{0.45,0.55}};
static constexpr std::array<Float,2> imageAreaY = {{0.45,0.55}};
static constexpr std::array<Float,2> imageAreaZ = {{0.45,0.55}}; //{{0.5-particleRadius,0.5+particleRadius}};// {{0.45,0.55}};
#endif

//
// Stats of the run
//
static unsigned statsNumBucketMoves = 0;
static unsigned statsNumCollisionsPP = 0;    // particle-particle collisions
static unsigned statsNumCollisionsPW = 0;

//
// Classes
//

class ParticlesIndex;

class Particle {
public: // data
  Vec3  pos;
  Vec3  v;
#if DBG_TRACK_PARTICLES
  unsigned pno;   // 1-based index
  bool     track;
#endif
public: // methods
  Particle() { }
  Particle(const Vec3 &newPos, const Vec3 &newV)
    : pos(newPos), v(newV) DBG_TRACK_PARTICLE_CODE(, track(false))
  {
  }
  ~Particle() { }
  Float energy() const {
    return m*(v*v)/2;
  }
  Float velocity() const {
    return v.len();
  }
  void move(Float dt) {
    auto collideWall = [](Float &partCoord, Float &partVelocity, Float wall1, Float wall2) {
      bool collided = false;
      if (partCoord < wall1) {
        partCoord = wall1 + (wall1-partCoord);
        partVelocity = -partVelocity;
        statsNumCollisionsPW++;
        collided = true;
      }
      if (partCoord > wall2) {
        partCoord = wall2 - (partCoord-wall2);
        partVelocity = -partVelocity;
        statsNumCollisionsPW++;
        collided = true;
      }
      return collided;
    };
    // move the point
    Vec3 pt(pos(X) + dt*v(X), pos(Y) + dt*v(Y), pos(Z) + dt*v(Z));
    // collide with walls
    bool collided = false;
    collided |= collideWall(pt(X), v(X), 0.+particleRadius,SZ(X)-particleRadius);
    collided |= collideWall(pt(Y), v(Y), 0.+particleRadius,SZ(Y)-particleRadius);
    collided |= collideWall(pt(Z), v(Z), 0.+particleRadius,SZ(Z)-particleRadius);
#if DBG_TRACK_PARTICLES
    if (track && collided)
      DBG_TRACK_PARTICLE_MSG(str(boost::format("particle#%1% collided with a wall") % pno)) // TODO report which wall
#endif
    // move to point
    move(pt);
  }
  void move(const Vec3 &newPos);
  static Float energyToVelocity(Float E) { // a reverse of energy()
    return std::sqrt(E*2/m);
  }
  void addEnergy(Float deltaE) {
    auto E1 = energy();
    auto E2 = E1 + deltaE;
    Float r = std::sqrt(E2/E1);
    v *= r;
  }
  static Float adjustCountPerUnit(Float V, unsigned cnt) {
    // count adjusted for unit size in the phase space
    return Float(cnt)/(4/3*V*V);
  }
  bool collisionCourse(const Particle &other) const {
    auto va = (v+other.v)/2;
    return (v-va)*(other.v-va) < 0; // when velocities are in opposite directions the scalar product is negative
  }
  Float distance2(const Particle &other) const {
    return (pos-other.pos).len2();
  }
  void collide(Particle &other) {
    //AOut() << "collide > this=" << *this << " other=" << other << "Etotal=" << (energy()+other.energy()) << std::endl;
    Vec3 dv = (other.v - v).project((other.pos - pos).normalizeZ());
    v += dv;
    other.v -= dv;
    //AOut() << "collide < this=" << *this << " other=" << other << "Etotal=" << (energy()+other.energy()) << std::endl;
#if DBG_TRACK_PARTICLES
    if (track)
      DBG_TRACK_PARTICLE_MSG(str(boost::format("particle#%1% collided with a particle #%2% (particle#%3% is the current particle in collision)") % pno % other.pno % pno))
    if (other.track)
      DBG_TRACK_PARTICLE_MSG(str(boost::format("particle#%1% collided with a particle #%2% (particle#%3% is the other particle in collision)") % other.pno % pno % other.pno))
#endif
  }
  friend std::ostream& operator<<(std::ostream &os, const Particle &p) {
    os << "{pos=" << p.pos << " v=" << p.v << "}";
    return os;
  }
}; // Particle

class ParticlesIndex { // index of particles in XYZ slot space
public:
  class Slot : public std::list<Particle*> {
    friend std::ostream& operator<<(std::ostream &os, const Slot &s) {
      unsigned ix = (&s - &instance.index[0][0][0])/(NSpaceSlots[1]*NSpaceSlots[2]);
      unsigned iy = (&s - &instance.index[ix][0][0])/NSpaceSlots[2];
      unsigned iz = (&s - &instance.index[ix][iy][0]);
      os << "{slot[" << ix << "][" << iy << "][" << iz << "] sz=" << s.size() << "}";
      return os;
    }
  }; // Slot
  static constexpr unsigned NSpaceSlots[3] = {constexpr_funcs::lim1((unsigned)(SZ(X)/constexpr_funcs::cbrt(SZ(X)*SZ(Y)*SZ(Z)/N*particlesPerBucket)+0.5)),
                                              constexpr_funcs::lim1((unsigned)(SZ(Y)/constexpr_funcs::cbrt(SZ(X)*SZ(Y)*SZ(Z)/N*particlesPerBucket)+0.5)),
                                              constexpr_funcs::lim1((unsigned)(SZ(Z)/constexpr_funcs::cbrt(SZ(X)*SZ(Y)*SZ(Z)/N*particlesPerBucket)+0.5))};
  static constexpr Float SZdivNSpaceSlots[3] = {SZa[0]/NSpaceSlots[0], SZa[1]/NSpaceSlots[1], SZa[2]/NSpaceSlots[2]};
  // leads to a fractional average occupancy of the bucket, seems to be the choice leading to the fastest computation
  static ParticlesIndex instance;
  static std::mutex instanceLock;
private:
  std::array<std::array<std::array<Slot,NSpaceSlots[2]>,NSpaceSlots[1]>,NSpaceSlots[0]> index;
public:
  ParticlesIndex() {
    std::cout << "ParticlesIndex: NSpaceSlots={" << NSpaceSlots[0] << "," << NSpaceSlots[1] << "," << NSpaceSlots[2] << "}"
              << " (total=" << formatUInt64(NSpaceSlots[0]*NSpaceSlots[1]*NSpaceSlots[2]) << ")" << std::endl;
  }
  auto& get() {return index;}
  void add(const Vec3 &pos, Particle *p) {
    add(findSlot(pos), p);
  }
  void add(Slot &slot, Particle *p) {
    slot.push_front(p);
  }
  void add(Particle *p) {
    add(findSlot(p), p);
  }
  void del(const Vec3 &pos, Particle *p) {
    del(findSlot(pos), p);
  }
  void del(Slot &slot, Particle *p) {
    for (auto i = slot.begin(); i != slot.end(); i++)
      if (*i == p) {
        slot.erase(i);
        return;
      }
    assert(0); // unreachable
  }
  void del(Particle *p) {
    del(findSlot(p), p);
  }
  Slot& findSlot(const Particle *p) {return findSlot(p->pos);}
  Slot& findSlot(const Vec3 &pos) {return index[coordToSlot(X, pos(X))][coordToSlot(Y, pos(Y))][coordToSlot(Z, pos(Z))];}
private:
  static unsigned coordToSlot(unsigned ci, Float c) {
    return c/SZdivNSpaceSlots[ci-1];
  }
}; // ParticlesIndex

ParticlesIndex ParticlesIndex::instance;
#if USE_PARALLELISM
std::mutex ParticlesIndex::instanceLock;
#endif

//
// Helper classes: iterators
//

namespace Iterate {

template<typename Fn>
static void CellNeighborsForward(int ix, int iy, int iz, Fn &&fn);

namespace Detail {

template<typename Fn>
static void ThroughOverlapsYZ(unsigned ix, Fn &&fn) {
  const auto& sx = ParticlesIndex::instance.get()[ix];
  for (unsigned iy = 0; iy < ParticlesIndex::NSpaceSlots[1]; iy++) {
    const auto& sy = sx[iy];
    for (unsigned iz = 0; iz < ParticlesIndex::NSpaceSlots[2]; iz++) {
      const auto& sz = sy[iz];
      // process all pairs: same bucket
      for (auto it1 = sz.begin(), ite = sz.end(); it1 != ite; it1++) {
        auto p1 = *it1;
        for (auto it2 = it1; ++it2 != ite;) {
          auto p2 = *it2;
#if DBG_TRACK_PARTICLES
          if (p1->track || p2->track)
            DBG_TRACK_PARTICLE_MSG(str(boost::format("checking same-bucket (%1%) for a collision between particle#%2% (%3%) and particle#%4% (%5%) ...")
                                                     % sz % p1->pno % *p1 % p2->pno % *p2))
#endif
          if (p1->distance2(*p2) <= particleRadius2)
            fn(p1, p2);
        }
      } // same bucket
      // process all pairs: cross-bucket
      CellNeighborsForward(ix,iy,iz, [&sz,fn](int ix, int iy, int iz) {
        auto &slot2 = ParticlesIndex::instance.get()[ix][iy][iz];
        for (auto it1 = sz.begin(), it1e = sz.end(); it1 != it1e; it1++) {
          auto p1 = *it1;
          for (auto it2 = slot2.begin(), it2e = slot2.end(); it2 != it2e; it2++) {
            auto p2 = *it2;
#if DBG_TRACK_PARTICLES
          if (p1->track || p2->track)
            DBG_TRACK_PARTICLE_MSG(str(boost::format("checking cross-bucket (%1%->%2%) for a collision between particle#%3% (%4%) and particle#%5% (%6%) ...")
                                                     % sz % slot2 % p1->pno % *p1 % p2->pno % *p2))
#endif
            if (p1->distance2(*p2) <= particleRadius2)
              fn(p1, p2);
          }
        }
      });
    }
  }
}

}; // Detail

template<typename Fn>
static void CellNeighborsForward(int ix, int iy, int iz, Fn &&fn) {
  for (auto [dix,diy,diz] : std::array<std::array<int,3>,13>( // all forward-oriented neighboring buckets: (3^3-1)/2 = 13
        {{{{ 0, 0,+1}}, {{ 0,+1,-1}}, {{ 0,+1, 0}}, {{ 0,+1,+1}},
          {{+1, 0,-1}}, {{+1, 0, 0}}, {{+1, 0,+1}}, {{+1,-1,-1}}, {{+1,-1, 0}}, {{+1,-1,+1}}, {{+1,+1,-1}}, {{+1,+1, 0}}, {{+1,+1,+1}}}}))
    if (/*0 <= ix+dix &&*/ ix+dix < (int)ParticlesIndex::NSpaceSlots[0] &&
        0 <= iy+diy &&     iy+diy < (int)ParticlesIndex::NSpaceSlots[1] &&
        0 <= iz+diz &&     iz+diz < (int)ParticlesIndex::NSpaceSlots[2])
      fn(ix+dix, iy+diy, iz+diz);
}

template<typename Fn>
static void ThroughOverlapsSer(Fn &&fn) {
  for (unsigned ix = 0; ix < ParticlesIndex::NSpaceSlots[0]; ix++)
    Detail::ThroughOverlapsYZ(ix, fn);
}

#if USE_PARALLELISM
template<typename Fn>
static void ThroughOverlapsPar(Fn &&fn) {
  unsigned dix = (ParticlesIndex::NSpaceSlots[0]+1)/NCPU;
  std::mutex lock;
  bool done[1+NCPU] = {false};
  auto fnStitchJob = [fn](unsigned xStitch) {
    //AOut() << TID << ">> stitch @ xStitch=" << xStitch << std::endl;
    Detail::ThroughOverlapsYZ(xStitch, fn);
    //AOut() << TID << "<< stitch @ xStitch=" << xStitch << std::endl;
  };
  for (unsigned ix = 0, cpu = 1; ix < ParticlesIndex::NSpaceSlots[0]; cpu++) {
    unsigned ix1 = ix + dix, ixe;
    if (cpu == NCPU && ix1 != ParticlesIndex::NSpaceSlots[0])
      ix1 = ParticlesIndex::NSpaceSlots[0];
    ixe = cpu < NCPU ? ix1-1 : ix1;
    threadPool->doJob([cpu, ix, ixe, fn, fnStitchJob, &done, &lock]() {
      //AOut() << TID << ">> {cpu#" << cpu << "} main @ ix=" << ix << ".." << ixe-1 << std::endl;
      for (unsigned x = ix; x < ixe; x++) 
        Detail::ThroughOverlapsYZ(x, fn);
      // done
      std::unique_lock<std::mutex> l(lock);
      done[cpu] = true;
      if (cpu > 1 && done[cpu-1])
        // stitch back
        threadPool->doJob([ix,fnStitchJob]() {
          fnStitchJob(ix-1);
        });
      if (cpu < NCPU && done[cpu+1])
        // stitch front
        threadPool->doJob([ixe,fnStitchJob]() {
          fnStitchJob(ixe);
        });
      //AOut() << TID << "<< main @ ix=" << ix << ".." << ixe-1 << std::endl;
    });
    ix = ix1;
  }
  // wait for all tasks to finish
  //AOut() << ">> waitForAll" << std::endl;
  threadPool->waitForAll();
  //AOut() << "<< waitForAll" << std::endl;
}
#endif

template<typename Fn>
static void ThroughOverlaps(Fn &&fn) {
#if USE_PARALLELISM
  if (2*NCPU-1 <= ParticlesIndex::NSpaceSlots[0]) // need to have at least 2 rows in each, except for the last one
    ThroughOverlapsPar(fn);
  else
    ThroughOverlapsSer(fn);
#else
  ThroughOverlapsSer(fn);
#endif
}

}; // Iterate

//
// Separate methods
//

void Particle::move(const Vec3 &newPos) {
  auto &slot1 = ParticlesIndex::instance.findSlot(this);
  pos = newPos;
  auto &slot2 = ParticlesIndex::instance.findSlot(this);
  if (&slot1 != &slot2) {
#if USE_PARALLELISM
    std::unique_lock<std::mutex> l(ParticlesIndex::instanceLock);
#endif
    ParticlesIndex::instance.del(slot1, this);
    ParticlesIndex::instance.add(slot2, this);
#if DBG_TRACK_PARTICLES
    if (track)
      DBG_TRACK_PARTICLE_MSG(str(boost::format("particle#%1% move between slots %2% -> %3%") % pno % slot1 % slot2)) // TODO find/keep slot numbers and print them
#endif
    //AOut() << "MOVE p=" << this << " x=" << x << " y=" << y << " z=" << z << std::endl;
    statsNumBucketMoves++;
  }
}

class DataBucket {
public:
  Float    V;    // velocity
  unsigned cnt;  // how many particles fall into this bucket
  DataBucket(Float newV) : V(newV), cnt(0) { }
}; // DataBucket

//
// Data
//

static std::array<Particle, Ncreate> particles;

//
// Random value generator
//

class Random {
  // types
  typedef std::uniform_real_distribution<Float> urdFloat;
  typedef std::uniform_int_distribution<unsigned> uidUnsigned;
  // fields
  unsigned                                 seed;
  std::mt19937                             generator;
  urdFloat     uniform01;
  uidUnsigned  uniform1N;
  urdFloat     uniformCoord[3]; // per dimension
  urdFloat     uniformEnergy;
public:
  Random()
  : seed(std::chrono::system_clock::now().time_since_epoch().count()),
    generator(seed),
    uniform01(0.0, 1.0),
    uniform1N(1,Ncreate),
    uniformCoord{urdFloat(0.0+particleRadius, SZ(X)-particleRadius),
                 urdFloat(0.0+particleRadius, SZ(Y)-particleRadius),
                 urdFloat(0.0+particleRadius, SZ(Z)-particleRadius)},
    uniformEnergy(initE1, initE2)
  { }
  auto rand01() {
    return uniform01(generator);
  }
  auto sphereAngles() {
    // from http://corysimon.github.io/articles/uniformdistn-on-sphere/
    Float theta = 2 * M_PI * rand01();
    Float phi = std::acos(1 - 2 * rand01());
    return std::array<Float,2>({{theta,phi}});
  }
  auto sphereCoords() {
    // from http://corysimon.github.io/articles/uniformdistn-on-sphere/
    auto [theta,phi] = sphereAngles();
    // convert to x/y/z
    double x = std::sin(phi)*std::cos(theta);
    double y = std::sin(phi)*std::sin(theta);
    double z = cos(phi);
    //
    return std::array<Float,3>({{x,y,z}});
  }
  auto coord(unsigned idx) {
    return uniformCoord[idx](generator);
  }
  auto energy() {
    return uniformEnergy(generator);
  }
  auto particlePair() {
    unsigned i1, i2;
    while ((i1 = uniform1N(generator)) == (i2 = uniform1N(generator))) { }
    return std::array<unsigned,2>({{i1,i2}});
  }
}; // Random

static Random rg;

static void generateParticles() {
  // TODO eliminate initial overlaps
  auto genParticleEnergyRange = []() {
    auto [sphX,sphY,sphZ] = rg.sphereCoords();
    auto E = rg.energy();
    auto V = Particle::energyToVelocity(E);
    return Particle(Vec3(rg.coord(0), rg.coord(1), rg.coord(2)), Vec3(sphX*V, sphY*V, sphZ*V));
  };
  for (unsigned i = 0; i < Ncreate; i++)
    particles[i] = genParticleEnergyRange();

#if DBG_TRACK_PARTICLES
  {
    //Float v = 0.001/*per-frame evolution percentage*/*0.08/*imageAreaX span*/ / 0.00001/*dt*/;
    //particles[0] = Particle(Vec3(0.5,0.5,0.5), Vec3(v,v,0));
    //particles[0].track = true;
    //particles[1] = Particle(Vec3(0.5+0.053+particleRadius,0.5+0.053,0.5), Vec3(-v,-v,0));
    //particles[2] = Particle(Vec3(-0.05+0.5,-0.05+0.5,0.5), Vec3(2*v,2*v,0));
  }
#endif

  // add to index
  for (auto &p : particles)
    ParticlesIndex::instance.add(p.pos, &p);

  // remove overlaps
  bool hadOverlaps;
  do {
    hadOverlaps = false;
    // find overlaps
    std::set<Particle*> overlaps;
    Iterate::ThroughOverlapsSer([&overlaps](Particle *p1, Particle *p2) {
      overlaps.insert(p2);
    });
    // regenerate collided particles
    hadOverlaps = !overlaps.empty();
    for (auto p : overlaps) {
      ParticlesIndex::instance.del(p);
      *p = genParticleEnergyRange();
      ParticlesIndex::instance.add(p);
    }
    AOut() << "generate.overlap.iteration: overlapsCount=" << overlaps.size() << std::endl;
  } while (hadOverlaps);
}

static std::vector<DataBucket> particlesToBuckets(const std::array<Particle,Ncreate> &particles) {
  Float minV = Particle::energyToVelocity(energyLimits[0]);
  Float maxV = Particle::energyToVelocity(energyLimits[1]);
  Float deltaV = (maxV - minV)/NumOutputBuckets;

  // generate initial buckets
  std::vector<DataBucket> buckets;
  for (unsigned i = 0; i < NumOutputBuckets; i++)
    buckets.push_back(DataBucket(minV + deltaV*i + deltaV/2));

  // count particles
  for (auto const &p : particles) {
    auto V = p.velocity();
    if (V < minV)
      continue;
    unsigned bucketNo = (V - minV)/deltaV;
    if (bucketNo >= NumOutputBuckets)
      continue; // discard the out-of-range energy
    buckets[bucketNo].cnt++;
  }

  return buckets;
}

template<typename Rec, typename FnPrn>
static void writeFile(const std::string &fname, const std::vector<Rec> &recs, FnPrn&& prn) {
  std::ofstream file;
  file.open(fname);
  for (auto const &r : recs)
    file << prn(r) << std::endl;
  file.close();
}

static void transferPercentageOfEnergy(Particle &p1, Particle &p2, Float frac) {
  // p1 transfers pct of enerty to p2
  //AOut() << "transferPercentageOfEnergy > p1.e=" << p1.energy() << " p2.e=" << p2.energy() << std::endl;
  auto deltaE = p1.energy()*frac;
  p1.addEnergy(-deltaE);
  p2.addEnergy(+deltaE);
  //AOut() << "transferPercentageOfEnergy < p1.e=" << p1.energy() << " p2.e=" << p2.energy() << std::endl;
}

static Float totalEnergy(std::array<Particle,Ncreate>::const_iterator it1, std::array<Particle,Ncreate>::const_iterator it2) {
  Float acc = 0.;
  auto it = it1;
  while (it != it2)
    acc += it++->energy();
  return acc / (it2 - it1);
}

//
// Images (save particles as an image for visual inspection)
//

#if DBG_SAVE_IMAGES

// tried CImg, but it is incompatible, because images having more than 1 bits per pixel are handled by calling external apps, like ImageMagick
#include <FreeImage.h>

class ImageSaver {
  // consts
  enum {SzX = 512, SzY = 512};
  // sub-area bounds that we image (we can't image all pixels because there's a million of them)
public:
  ImageSaver() { // have a constructor primarily for initialize/deinitialize functions
    FreeImage_Initialise();
    // delete previous images
    system("rm -f image-t=*");
  }
  ~ImageSaver() {
    FreeImage_DeInitialise();
  }
  void save(Float t, unsigned digitsInTime) const {
    // fname
    char fname[32];
    sprintf(fname, "image-t=%.*lf.png", digitsInTime, t);
    // write images
    FIBITMAP *bitmap = FreeImage_Allocate(SzX,SzY,3*8/*BitsPerPixel*/);
    RGBQUAD pixel;
    pixel = {255,255,255};
    FreeImage_FillBackground(bitmap, &pixel);
    for (auto &p : particles) {
      if (inBounds(p)) {
        auto [x,y,zc/*0..199*/] = mapToImage(p);
        //pixel = {255,BYTE(255-zc-56),BYTE(255-zc-56)}; // BGR
        pixel = {0,0,0}; // BGR
#if DBG_TRACK_PARTICLES
        if (p.track) {
          pixel = {0,0,255};
          //AOut() << "Display tracked particle @t=" << t << ": " << p << std::endl;
        }
#endif
        FreeImage_SetPixelColor(bitmap,x,y,&pixel);
      }
    }
    FreeImage_Save(FIF_PNG, bitmap, fname, 0);
    FreeImage_Unload(bitmap);
  }
private: // internals
  static bool inBounds(const Particle &p) {
    return inBounds(p.pos(X), imageAreaX) &&
           inBounds(p.pos(Y), imageAreaY) &&
           inBounds(p.pos(Z), imageAreaZ);
  }
  static bool inBounds(Float particleCoord, const std::array<Float,2> &bounds) {
    return bounds[0] <= particleCoord && particleCoord < bounds[1];
  }
  static std::tuple<unsigned/*X: 0..SzX-1*/,unsigned/*Y: 0..SzY-1*/,unsigned/*Z: 0..199*/> mapToImage(const Particle &p) {
    return {mapToImage(p.pos(X), imageAreaX, SzX),
            mapToImage(p.pos(Y), imageAreaY, SzY),
            mapToImage(p.pos(Z), imageAreaZ, 200)};
  }
  static unsigned mapToImage(Float particleCoord, const std::array<Float,2> &bounds, unsigned Sz) {
    return (particleCoord - bounds[0])/(bounds[1] - bounds[0])*Sz;
  }
}; // ImageSaver

static ImageSaver imageSaver;
#endif

//
// Serializer/unserializer
//

class Serializer {
  using json = nlohmann::json;
public:
  static auto serialize() {
    auto serializeVec = [](const Vec3 &v){
      json j = json::array();
      j.insert(j.end(), v(X));
      j.insert(j.end(), v(Y));
      j.insert(j.end(), v(Z));
      return j;
    };
    json arr = json::array();
    DBG_TRACK_PARTICLE_CODE(unsigned n = 1;)
    for (auto &p : particles) {
      json particle = json::object();
      DBG_TRACK_PARTICLE_CODE(particle["n"] = n++;)
      particle["x"] = serializeVec(p.pos);
      particle["v"] = serializeVec(p.v);
      arr.insert(arr.end(), particle);
    }
    return arr;
  }
  static void unserialize(const json &j) {
    auto unserializeVec = [](json &j){
      Vec3 v;
      v(X) = j[0];
      v(Y) = j[1];
      v(Z) = j[2];
      return v;
    };
    unsigned i = 0;
    for (auto &p : particles) {
      auto particle = j[i++];
      p = Particle(unserializeVec(particle["x"]), unserializeVec(particle["v"]));
    }
  }
}; // Serializer

//
// Evolve
//

class Evolver {
public:
  static void evolvePairwiseVelocity() {
    for (unsigned i = 0; i < Pairwise_RandomizeRounds; i++) {
      auto pp = rg.particlePair();
      transferPercentageOfEnergy(particles[pp[0]-1], particles[pp[1]-1], InteractionPct);
    }
  }
  static void evolvePhysically() {
    uint64_t cpuCycles0 = xasm::getCpuCycles();
#if DBG_SAVE_IMAGES
    imageSaver.save(0./*t*/, 5/*digits in time*/);
#endif
    // evolve
    Float t = 0.;
    for (unsigned cycle = 0; cycle < numCycles; cycle++) {
      //
      // move
      //

      moveIterateByParticle();
      //moveIterateByBucket();
      t += dt;

      //
      // iter-particle collisions
      //
      auto prevStatsNumCollisions = statsNumCollisionsPP;
      Iterate::ThroughOverlaps([](Particle *p1, Particle *p2) {
        assert(p1->collisionCourse(*p2)); // overshoot the center due to too high speed?
        p1->collide(*p2);
        statsNumCollisionsPP++;
      });
#if DBG_SAVE_IMAGES
      imageSaver.save(t, 5/*digits in time*/);
#endif
      // print tick and stats
      if ((cycle+1) % cyclePrintPeriod == 0) {
        uint64_t cpuCyclesNow = xasm::getCpuCycles();
        AOut() << "tick#" << cycle+1 << ":evolvePhysically:"
               << " avgCpuCyclesPerTick=" << formatUInt64((cpuCyclesNow - cpuCycles0)/uint64_t(cycle+1))
               << " statsNumCollisionsPP=" << statsNumCollisionsPP << " (+" << (statsNumCollisionsPP-prevStatsNumCollisions) << ")"
               << " statsNumCollisionsPW=" << statsNumCollisionsPW
               << " statsNumBucketMoves=" << statsNumBucketMoves
               << std::endl;
      }
    }
  }
private:
  static void moveIterateByParticle() { // faster when occupancy percentage is low
#if USE_PARALLELISM
    if (particles.size() >= NCPU)
      moveIterateByParticlePar();
    else
      moveIterateByParticleSeq();
#else
    moveIterateByParticleSeq();
#endif
  }
#if USE_PARALLELISM
  static void moveIterateByParticlePar() {
    unsigned ie = particles.size(), di = ie / NCPU;
    for (unsigned idx1 = 0, cpu = 1; idx1 < ie; cpu++) {
      unsigned idx2 = idx1 + di;
      if (idx2 > ie)
        idx2 = ie;
      if (cpu == NCPU && idx2 < ie)
        idx2 = ie;
      threadPool->doJob([idx1, idx2]() {
        for (unsigned i = idx1; i < idx2; i++)
          particles[i].move(dt);
      }); 
      idx1 = idx2;
    }
    // wait for all tasks to finish
    threadPool->waitForAll();
  }
#endif
  static void moveIterateByParticleSeq() {
    for (auto &p : particles)
      p.move(dt);
  }
  static void moveIterateByBucketHelper(ParticlesIndex::Slot::iterator it, const ParticlesIndex::Slot::iterator &ite) {
    auto p = *it;
    it++;
    if (it != ite)
      moveIterateByBucketHelper(it, ite);
    p->move(dt);
  }
  static void moveIterateByBucket() { // slower when occupancy percentage is low
    for (auto slot = &ParticlesIndex::instance.get()[0][0][0], slote = slot + ParticlesIndex::NSpaceSlots[0]*ParticlesIndex::NSpaceSlots[1]*ParticlesIndex::NSpaceSlots[2]; slot < slote; slot++)
      if (!slot->empty()) {
        moveIterateByBucketHelper(slot->begin(), slot->end());
        //sz += slot->size();
      }
  }
}; // Evolver

//
// MAIN
//

int mainGuarded(int argc, char *argv[]) {

  if (DBG_TRACK_PARTICLES) {
    std::cout << "!!!DEBUG CODE!!!" << std::endl;
    std::cout << "!!!DEBUG CODE!!! (DBG_TRACK_PARTICLES is enabled)" << std::endl;
    std::cout << "!!!DEBUG CODE!!!" << std::endl;
  }

#if USE_PARALLELISM
  threadPool = new ThreadPool(NCPU);
#endif

  //
  // cycles
  //

  auto cpuCycles0 = xasm::getCpuCycles();

  //
  // parse arguments
  //

#if DBG_TRACK_PARTICLES
  auto splitToUInt = [](std::string strToSplit, char delimeter) {
    std::stringstream ss(strToSplit);
    std::string item;
    std::vector<unsigned> split;
    while (std::getline(ss, item, delimeter))
       split.push_back(std::stoul(item));
    return split;
  };
#endif

  cxxopts::Options options("Particle Simulator", "Simulator of the gas particles motion");
  options.add_options()
    ("r,restart", "Restart using the snapshot of particle positions/velocities", cxxopts::value<std::string>())
    ("o,output",  "Write the final particle state into this file (default would be set to particles-{time}.json)", cxxopts::value<std::string>()->default_value(
      str(boost::format("particles-%1%.json") % std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count())))
    ("c,cycles",  "How many cycles to perform (default is 4000)", cxxopts::value<unsigned>()->default_value("4000"))
    ("p,print",  "How frequently to print cycles (default is 1, which means print every cycle)", cxxopts::value<unsigned>()->default_value("1"))
#if DBG_TRACK_PARTICLES
    ("t,track",   "Track particles: 1-based, colon-separated values", cxxopts::value<std::string>())
      // could use std::vector<unsigned> here to do '-t pno1 -t pno 2 ...', but maybe it's better (more compact) this way with colon-separated values
#endif
    ;

  auto result = options.parse(argc, argv);
  bool optRestart = result.count("restart") > 0;
  std::string optRestartFile = optRestart ? result["restart"].as<std::string>() : "";
  std::string optOutputFile = result["output"].as<std::string>();
  numCycles = result["cycles"].as<unsigned>();
  cyclePrintPeriod = result["print"].as<unsigned>();
#if DBG_TRACK_PARTICLES
  std::vector<unsigned> optTrack;
  if (result.count("track") > 0)
    optTrack = splitToUInt(result["track"].as<std::string>(), ':');
  { // assign pno in particles
    unsigned pno = 1;
    for (auto &p : particles)
      p.pno = pno++;
  }
#endif

  //
  // checks
  //
  for (auto c : {X,Y,Z})
    assert(2*particleRadius < Float(SZ(c))/ParticlesIndex::NSpaceSlots[c-1]); // particle diameter should be < index slot size, because otherwise collision detection
                                                                              // needs to look out more than in 1 slot away which makes it impractical
  //
  // params & stats
  //

  AOut() << "params(init): dt=" << dt
                      << " Vthermal=" << Vthermal
                      << " numCycles=" << numCycles
                      << std::endl;
  AOut() << "stats(init): spacePercentageOccupiedByParticles=" << Ncreate*(4./3.*M_PI*std::pow(particleRadius,3))/(SZ(X)*SZ(Y)*SZ(Z))*100. << "%"
                     << " avgParticlePerBucket=" << Float(Ncreate)/(ParticlesIndex::NSpaceSlots[0]*ParticlesIndex::NSpaceSlots[1]*ParticlesIndex::NSpaceSlots[2])
                     << " P/Patm (at Troom)=" << ((Ncreate/Na)*R*Troom/(SZ(X)*SZ(Y)*SZ(Z))/Patm)
                     << std::endl;

  //
  // generate or restart
  //
  if (optRestart) {
    std::cout << "unserializing from " << optRestartFile << " ..." << std::endl;
    std::ifstream file;
    file.open(optRestartFile, std::ifstream::in);
    Serializer::unserialize(nlohmann::json::parse(file));
    file.close();
    for (auto &p : particles)
      ParticlesIndex::instance.add(&p);
    std::cout << "done unserializing from " << optRestartFile << " ..." << std::endl;
  } else {
    generateParticles();
  }
#if DBG_TRACK_PARTICLES
  for (auto track : optTrack)
    if (track != 0)
      particles[track-1].track = true;
#endif

  //
  // initial log
  //
  std::cout << "log(init): energy-before=" << totalEnergy(particles.begin(), particles.end()) << " numCycles=" << numCycles << std::endl;

  //
  // evolve
  //
  if (0)
    Evolver::evolvePairwiseVelocity();
  if (1)
    Evolver::evolvePhysically();

  //
  // final log & stats
  //
  std::cout << "log(fini): energy-after=" << totalEnergy(particles.begin(), particles.end()) << std::endl;
  std::cout << "stats(fini): collisionsPerParticle(PP)=" << Float(statsNumCollisionsPP)/N
                        << " collisions(PP)=" << formatUInt64(statsNumCollisionsPP)
                        << " collisions(PW)=" << formatUInt64(statsNumCollisionsPW)
                        << std::endl;

  //
  // output
  //
  writeFile<DataBucket>("buckets.txt", particlesToBuckets(particles), [](const DataBucket &b) {
    std::ostringstream os;
    os << b.V << " " << b.cnt;
    return os.str();
  });

  //
  // serialize for later reuse
  //

  {
    std::cout << "serializing to " << optOutputFile << " ..." << std::endl;
    std::ofstream file;
    file.open(optOutputFile, std::ofstream::out);
    file << Serializer::serialize();
    file.close();
    std::cout << "done serializing to " << optOutputFile << " ..." << std::endl;
  }

  //
  // cycles
  //
  auto cpuCyclesEnd = xasm::getCpuCycles();
  std::cout << "consumed cycles: " << formatUInt64(cpuCyclesEnd-cpuCycles0) << " {now at " << formatUInt64(cpuCyclesEnd) << "}" << std::endl;

#if USE_PARALLELISM
  delete threadPool;
#endif

  return 0;
}

int main(int argc, char *argv[]) {
  try {
    return mainGuarded(argc, argv);
  } catch (cxxopts::OptionException e) {
    std::cerr << "Option parsing error: " << e.what() << std::endl;
    return 1;
  } catch (std::ifstream::failure e) {
    std::cerr << "File error: " << e.what() << std::endl;
    return 1;
  } catch (nlohmann::json::parse_error e) {
    std::cerr << "json error: " << e.what() << std::endl;
    return 1;
  } catch (std::exception const& e) {
    std::cerr << "Unknown exception of type '" << typeid(e).name() << "' caught" << std::endl;
    return 1;
  }
}
