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
#include <csignal>
#include <exception>
#include <functional>
#include <stdexcept> // runtime_error

#include <cxxopts.hpp>
#include <nlohmann/json.hpp>
#include <boost/format.hpp>

#if USE_PARALLELISM
#include <thread>
#include "ThreadPool.h"
#endif

#include "Vec3.h"

//
// TODO
// * gravity
// * thermal flow
//

//
// types
//

typedef double Float;
typedef uint64_t Count;
typedef uint64_t CpuCycles;

//
// physics constants and laws
//

namespace PhysicsConsts {

static constexpr Float Na = 6.022e+23;  // Avogadro constant
static constexpr Float R = 8.3144598;   // 8.3144598(48) in kg m² s⁻² K−1 mol−1 , see https://en.wikipedia.org/wiki/Gas_constant
static constexpr Float k = R/Na;        // Boltzmann constant, see https://en.wikipedia.org/wiki/Boltzmann_constant
static constexpr Float Troom = 293.15;  // room temperature in Kelvins: 20°C=293.15K
static constexpr Float Patm = 101325;   // 101.325 kPa

}; // PhysicsConsts

namespace PhysicsFormulas {

using namespace PhysicsConsts;

static Float constexpr energyToTemperature(Float E) {return E/k/(3./2.);}
static Float constexpr temperatureToEnergy(Float T) {return k*T*(3./2.);}

}; // PhysicsFormulas

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
#define DBG_TRACK_PARTICLE_CODE(cmds ...) /*do nothing*/
#endif

#if USE_PARALLELISM
#define TID "[" << ThreadPool::tid() << "]"
#else
#define TID ""
#endif

template<typename Tsrc, typename Tdst> Tdst convertType(const Tsrc &v);
template<> std::string convertType<std::string,std::string>(const std::string &v) {return v;}
template<> uint64_t convertType<std::string,uint64_t>(const std::string &v) {return std::stoull(v);}
template<> Float convertType<std::string,Float>(const std::string &v) {return std::stod(v);}

template<typename T>
static std::vector<T> splitString(const std::string &strToSplit, char delimeter) {
  std::stringstream ss(strToSplit);
  std::string item;
  std::vector<T> split;
  while (std::getline(ss, item, delimeter))
     split.push_back(convertType<std::string,T>(item));
  return split;
}


//
// helper classes
//

class exception : public std::exception {
  std::runtime_error m;
public:
  exception(const std::string &newMsg) : m(newMsg) { } // we use exception directly so constructor is public
  const char* what() const noexcept override {
    return m.what();
  }
}; // exception

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

}; // constexpr_funcs

namespace xasm {

  static inline CpuCycles getCpuCycles() {
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((CpuCycles)hi << 32) | lo;
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
// E = ½mv² = 3/2 kT      // relation between the average kinetic energy and temperature

//
// params
//
static constexpr Count N = 1000000;
static constexpr Count Ncreate = N; // (DEBUG) should be =N, but allows to override the number of particles for debugging
static constexpr Float particlesPerBucket = 1.; // average number of particles per bucket. Performance drops when the number is >1 (like 2)
static constexpr Float m = 6.646476410e-27; // He atomic mass in kg
static constexpr unsigned Pairwise_RandomizeRounds = 1500*N;
static constexpr Float SZa[3] = {4e-6,1e-7,1e-7};    // size of the box, when the particles are outside this box they reflect to the other side
static constexpr Float SZ(int i) {return SZa[i-1];}
static constexpr Float particleRadius = 140e-12; // He atomic radius
static constexpr Float particleRadius2 = (2*particleRadius)*(2*particleRadius);
static constexpr Float Teff = PhysicsConsts::Troom;  // effective temperature (same as default temperature)
static constexpr Float Vthermal = constexpr_funcs::sqrt(PhysicsFormulas::temperatureToEnergy(Teff)*2/m);  // thermal velocity
static constexpr Float Vcutoff = Vthermal*4.; // velocity over which there is a neglible number of particles XXX TODO 5. should be percentage estimate based
static constexpr Float penetrationCoefficient = 0.3; // fraction of the particle radius that we allow to be penetrated at worst, considering Vcutoff
static constexpr Float dt = particleRadius*penetrationCoefficient/Vcutoff;
static constexpr unsigned NumOutputBuckets = 200; // how many buckets we build
static constexpr Float InteractionPct = 0.001; // how much energy is transferred
#if DBG_SAVE_IMAGES
static constexpr std::array<Float,2> imageAreaX = {{0.45,0.55}};
static constexpr std::array<Float,2> imageAreaY = {{0.45,0.55}};
static constexpr std::array<Float,2> imageAreaZ = {{0.45,0.55}}; //{{0.5-particleRadius,0.5+particleRadius}};// {{0.45,0.55}};
#endif

// XXX XXX Where does this belong?
static Float maxwellSigma(Float T) {return std::sqrt(PhysicsConsts::k*T/m);}

//
// Misc variables
//

static bool signalOccurred = false;

//
// Stats of the run
//
#if USE_PARALLELISM
static std::mutex statsLock;  // protects all but statsNumBucketMoves that is protected by ParticlesIndex::instanceLock
#endif
static Count statsNumBucketMoves = 0;
static Count statsNumCollisionsPP = 0;    // particle-particle collisions
static Count statsNumCollisionsPW = 0;
static Count statsNumCollisionsPWtemp = 0; // PW collisions with temperature adjustment
static Float statsWallDeltaEnergy = 0;     // how much energy is injected through the wall (see the 'temperature' option, negative value means that the energy is taken)

//
// Distributions
//
/*
namespace Distributions {
namespace Maxwell {

using namespace PhysicsConsts;

static Float distributionPerVelocity(Float T, Float V) {
  auto V2 = V*V;
  return std::pow(m/(2*M_PI*k*T), 3./2.)*(4*M_PI*V2)*std::exp(-m*V2/(2*k*T));
}

}; // Maxwell
}; // Distributions
*/

//
// Wall functors
//

static Float wallFunctorElastic(Float v) {return -v;}
static std::function<Float(Float)> wallFunctors[3][2];

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
    return m*v.len2()/2;
  }
  Float velocity() const {
    return v.len();
  }
  void move(Float dt) {
    auto collideWall = [](unsigned coord, Float &partCoord, Float &partVelocity, Float wall1, Float wall2) {
      bool collided = false;
      if (partCoord < wall1) {
        partCoord = wall1 + (wall1-partCoord);
        partVelocity = wallFunctors[coord-1][0](partVelocity);
        statsNumCollisionsPW++;
        collided = true;
      }
      if (partCoord > wall2) {
        partCoord = wall2 - (partCoord-wall2);
        partVelocity = wallFunctors[coord-1][1](partVelocity);
        statsNumCollisionsPW++;
        collided = true;
      }
      return collided;
    };
    // move the point
    Vec3 pt(pos(X) + dt*v(X), pos(Y) + dt*v(Y), pos(Z) + dt*v(Z));
    // collide with walls
    bool collided = false;
    collided |= collideWall(X, pt(X), v(X), 0.+particleRadius, SZ(X)-particleRadius);
    collided |= collideWall(Y, pt(Y), v(Y), 0.+particleRadius, SZ(Y)-particleRadius);
    collided |= collideWall(Z, pt(Z), v(Z), 0.+particleRadius, SZ(Z)-particleRadius);
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
#if USE_PARALLELISM
  static std::mutex instanceLock;
#endif
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
// History
//

class History {
  typedef std::chrono::time_point<std::chrono::system_clock> Tm;
  using json = nlohmann::json;
  class Record {
    std::string action;
    Tm tm1; // time it started
    Tm tm2; // time it ended
    std::map<std::string,std::string> props; // arbitrary properties of the run
    auto serialize() const {
      auto tmToStr = [](Tm tm) {
        std::time_t tm_c = std::chrono::system_clock::to_time_t(tm);
        std::ostringstream ss;
        ss << std::put_time(std::localtime(&tm_c), "%F %T");
        return ss.str();
      };
      json j = json::object();
      j["action"] = action;
      j["tm1"] = std::chrono::duration_cast<std::chrono::milliseconds>(tm1.time_since_epoch()).count();
      j["tm1-str"] = tmToStr(tm1);
      j["tm2"] = std::chrono::duration_cast<std::chrono::milliseconds>(tm2.time_since_epoch()).count();
      j["tm2-str"] = tmToStr(tm2);
      if (!props.empty()) {
        auto jprops = json::object();
        for (auto const &kv : props)
          jprops[kv.first] = kv.second;
        j["props"] = jprops;
      }
      return j;
    }
    static Record unserialize(const json &j) {
      Record rec;
      rec.action = j["action"];
      rec.tm1 = Tm(std::chrono::milliseconds(j["tm1"]));
      rec.tm2 = Tm(std::chrono::milliseconds(j["tm2"]));
      if (j.find("props") != j.end()) {
        auto jprops = j["props"];
        for (auto it = jprops.begin(), ite = jprops.end(); it != ite; it++)
          rec.props[it.key()] = it.value();
      }
      return rec;
    }
    friend class History;
  }; // Record
  std::vector<Record>   lst;
public:
  History() { }
  void begin(const std::string &action) {
    Record rec;
    rec.action = action;
    rec.tm1 = std::chrono::system_clock::now();
    lst.push_back(rec);
  }
  void addProperty(const std::string &key, const std::string &val) {
    lst.rbegin()->props[key] = val;
  }
  void end() {
    lst.rbegin()->tm2 = std::chrono::system_clock::now();
  }
  auto serialize() const {
    json j = json::array();
    for (auto rec : lst)
      j.insert(j.end(), rec.serialize());
    return j;
  }
  static History unserialize(const json &j) {
    History h;
    for (auto jr : j)
      h.lst.push_back(Record::unserialize(jr));
    return h;
  }
}; // History

static History history;

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
          //if (p1->track || p2->track)
          //  DBG_TRACK_PARTICLE_MSG(str(boost::format("checking same-bucket (%1%) for a collision between particle#%2% (%3%) and particle#%4% (%5%) ...")
          //                                           % sz % p1->pno % *p1 % p2->pno % *p2))
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
          //if (p1->track || p2->track)
          //  DBG_TRACK_PARTICLE_MSG(str(boost::format("checking cross-bucket (%1%->%2%) for a collision between particle#%3% (%4%) and particle#%5% (%6%) ...")
          //                                           % sz % slot2 % p1->pno % *p1 % p2->pno % *p2))
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
  Count    cnt;  // how many particles fall into this bucket
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
  typedef std::normal_distribution<Float> normFloat;
  // fields
  unsigned                                 seed;
  std::mt19937                             generator;
  urdFloat     uniformFloat01;
  uidUnsigned  uniform1N;
  uidUnsigned  uniformUInt01;
  urdFloat     uniformCoord[3]; // per dimension
  urdFloat     uniformEnergy;
  normFloat    normal;
public:
  Random()
  : seed(std::chrono::system_clock::now().time_since_epoch().count()),
    generator(seed),
    uniformFloat01(0.0, 1.0),
    uniform1N(1,Ncreate),
    uniformUInt01(0,1),
    uniformCoord{urdFloat(0.0+particleRadius, SZ(X)-particleRadius),
                 urdFloat(0.0+particleRadius, SZ(Y)-particleRadius),
                 urdFloat(0.0+particleRadius, SZ(Z)-particleRadius)},
    uniformEnergy(PhysicsFormulas::temperatureToEnergy(Teff)*0.95, PhysicsFormulas::temperatureToEnergy(Teff)*1.05)
  { }
  auto rand01() {
    return uniformFloat01(generator);
  }
  auto boolean() {
    return uniformUInt01(generator) == 1 ? true : false;
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
  auto energyUniformArea() {
    // one area
    return uniformEnergy(generator);
  }
  auto energySpike50_150() {
    // more dissipated pattern
    return PhysicsFormulas::temperatureToEnergy(Teff)*(boolean() ? 0.5 : 1.5);
  }
  auto particlePair() {
    unsigned i1, i2;
    while ((i1 = uniform1N(generator)) == (i2 = uniform1N(generator))) { }
    return std::array<unsigned,2>({{i1,i2}});
  }
  Vec3 maxwell3(Float sigma) {
    return Vec3(sigma*normal(generator), sigma*normal(generator), sigma*normal(generator));
  }
  Float maxwell1Plus(Float sigma) {
    auto v = sigma*normal(generator);
    return v >= 0 ? v : -v;
  }
}; // Random

static Random rg;

template<typename VGen>
static void generateParticles(VGen &&vgen) {
  // record in history
  history.begin(str(boost::format("generate %1% particles") % Ncreate));

  // TODO eliminate initial overlaps
  auto genParticleEnergyRange = [vgen]() {
    return Particle(Vec3(rg.coord(0), rg.coord(1), rg.coord(2)), vgen());
  };
  for (Count i = 0; i < Ncreate; i++)
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

  // record in history
  history.end();
}

static std::vector<DataBucket> particlesToBuckets(const std::array<Particle,Ncreate> &particles) {
  std::array<Float,2> energyLimits = {{0., 6.*PhysicsFormulas::temperatureToEnergy(Teff)}}; // TODO replace the ad-hoc 6. with the analytical expression based on the percentage of ignored particles
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
      continue; // discard out-of-range particles
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

namespace ParticlesTotal {

static Float energy() {
  Float e = 0.;
  for (auto &p : particles)
    e += p.energy();
  return e;
}

static Float temperature() {
  return PhysicsFormulas::energyToTemperature(energy()/Ncreate);
}

}; // ParticlesTotal

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
    json j = json::object();
    j["particles"] = serializeParticles();
    j["history"] = history.serialize();
    return j;
  }
  static void unserialize(const json &j) {
    unserializeParticles(j["particles"]);
    history = History::unserialize(j["history"]);
  }
private:
  static json serializeParticles() {
    auto serializeVec = [](const Vec3 &v){
      json j = json::array();
      j.insert(j.end(), v(X));
      j.insert(j.end(), v(Y));
      j.insert(j.end(), v(Z));
      return j;
    };

    json arrParticles = json::array();
    DBG_TRACK_PARTICLE_CODE(unsigned n = 1;)
    for (auto &p : particles) {
      json particle = json::object();
      DBG_TRACK_PARTICLE_CODE(particle["n"] = n++;)
      particle["x"] = serializeVec(p.pos);
      particle["v"] = serializeVec(p.v);
      arrParticles.insert(arrParticles.end(), particle);
    }

    return arrParticles;
  }
  static void unserializeParticles(const json &j) {
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
#if DBG_TRACK_PARTICLES
      p.track = false; // not serialized
#endif
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
  template<typename FnBefore, typename FnAfter>
  static void evolvePhysically(Count numCycles, FnBefore &&fnBeforeCycle, FnAfter &&fnAfterCycle) {
    // record in history
    history.begin(str(boost::format("evolve %1% particles, numCycles=%2%") % Ncreate % numCycles));
    history.addProperty("temperatureBefore",  str(boost::format("%1%") % ParticlesTotal::temperature()));
    // evolve
#if DBG_SAVE_IMAGES
    imageSaver.save(0./*t*/, 5/*digits in time*/);
#endif
    // evolve
    Float t = 0.;
    for (Count cycle = 1; cycle <= numCycles; cycle++) {
      fnBeforeCycle(cycle);

      //
      // move
      //

      moveIterateByParticle();
      //moveIterateByBucket();
      t += dt;

      //
      // iter-particle collisions
      //
      Iterate::ThroughOverlaps([](Particle *p1, Particle *p2) {
        assert(p1->collisionCourse(*p2)); // overshoot the center due to too high speed?
        p1->collide(*p2);
        statsNumCollisionsPP++;
      });
#if DBG_SAVE_IMAGES
      imageSaver.save(t, 5/*digits in time*/);
#endif
      // cycle is done: allow the caller to perform custom actions
      fnAfterCycle(cycle);
    }
    // record in history
    history.addProperty("numCycles",         str(boost::format("%1%") % numCycles));
    history.addProperty("dt",                str(boost::format("%1%") % dt));
    history.addProperty("temperatureAfter",  str(boost::format("%1%") % ParticlesTotal::temperature()));
    history.end();
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
    Count ie = particles.size(), di = ie / NCPU;
    for (Count idx1 = 0, cpu = 1; idx1 < ie; cpu++) {
      Count idx2 = idx1 + di;
      if (idx2 > ie)
        idx2 = ie;
      if (cpu == NCPU && idx2 < ie)
        idx2 = ie;
      threadPool->doJob([idx1,idx2]() {
        for (Count i = idx1; i < idx2; i++)
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
  //static void moveIterateByBucketHelper(ParticlesIndex::Slot::iterator it, const ParticlesIndex::Slot::iterator &ite) {
  //  auto p = *it;
  //  it++;
  //  if (it != ite)
  //    moveIterateByBucketHelper(it, ite);
  //  p->move(dt, fnWallX1, fnWallX2, fnWallY1, fnWallY2, fnWallZ1, fnWallZ2);
  //}
  //static void moveIterateByBucket() { // slower when occupancy percentage is low
  //  for (auto slot = &ParticlesIndex::instance.get()[0][0][0], slote = slot + ParticlesIndex::NSpaceSlots[0]*ParticlesIndex::NSpaceSlots[1]*ParticlesIndex::NSpaceSlots[2]; slot < slote; slot++)
  //    if (!slot->empty()) {
  //      moveIterateByBucketHelper(slot->begin(), slot->end());
  //      //sz += slot->size();
  //    }
  //}
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
  // local functions
  //

  auto fnSaveOutput = [](const std::string &fname) {
    AOut() << "serializing to " << fname << " ..." << std::endl;
    std::ofstream file;
    file.open(fname, std::ofstream::out);
    file << Serializer::serialize();
    file.close();
    AOut() << "done serializing to " << fname << " ..." << std::endl;
  };
  auto fnSaveBuckets = [](Count cycle) {
    writeFile<DataBucket>(str(boost::format("buckets-cycle%1%.txt") % cycle), particlesToBuckets(particles), [](const DataBucket &b) {
      std::ostringstream os;
      os << b.V << " " << b.cnt;
      return os.str();
    });
  };

  //
  // cycles
  //

  auto cpuCycles0 = xasm::getCpuCycles();

  //
  // parse arguments
  //

  cxxopts::Options options("Particle Simulator", "Simulator of the gas particles motion");
  options.add_options()
    ("d,distribution", "Velocity distribution used to generate particles, and to inject energy (maxwell,spike-50-150,uniform-area)", cxxopts::value<std::string>()->default_value("maxwell"))
    ("r,restart",      "Restart using the snapshot of particle positions/velocities", cxxopts::value<std::string>())
    ("o,output",       "Write the final particle state into this file (default would be set to particles-{time}.json)", cxxopts::value<std::string>()->default_value(
      str(boost::format("particles-%1%.json") % std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count())))
    ("c,cycles",       "How many cycles to perform (default is 4000)", cxxopts::value<Count>()->default_value("4000"))
    ("p,print",        "How frequently to print cycles (default is 1, which means print every cycle)", cxxopts::value<unsigned>()->default_value("1"))
    ("b,buckets",      "How frequently to print buckets (default is 0, which means print only in the end)", cxxopts::value<unsigned>()->default_value("0"))
    ("T,temperature",  "Set temperature of the walls, ex. -X:271:+X:272", cxxopts::value<std::string>()->default_value(""))
#if DBG_TRACK_PARTICLES
    ("t,track",        "Track particles: 1-based, colon-separated values", cxxopts::value<std::string>())
      // could use std::vector<unsigned> here to do '-t pno1 -t pno 2 ...', but maybe it's better (more compact) this way with colon-separated values
#endif
    ;

  auto result = options.parse(argc, argv);
  std::string optDistribution = result["distribution"].as<std::string>();
  if (optDistribution!="maxwell" && optDistribution!="spike-50-150" && optDistribution!="uniform-area")
    throw exception(str(boost::format("unknown velocity distribution '%1%' supplied to --generate argument") % optDistribution));
  bool optRestart = result.count("restart") > 0;
  std::string optRestartFile = optRestart ? result["restart"].as<std::string>() : "";
  std::string optOutputFile = result["output"].as<std::string>();
  Count numCycles = result["cycles"].as<Count>();
  unsigned cyclePrintPeriod = result["print"].as<unsigned>();
  unsigned cycleWriteBuckets = result["buckets"].as<unsigned>();
  auto optWallTemperatures = splitString<std::string>(result["temperature"].as<std::string>(), ':');
  if (optWallTemperatures.size()%2 != 0)
    throw exception(str(boost::format("wall temperature specifier is expected to consist of pairs, supplied '%1%' parameters") % optWallTemperatures.size()));
  
  Float optWallTempSigma[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
  for (unsigned i = 0, ie = optWallTemperatures.size(); i!=ie; i++) {
    auto val = optWallTemperatures[i];
    if (val.size() != 2 || (val[0]!='-' && val[0]!='+') || (val[1]!='X' && val[1]!='Y' && val[1]!='Z'))
      throw exception(str(boost::format("wall temperature specifier is expected to have {sign}{Coord}, supplied '%1%'") % val));
    switch (val[1]) {
    case 'X': *(val[0]=='-' ? &optWallTempSigma[X-1][0] : &optWallTempSigma[X-1][1]) = maxwellSigma(convertType<std::string,Float>(optWallTemperatures[++i])); break;
    case 'Y': *(val[0]=='-' ? &optWallTempSigma[Y-1][0] : &optWallTempSigma[Y-1][1]) = maxwellSigma(convertType<std::string,Float>(optWallTemperatures[++i])); break;
    case 'Z': *(val[0]=='-' ? &optWallTempSigma[Z-1][0] : &optWallTempSigma[Z-1][1]) = maxwellSigma(convertType<std::string,Float>(optWallTemperatures[++i])); break;
    }
  }
#if DBG_TRACK_PARTICLES
  std::vector<Count> optTrack;
  if (result.count("track") > 0)
    optTrack = splitString<Count>(result["track"].as<std::string>(), ':');
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
                     << " P/Patm (at Teff)=" << ((Ncreate/PhysicsConsts::Na)*PhysicsConsts::R*Teff/(SZ(X)*SZ(Y)*SZ(Z))/PhysicsConsts::Patm)
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
    if (optDistribution == "maxwell") {
      Float sigma = maxwellSigma(Teff);
      generateParticles([sigma]() {
        return rg.maxwell3(sigma);
      });
    } else if (optDistribution == "spike-50-150") {
      generateParticles([]() {
        auto [sphX,sphY,sphZ] = rg.sphereCoords();
        auto E = rg.energySpike50_150();
        auto V = Particle::energyToVelocity(E);
        return Vec3(sphX*V, sphY*V, sphZ*V);
      });
    } else if (optDistribution == "area") {
      generateParticles([]() {
        auto [sphX,sphY,sphZ] = rg.sphereCoords();
        auto E = rg.energyUniformArea();
        auto V = Particle::energyToVelocity(E);
        return Vec3(sphX*V, sphY*V, sphZ*V);
      });
    }
  }
#if DBG_TRACK_PARTICLES
  for (auto track : optTrack)
    if (track != 0)
      particles[track-1].track = true;
#endif

  //
  // initial log
  //
  {
    auto energy = ParticlesTotal::energy();
    std::cout << "log(init): energy(before)=" << energy
                        << " temperature(before)=" << PhysicsFormulas::energyToTemperature(energy/Ncreate)
                        << " numCycles=" << numCycles
                        << std::endl;
  }

  //
  // signal processing
  //

  std::signal(SIGHUP, [](int signum) { signalOccurred = true; });

  //
  // evolve
  //
  if (0)
    Evolver::evolvePairwiseVelocity();
  if (1) {
    unsigned prevStatsNumCollisionsPP = 0;
    auto cpuCyclesBefore = xasm::getCpuCycles();
    // fill wall collision functors
    for (auto c : std::array<unsigned,3>({{0,1,2}})) {
      if (optWallTempSigma[c][0]==0.)
        wallFunctors[c][0] = &wallFunctorElastic;
      else
        wallFunctors[c][0] = [=](Float v) {
          auto newV = rg.maxwell1Plus(optWallTempSigma[c][0]);
          Float deltaEnergy = m*(newV*newV - v*v)/2;
          {
#if USE_PARALLELISM
            std::unique_lock<std::mutex> l(statsLock);
#endif
            statsNumCollisionsPWtemp++;
            statsWallDeltaEnergy += deltaEnergy;
          }
          return +newV;
        };
      if (optWallTempSigma[c][1]==0.)
        wallFunctors[c][1] = &wallFunctorElastic;
      else
        wallFunctors[c][1] = [=](Float v) {
          auto newV = rg.maxwell1Plus(optWallTempSigma[c][1]);
          Float deltaEnergy = m*(newV*newV - v*v)/2;
          {
#if USE_PARALLELISM
            std::unique_lock<std::mutex> l(statsLock);
#endif
            statsNumCollisionsPWtemp++;
            statsWallDeltaEnergy += deltaEnergy;
          }
          return -newV;
        };
    }
    // evolve
    Evolver::evolvePhysically(numCycles, [&prevStatsNumCollisionsPP](unsigned cycle) {
        prevStatsNumCollisionsPP = statsNumCollisionsPP;
      }, [numCycles,cyclePrintPeriod,cycleWriteBuckets,&optOutputFile,&prevStatsNumCollisionsPP,&cpuCyclesBefore,fnSaveOutput,fnSaveBuckets](unsigned cycle) {
        // print tick and stats
        if (cycle%cyclePrintPeriod == 0) {
          Count cpuAtCycle = xasm::getCpuCycles();
          AOut() << "tick#" << cycle << ":evolvePhysically:"
                 << " avgCpuCyclesPerTick=" << formatUInt64((cpuAtCycle - cpuCyclesBefore)/Count(cycle))
                 << " statsNumCollisionsPP=" << statsNumCollisionsPP << " (+" << (statsNumCollisionsPP-prevStatsNumCollisionsPP) << ")"
                 << " statsNumCollisionsPW=" << statsNumCollisionsPW
                 << (statsNumCollisionsPWtemp>0 ? str(boost::format(" statsNumCollisionsPWtemp=%1%") % statsNumCollisionsPWtemp) : "")
                 << (statsWallDeltaEnergy!=0. ? str(boost::format(" statsWallDeltaEnergy=%1% (ΔT=%2%)")
                                                                  % statsWallDeltaEnergy
                                                                  % PhysicsFormulas::energyToTemperature(statsWallDeltaEnergy/Ncreate))
                                              : "")
                 << " statsNumBucketMoves=" << statsNumBucketMoves
                 << std::endl;
        }

        // save output if requested
        bool savedOnSignal = false;
        if (signalOccurred && cycle != numCycles) {
          fnSaveOutput(str(boost::format("%1%.cycle%2%.json") % optOutputFile % cycle));
          fnSaveBuckets(cycle);
          signalOccurred = false;
          savedOnSignal = true;
        }

        // write buckets
        if (!savedOnSignal && ((cycleWriteBuckets!=0 && cycle%cycleWriteBuckets==0) || cycle == numCycles))
          fnSaveBuckets(cycle);
      }
    );
  }

  //
  // final log & stats
  //
  {
    auto energy = ParticlesTotal::energy();
    std::cout << "log(fini): energy(after)=" << energy << std::endl;
    std::cout << "stats(fini): collisionsPerParticle(PP)=" << Float(statsNumCollisionsPP)/N
                          << " temperature(after)=" << PhysicsFormulas::energyToTemperature(energy/Ncreate)
                          << " collisions(PP)=" << formatUInt64(statsNumCollisionsPP)
                          << " collisions(PW)=" << formatUInt64(statsNumCollisionsPW)
                          << std::endl;
  }

  //
  // serialize for later reuse
  //

  fnSaveOutput(optOutputFile);
  if (numCycles == 0)
    fnSaveBuckets(0);

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
  } catch (cxxopts::OptionException const& e) {
    std::cerr << "option parsing error: " << e.what() << std::endl;
    return 1;
  } catch (std::ifstream::failure const& e) {
    std::cerr << "file error: " << e.what() << std::endl;
    return 1;
  } catch (nlohmann::json::exception const& e) {
    std::cerr << "json error: " << e.what() << std::endl;
    return 1;
  } catch (exception const& e) {
    std::cerr << "application error: " << e.what() << std::endl;
    return 1;
  } catch (std::exception const& e) {
    std::cerr << "unknown exception of type '" << typeid(e).name() << "' caught: " << e.what() << std::endl;
    return 1;
  }
}
