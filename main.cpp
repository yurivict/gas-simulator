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

#include "Vec3.h"

// types
typedef double Float;

#define DBG_SAVE_IMAGES 0    // save images representing evolution of one particluar area
#define DBG_TRACK_PARTICLE 0 // track one particle

namespace math_constexpr {

// borrowed from the answer in https://stackoverflow.com/questions/23781506/compile-time-computing-of-number-of-bits-needed-to-encode-n-different-states/23795363
constexpr unsigned floorlog2(unsigned x) {
  return x == 1 ? 0 : 1+floorlog2(x >> 1);
}

constexpr unsigned ceillog2(unsigned x) {
  return x == 1 ? 0 : floorlog2(x - 1) + 1;
}

}; // math_constexpr

//
// params
//
static constexpr unsigned N = 1000000;
static constexpr unsigned Ncreate = N; // (DEBUG) should be =N, but allows to override the number of particles for debugging
static constexpr unsigned Pairwise_RandomizeRounds = 1500*N;
static constexpr unsigned SZ = 1.;    // size of the box, when the particles are outside this box they reflect to the other side
static constexpr unsigned numCycles = 4000;
static constexpr Float dt = 0.00001;
static constexpr Float initE1 = 1.9;
static constexpr Float initE2 = 2.1;
static constexpr std::array<Float,2> energyLimits = {{0., 2.*initE2}}; // consider energies up to (..energyUpTo)
static constexpr unsigned NumOutputBuckets = 200; // how many buckets we build
static constexpr Float m = 1.; // nominal mass, it shouldn't matter here
static constexpr Float InteractionPct = 0.001; // how much energy is transferred
static constexpr Float particleRadius = 0.003; //SZ/ParticlesIndex::NSpaceSlots/15.; // XXX arbitrary coefficient
static constexpr Float particleRadius2 = (2*particleRadius)*(2*particleRadius);

static constexpr std::array<Float,2> imageAreaX = {{0.45,0.55}};
static constexpr std::array<Float,2> imageAreaY = {{0.45,0.55}};
static constexpr std::array<Float,2> imageAreaZ = {{0.45,0.55}}; //{{0.5-particleRadius,0.5+particleRadius}};// {{0.45,0.55}};

//
// Stats of the run
//
static unsigned statsNumBucketMoves = 0;
static unsigned statsNumCollisions = 0;

//
// Classes
//

class ParticlesIndex;

class Particle {
public: // data
  Vec3  pos;
  Vec3  v;
#if DBG_TRACK_PARTICLE
  bool track;
#endif
public: // methods
  Particle() { }
  Particle(const Vec3 &newPos, const Vec3 &newV)
    : pos(newPos), v(newV)
#if DBG_TRACK_PARTICLE
      , track(false)
#endif
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
    auto collideWall = [](Float &partCoord, Float &partVelocity, Float wall1, Float wall2){
      if (partCoord < wall1) {
        partCoord = wall1 + (wall1-partCoord);
        partVelocity = -partVelocity;
      }
      if (partCoord > wall2) {
        partCoord = wall2 - (partCoord-wall2);
        partVelocity = -partVelocity;
      }
    };
    // move the point
    Vec3 pt(pos(X) + dt*v(X), pos(Y) + dt*v(Y), pos(Z) + dt*v(Z));
    // collide with walls
    collideWall(pt(X), v(X), 0.+particleRadius,SZ-particleRadius);
    collideWall(pt(Y), v(Y), 0.+particleRadius,SZ-particleRadius);
    collideWall(pt(Z), v(Z), 0.+particleRadius,SZ-particleRadius);
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
    //std::cout << "collide > this=" << *this << " other=" << other << "Etotal=" << (energy()+other.energy()) << std::endl;
    Vec3 dv = (other.v - v).project((other.pos - pos).normalizeZ());
    v += dv;
    other.v -= dv;
    //std::cout << "collide < this=" << *this << " other=" << other << "Etotal=" << (energy()+other.energy()) << std::endl;
  }
  static Float clp(Float c) {
    while (c < 0.)
      c += SZ;
    while (c >= SZ)
      c -= SZ;
    return c;
  }
  friend std::ostream& operator<<(std::ostream &os, const Particle &p) {
    os << "{pos=" << p.pos << " v=" << p.v << "}";
    return os;
  }
}; // Particle

class ParticlesIndex { // index of particles in XYZ slot space
public:
  typedef std::list<Particle*> Slot;
  enum {NSpaceSlots = 1<<(math_constexpr::ceillog2(N)/3+1)}; // leads to a fractional average occupancy of the bucket, seems to be the choice leading to the fastest computation
private:
  std::array<std::array<std::array<Slot,NSpaceSlots>,NSpaceSlots>,NSpaceSlots> index;
public:
  ParticlesIndex() {
    std::cout << "ParticlesIndex: NSpaceSlots=" << NSpaceSlots << std::endl;
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
  Slot& findSlot(const Vec3 &pos) {return index[coordToSlot(pos(X))][coordToSlot(pos(Y))][coordToSlot(pos(Z))];}
private:
  static unsigned coordToSlot(Float c) {
    return c / (1./NSpaceSlots);
  }
}; // ParticlesIndex

static ParticlesIndex particlesIndex;

//
// Helper classes: iterators
//
template<unsigned Max, typename Fn>
static void IterateCellNeighborsForward(int ix, int iy, int iz, Fn &&fn) {
  for (auto [dix,diy,diz] : std::array<std::array<int,3>,13>( // all forward-oriented neighboring buckets: (3^3-1)/2 = 13
        {{{{ 0, 0,+1}}, {{ 0,+1,-1}}, {{ 0,+1, 0}}, {{ 0,+1,+1}},
          {{+1, 0,-1}}, {{+1, 0, 0}}, {{+1, 0,+1}}, {{+1,-1,-1}}, {{+1,-1, 0}}, {{+1,-1,+1}}, {{+1,+1,-1}}, {{+1,+1, 0}}, {{+1,+1,+1}}}}))
    if (0 <= ix+dix && ix+dix < Max && 0 <= iy+diy && iy+diy < Max && 0 <= iz+diz && iz+diz < Max)
      fn(ix+dix, iy+diy, iz+diz);
}

template<typename Fn>
static void IterateThroughOverlaps(Fn &&fn) {
  const unsigned NSpaceSlots = ParticlesIndex::NSpaceSlots;
  for (int ix = 0; ix < NSpaceSlots; ix++) {
    const auto& sx = particlesIndex.get()[ix];
    for (int iy = 0; iy < NSpaceSlots; iy++) {
      const auto& sy = sx[iy];
      for (int iz = 0; iz < NSpaceSlots; iz++) {
        const auto& sz = sy[iz];
        // process all pairs: same bucket
        for (auto it1 = sz.begin(), ite = sz.end(); it1 != ite; it1++) {
          auto p1 = *it1;
          for (auto it2 = it1; ++it2 != ite;) {
            auto p2 = *it2;
            if (p1->distance2(*p2) <= particleRadius2)
              fn(p1, p2);
          }
        } // same bucket
        // process all pairs: cross-bucket
        IterateCellNeighborsForward<NSpaceSlots>(ix,iy,iz, [&sz,fn](int ix, int iy, int iz) {
          auto &slot2 = particlesIndex.get()[ix][iy][iz];
          for (auto it1 = sz.begin(), it1e = sz.end(); it1 != it1e; it1++) {
            auto p1 = *it1;
            for (auto it2 = slot2.begin(), it2e = slot2.end(); it2 != it2e; it2++) {
              auto p2 = *it2;
              if (p1->distance2(*p2) <= particleRadius2)
                fn(p1, p2);
            }
          }
        });
      }
    }
  }
}

//
// Separate methods
//

void Particle::move(const Vec3 &newPos) {
  auto &slot1 = particlesIndex.findSlot(this);
  pos = newPos;
  auto &slot2 = particlesIndex.findSlot(this);
  if (&slot1 != &slot2) {
    particlesIndex.del(slot1, this);
    particlesIndex.add(slot2, this);
    //std::cout << "MOVE p=" << this << " x=" << x << " y=" << y << " z=" << z << std::endl;
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
  unsigned                                 seed;
  std::mt19937                             generator;
  std::uniform_real_distribution<Float>    uniform01;
  std::uniform_int_distribution<unsigned>  uniform1N;
  std::uniform_real_distribution<Float>    uniformCoord;
  std::uniform_real_distribution<Float>    uniformEnergy;
public:
  Random()
  : seed(std::chrono::system_clock::now().time_since_epoch().count()),
    generator(seed),
    uniform01(0.0, 1.0),
    uniform1N(1,Ncreate),
    uniformCoord(0.0+particleRadius, SZ-particleRadius),
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
  auto coord() {
    return uniformCoord(generator);
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
    return Particle(Vec3(rg.coord(), rg.coord(), rg.coord()), Vec3(sphX*V, sphY*V, sphZ*V));
  };
  for (int i = 0; i < Ncreate; i++)
    particles[i] = genParticleEnergyRange();

#if DBG_TRACK_PARTICLE
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
    particlesIndex.add(p.pos, &p);

  // remove overlaps
  bool hadOverlaps;
  do {
    hadOverlaps = false;
    // find overlaps
    std::set<Particle*> overlaps;
    IterateThroughOverlaps([&overlaps](Particle *p1, Particle *p2) {
      overlaps.insert(p2);
    });
    // regenerate collided particles
    hadOverlaps = !overlaps.empty();
    for (auto p : overlaps) {
      particlesIndex.del(p);
      *p = genParticleEnergyRange();
      particlesIndex.add(p);
    }
    std::cout << "generate.overlap.iteration: overlapsCount=" << overlaps.size() << std::endl;
  } while (hadOverlaps);
}

static std::vector<DataBucket> particlesToBuckets(const std::array<Particle,Ncreate> &particles) {
  Float minV = Particle::energyToVelocity(energyLimits[0]);
  Float maxV = Particle::energyToVelocity(energyLimits[1]);
  Float deltaV = (maxV - minV)/NumOutputBuckets;

  // generate initial buckets
  std::vector<DataBucket> buckets;
  for (int i = 0; i < NumOutputBuckets; i++)
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
  //std::cout << "transferPercentageOfEnergy > p1.e=" << p1.energy() << " p2.e=" << p2.energy() << std::endl;
  auto deltaE = p1.energy()*frac;
  p1.addEnergy(-deltaE);
  p2.addEnergy(+deltaE);
  //std::cout << "transferPercentageOfEnergy < p1.e=" << p1.energy() << " p2.e=" << p2.energy() << std::endl;
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
#if DBG_TRACK_PARTICLE
        if (p.track) {
          pixel = {0,0,255};
          //std::cout << "Display tracked particle @t=" << t << ": " << p << std::endl;
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
// Evolve
//

static void evolvePairwiseVelocity() {
  for (int i = 0; i < Pairwise_RandomizeRounds; i++) {
    auto pp = rg.particlePair();
    transferPercentageOfEnergy(particles[pp[0]-1], particles[pp[1]-1], InteractionPct);
  }
}

static void evolvePhysically() {
#if DBG_SAVE_IMAGES
  imageSaver.save(0./*t*/, 5/*digits in time*/);
#endif
  // evolve
  Float t = 0.;
  for (int cycle = 0; cycle < numCycles; cycle++) {
    // move
    for (auto &p : particles)
      p.move(dt);
    t += dt;
    // iter-particle collisions
    auto prevStatsNumCollisions = statsNumCollisions;
    IterateThroughOverlaps([](Particle *p1, Particle *p2) {
      assert(p1->collisionCourse(*p2)); // overshoot the center due to too high speed?
      p1->collide(*p2);
      statsNumCollisions++;
    });
#if DBG_SAVE_IMAGES
    imageSaver.save(t, 5/*digits in time*/);
#endif
    // print tick and stats
    std::cout << "tick: evolvePhysically: t=" << t << " statsNumCollisions=" << statsNumCollisions << " (+" << (statsNumCollisions-prevStatsNumCollisions) << ") statsNumBucketMoves=" << statsNumBucketMoves << std::endl;
  }
}

//
// MAIN
//

int main(int argc, const char *argv[]) {

  //
  // generate
  //
  generateParticles();

  //
  // initial log & stats
  //
  std::cout << "log(init): energy-before=" << totalEnergy(particles.begin(), particles.end()) << " dt=" << dt << " numCycles=" << numCycles << std::endl;
  std::cout << "stats(init): spacePercentageOccupiedByParticles=" << Ncreate*(4./3.*M_PI*std::pow(particleRadius,3))/(SZ*SZ*SZ)*100. << "%"
                        << " avgParticlePerBucket=" << Float(Ncreate)/std::pow(Float(ParticlesIndex::NSpaceSlots),3)
                        << std::endl;

  //
  // evolve
  //
  if (0)
    evolvePairwiseVelocity();
  if (1)
    evolvePhysically();

  //
  // final log & stats
  //
  std::cout << "log(fini): energy-after=" << totalEnergy(particles.begin(), particles.end()) << std::endl;
  std::cout << "stats(fini): collisionsPerParticle=" << Float(statsNumCollisions)/N << std::endl;

  //
  // output
  //
  writeFile<DataBucket>("buckets.txt", particlesToBuckets(particles), [](const DataBucket &b) {
    std::ostringstream os;
    os << b.V << " " << b.cnt;
    return os.str();
  });
}

