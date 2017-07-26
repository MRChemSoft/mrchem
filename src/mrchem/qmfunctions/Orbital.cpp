//#include <fstream>

#include "Orbital.h"
#include "FunctionTree.h"
#include "SerialFunctionTree.h"
#include "parallel.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace std;

Orbital workOrb(2, Orbital::Paired);

Orbital::Orbital(int occ, int s)
        : QMFunction(0, 0),
          spin(s),
          occupancy(occ),
          error(1.0) {
}

Orbital::Orbital(const Orbital &orb)
        : QMFunction(orb),
          spin(orb.spin),
          occupancy(orb.occupancy),
          error(1.0) {
}

void Orbital::clear(bool free) {
    clearReal(free);
    clearImag(free);
}

void Orbital::compare(const Orbital &orb) const {
    if (this->compareOccupancy(orb) < 0) {
        MSG_WARN("Different occupancy");
    }
    if (this->compareSpin(orb) < 0) {
        MSG_WARN("Different spin");
    }
}

int Orbital::compareOccupancy(const Orbital &orb) const {
    if (this->getOccupancy() == orb.getOccupancy()) {
        return this->getOccupancy();
    }
    return -1;
}

int Orbital::compareSpin(const Orbital &orb) const {
    if (this->getSpin() == orb.getSpin()) {
        return this->getSpin();
    }
    return -1;
}

double Orbital::getExchangeFactor(const Orbital &orb) const {
    if (orb.getSpin() == Orbital::Paired) {
        return 1.0;
    } else if (this->getSpin() == Orbital::Paired) {
        return 0.5;
    } else if (this->getSpin() == orb.getSpin()) {
        return 1.0;
    }
    return 0.0;
}

bool Orbital::isConverged(double prec) const {
    bool converged = false;
    if (prec < 0.0) converged = true;
    if (getError() < prec) converged = true;
    return converged;
}

complex<double> Orbital::dot(Orbital &ket) {
    Orbital &bra = *this;
    if ((bra.getSpin() == Orbital::Alpha) and (ket.getSpin() == Orbital::Beta)) return 0.0;
    if ((bra.getSpin() == Orbital::Beta) and (ket.getSpin() == Orbital::Alpha)) return 0.0;
    return QMFunction::dot(ket);
}

//send an orbital with MPI
void Orbital::send_Orbital(int dest, int tag){
#ifdef HAVE_MPI
  MPI_Status status;
  MPI_Comm comm=MPI_COMM_WORLD;

  struct Metadata{
    int spin;
    int occupancy;
    int NchunksReal;
    int NchunksImag;
    double error;
  };

  Metadata Orbinfo;

  Orbinfo.spin=this->getSpin();
  Orbinfo.occupancy=this->getOccupancy();
  Orbinfo.error=this->getError();
  if(this->hasReal()){
    Orbinfo.NchunksReal = this->real().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
  }else{Orbinfo.NchunksReal = 0;}
  if(this->hasImag()){
    Orbinfo.NchunksImag = this->imag().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
  }else{Orbinfo.NchunksImag = 0;}
  
  int count=sizeof(Metadata);
  MPI_Send(&Orbinfo, count, MPI_BYTE, dest, 0, comm);
  
  if(this->hasReal())Send_SerialTree(&this->real(), Orbinfo.NchunksReal, dest, tag, comm);
  if(this->hasImag())Send_SerialTree(&this->imag(), Orbinfo.NchunksImag, dest, tag*10000, comm);
  
#endif
}

//send an orbital with MPI
void Orbital::Isend_Orbital(int dest, int tag){
#ifdef HAVE_MPI
  MPI_Status status;
  MPI_Comm comm=MPI_COMM_WORLD;

  struct Metadata{
    int spin;
    int occupancy;
    int NchunksReal;
    int NchunksImag;
    double error;
  };

  Metadata Orbinfo;

  Orbinfo.spin=this->getSpin();
  Orbinfo.occupancy=this->getOccupancy();
  Orbinfo.error=this->getError();
  if(this->hasReal()){
    Orbinfo.NchunksReal = this->real().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
  }else{Orbinfo.NchunksReal = 0;}
  if(this->hasImag()){
    Orbinfo.NchunksImag = this->imag().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
  }else{Orbinfo.NchunksImag = 0;}
  
  MPI_Request request;
  int count=sizeof(Metadata);
  MPI_Isend(&Orbinfo, count, MPI_BYTE, dest, 0, comm, &request);
  
  if(this->hasReal())ISend_SerialTree(&this->real(), Orbinfo.NchunksReal, dest, tag, comm);
  if(this->hasImag())ISend_SerialTree(&this->imag(), Orbinfo.NchunksImag, dest, tag*10000, comm);
  
#endif
}
//receive an orbital with MPI
void Orbital::Rcv_Orbital(int source, int tag){
#ifdef HAVE_MPI
  MPI_Status status;
  MPI_Comm comm=MPI_COMM_WORLD;

  struct Metadata{
    int spin;
    int occupancy;
    int NchunksReal;
    int NchunksImag;
    double error;
  };

  Metadata Orbinfo;

  int count=sizeof(Metadata);
  MPI_Recv(&Orbinfo, count, MPI_BYTE, source, 0, comm, &status);
  this->setSpin(Orbinfo.spin);
  this->setOccupancy(Orbinfo.occupancy);
  this->setError(Orbinfo.error);
  
  if(Orbinfo.NchunksReal>0){
    if(not this->hasReal()){
      //We must have a tree defined for receiving nodes. Define one:
      this->allocReal();
    }
    Rcv_SerialTree(&this->real(), Orbinfo.NchunksReal, source, tag, comm);}

  if(Orbinfo.NchunksImag>0){
    if(not this->hasImag()){
      //We must have a tree defined for receiving nodes. Define one:
      this->allocImag();
    }
    Rcv_SerialTree(&this->imag(), Orbinfo.NchunksImag, source, tag*10000, comm);
  }else{
    //&(this->imag())=0;
  }
  
#endif

}

//receive an orbital with MPI
void Orbital::IRcv_Orbital(int source, int tag){
#ifdef HAVE_MPI
  MPI_Status status;
  MPI_Comm comm=MPI_COMM_WORLD;

  struct Metadata{
    int spin;
    int occupancy;
    int NchunksReal;
    int NchunksImag;
    double error;
  };

  Metadata Orbinfo;
  MPI_Request request;

  int count=sizeof(Metadata);
  MPI_Irecv(&Orbinfo, count, MPI_BYTE, source, 0, comm, &request);

  this->setSpin(Orbinfo.spin);
  this->setOccupancy(Orbinfo.occupancy);
  this->setError(Orbinfo.error);
  
  if(Orbinfo.NchunksReal>0){
    if(not this->hasReal()){
      //We must have a tree defined for receiving nodes. Define one:
      this->allocReal();
    }
    IRcv_SerialTree(&this->real(), Orbinfo.NchunksReal, source, tag, comm);}

  if(Orbinfo.NchunksImag>0){
    if(not this->hasImag()){
      //We must have a tree defined for receiving nodes. Define one:
      this->allocImag();
    }
    IRcv_SerialTree(&this->imag(), Orbinfo.NchunksImag, source, tag*10000, comm);
  }else{
    //&(this->imag())=0;
  }
  
#endif

}
