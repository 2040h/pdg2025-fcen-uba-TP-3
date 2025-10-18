//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-05 23:13:01 taubin>
//------------------------------------------------------------------------
//
// HalfEdges.cpp (Assignment 3)
//
// Written by: Rundong He
//
// Software developed for the course
// Digital Geometry Processing
// Copyright (c) 2025, Gabriel Taubin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
//     * Redistributions of source code must retain the above
//       copyright notice, this list of conditions and the following
//       disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials
//       provided with the distribution.
//     * Neither the name of the Brown University nor the names of its
//       contributors may be used to endorse or promote products
//       derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GABRIEL
// TAUBIN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
// OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
// USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#include <math.h>
#include "HalfEdges.hpp"
// #include "Graph.hpp"

// TODO Mon Mar 6 2023
// - merge your code from Assignment 2

HalfEdges::HalfEdges(const int nVertices, const vector<int>&  coordIndex):
  Edges(nVertices), // a graph with no edges is created here
  _coordIndex(coordIndex),
  _twin(coordIndex.size(), -1),
  _face(coordIndex.size(), -1),
  _firstCornerEdge(),
  _cornerEdge()
{
  // TODO

  // - both the _twin array and the _face array should end up being of
  //   the same size as the _coordIndex array
  // - for each corner index iC contained in face iF
  // - the half edge src is iC, and the half edge dst is iC+1 if iC is
  //   not the last corner of the face; otherwise is iC0
  //
  //   if _coordIndex[iC]>=0 then
  //     _face[iC] should be equal to iF
  //     if 
  //     _twin[iC] should be equal to the corner index of the twin half edge 
  //   if _coordIndex[iC]<0 then
  //   _face[
  
  //          ⬆ ???
  
  int nV = nVertices;
  int nC = static_cast<int>(_coordIndex.size()); // number of corners

  // 0) just to be safe, verify that for each corner iC that
  //    -1<=iV && iV<nV, where iV=coordIndex[iC]
  //
  //    if you find a violation, you can increment and the variable
  //    nV, and then use the method Edges::_reset() to adjust the
  //    number of vertices of the graph, if necessary; or you can
  //    abort throwing an StrException

  for (int iC = 0; iC < nC; ++iC) {
    int iV = coordIndex[iC];
    if (-1 > iV or iV >= nV) _reset(iV + 1);
  }
  
  // 1) create an empty vector<int> to count the number of incident
  //    faces per edge; size is not known at this point because the
  //    edges have not been created yet
  vector<int> nFacesEdge;

  // 2) insert all the edges in the graph; at the same time initialize
  //    the _twin array so that all the half edges are boundary, count
  //    the number of incident faces per edge, fill the _face
  //    array, and count the number of faces incident to each edge
  int iV0,iV1,iF,iE,iC,iC0,iC1;
  for(iF=iC0=iC1=0;iC1<nC;iC1++) {
    if(_coordIndex[iC1]>=0) continue;
    // face iF comprises corners iC0<=iC<iC1
    // - each corner in this range corresponds to one half edge
    // - find the next corner within the face
    // - get the two vertex indices and insert an edge in the graph if
    //   not already there
    for (iC = iC0; iC < iC1; ++iC) {
      _face[iC] = iF;
      int x = iC;
      int y = (iC + 1 == iC1) ? iC0 : iC + 1;

      iV0 = _coordIndex[x];
      iV1 = _coordIndex[y];

      // - note that Edges::_insertEdge return the edge index number of
      //   a newly created edge, or the index of an extisting edge
      iE = _insertEdge(iV0,iV1); // Edges method
      // - note that iE might be >= nFacesEdge.size() at this point, and
      //   you may need to increase the size of nFacesEdge first
      // - ...
      if (iE >= nFacesEdge.size()) {
        if (iE >= nFacesEdge.capacity()) nFacesEdge.reserve(2 * nFacesEdge.capacity());
        nFacesEdge.resize(iE + 1, 0);
      }
      nFacesEdge[iE]++;
    }

    // increment variables to continue processing next face
    iC0 = iC1+1; iF++;
  }
  int nF = iF;

  int nE = getNumberOfEdges();
  
  // 3) create an array to hold the first twin corner for each edge
  vector<int> twinCorner(nE, -1);
  // - the size of this array should be equal to the number of edges
  // - initialize it with -1's

  // 4) fill the _twin array
  // - visit all the half edges using a loop similar to the one used in step 2)
  // - for each half edge iC, get the src and dst vertex indices, and
  //   from them the index iE of the corresponding edge
  // - if twinCorner[iE]<1 save iC in twinCorner[iE]
  //   otherwise save the value stored in twinCorner[iE] in _twin[iC]
  //   and iC in _twin[_twin[iC]]
  for(iF=iC0=iC1=0;iC1<nC;iC1++) {
    if(_coordIndex[iC1]>=0) continue;

    _twin[iC1] = iC1 - iC0;

    for (iC = iC0; iC < iC1; ++iC) {
      int src = iC;
      int dst = (iC + 1 == iC1) ? iC0 : iC + 1;

      iE = getEdge(src, dst);
      if (twinCorner[iE] < 0)
        twinCorner[iE] = iC;
      else {
        _twin[iC] = twinCorner[iE];
        _twin[_twin[iC]] = iC;
      }
    }

    iC0 = iC1+1; iF++;
  }

  // consistently oriented
  /* \                  / */
  /*  \ iC01 <-- iC00  /  */
  /*   X ---- iE ---- X   */
  /*  / iC10 --> iC11  \  */
  /* /                  \ */

  // oposite orientation
  /* \                  / */
  /*  \ iC01 --> iC00  /  */
  /*   X ---- iE ---- X   */
  /*  / iC10 --> iC11  \  */
  /* /                  \ */

  // a decision has to be made about inconsistently oriented half
  // edges incident to the same edge, as well as how to deal with
  // singular edges; for the moment let's assume that the mesh does
  // not have singular edges, but inconsistently oriented half edges
  // incident to the same edge are made twins (i.e. we do not have to
  // check for orientation here); later on we may want to modify this
  // class to have the option to do one thing or the other, and
  // methods to indicate which case we have.

  // get everything up to here implemented, debugged, and commited
  // before continuing

  // 5) initialize the array of arrays representing the half-edge to
  //    edge incident relationships
  //    _firstCornerEdge, and _cornerEdge
  //    - the size of _firstCornerEdge should be equal to nE+1
  //    - the size of _cornerEdge should be equal to the number of valid corners
  //      nC-nF
  //    - set boundaries
  //      _firstCornerEdge[0]=0
  //      _firstCornerEdge[iE+1] = _firstCornerEdge[iE]+nFacesEdge[iE] (1<=iE<nE)
  _firstCornerEdge = vector<int>(getNumberOfEdges());
  _cornerEdge = vector<int>(nC - nF, -1);

  _firstCornerEdge[0] = 0;
  for (iE = 1; iE < getNumberOfEdges(); ++iE)
    _firstCornerEdge[iE] = _firstCornerEdge[iE - 1] + nFacesEdge[iE];


  // 6) fill the array of arrays - the indices of corners incident to
  //    edge iE (1 if boundary, 2 if regular, >2 if singular) should
  //    be stored consecutively in _cornerEdge starting at the
  //    location _firstCornerEdge[iE]

  for(iF=iC0=iC1=0;iC1<nC;iC1++) {
    if(_coordIndex[iC1]>=0) continue;

    for (iC = iC0; iC < iC1; ++iC) {
      int src = iC;
      int dst = (iC + 1 == iC1) ? iC0 : iC + 1;

      iE = getEdge(src, dst);
      int* cell = &_cornerEdge[_firstCornerEdge[iE]];
      if (*cell == -1) *cell = iC;
      else *(cell + 1) = iC;
      // assuming no singular edge
    }

    iC0 = iC1+1; iF++;
  }


  _hasBoundaryEdges = _hasRegularEdges = _hasSingularEdges = false;

  for (int iE = 0; iE < getNumberOfEdges(); ++iE) {
    if (isBoundaryEdge(iE)) _hasBoundaryEdges = true;
    else if (isRegularEdge(iE)) _hasRegularEdges = true;
    else if (isSingularEdge(iE)) _hasSingularEdges = true;
  }
}

int HalfEdges::getNumberOfCorners() const {
  return static_cast<int>(_coordIndex.size());
}

// in all subsequent methods check that the arguments are valid, and
// return -1 if any argument is out of range
#define check(condition) if (not (condition)) return -1;
#define nC (_coordIndex.size())
#define iC_in_range (0 <= iC and iC < nC)
#define valid_iC (iC_in_range and (_coordIndex[iC] >= 0))
#define nE (getNumberOfEdges())
#define valid_iE (0 <= iE and iE < nE)

// half-edge method srcVertex()
int HalfEdges::getFace(const int iC) const {
  check(valid_iC);
  return _face[iC];
}

// half-edge method srcVertex()
int HalfEdges::getSrc(const int iC) const {
  check(valid_iC);
  return _coordIndex[iC];
}

// half-edge method dstVertex()
int HalfEdges::getDst(const int iC) const {
  check(valid_iC);
  return _coordIndex[getNext(iC)];
}

// half-edge method next()
int HalfEdges::getNext(const int iC) const {
  // TODO
  // if iC is the last corner of its face, use the face size
  // stored in _twin[iC+1] to locate the first corner of the face
  check(valid_iC);
  if (_coordIndex[iC + 1] >= 0) return iC + 1;

  int face_size = _twin[iC + 1];
  int jC = iC - face_size + 1;
  return jC;
}

// half-edge method prev()
int HalfEdges::getPrev(const int iC) const {
  // TODO
  
  // if iC is the first corner of its face, since the face size is
  // stored at the end of the face in the _twin array, you will have
  // to search forward for the last corner of the face; you can use
  // the fact that all the faces have at least 3 corners to start the
  // search for the face separator at iC+3
  check(valid_iC);
  if (iC - 1 >= 0 and _coordIndex[iC - 1] >= 0) return iC - 1;

  int jC = iC + 3;
  while (_coordIndex[jC] >= 0) jC++;
  return jC - 1;
}

int HalfEdges::getTwin(const int iC) const {
  check(valid_iC);
  return _twin[iC];
}

// represent the half edge as an array of lists, with one list
// associated with each edge

int HalfEdges::getNumberOfEdgeHalfEdges(const int iE) const {
  if (not valid_iE) return 0;
  if (iE + 1 == nE) return _firstCornerEdge.size() - _firstCornerEdge[iE];
  return _firstCornerEdge[iE + 1] - _firstCornerEdge[iE];
}

int HalfEdges::getEdgeHalfEdge(const int iE, const int j) const {
  check(valid_iE);
  if (0 > j or j >= getNumberOfEdgeHalfEdges(iE)) return -1;

  return _cornerEdge[_firstCornerEdge[iE] + j];
}


// TODO Mon Mar 6 2023
// - new functions to implement

bool HalfEdges::isOriented(const int iC) const {
  if (not valid_iC) return false;
  // iC     : iV00->iV01
  // iCtwin : iV10->iV11

  int iCtwin = getTwin(iC);
  int iV00 = getSrc(iC), iV01 = getDst(iC),
      iV10 = getSrc(iCtwin), iV11 = getDst(iCtwin);
  
  // if (twin half edges are consistently oriented)
  /* \                  / */
  /*  \ iV01 <-- iV00  /  */
  /*   X ------------ X   */
  /*  / iV10 --> iV11  \  */
  /* /                  \ */
  // return true;
  
  // else if (twin half edges are not consistently oriented)
  /* \                  / */
  /*  \ iV01 <-- iV00  /  */
  /*   X ------------ X   */
  /*  / iV11 <-- iV10  \  */
  /* /                  \ */
  // return false;

  if (iV00 == iV11 and iV01 == iV10) return true;
  else if (iV00 == iV10 and iV01 == iV11) return false;

  cout << "This line shouldn't be reached" << endl;
  return false;
}

// half-edge method getFaceSize()
int HalfEdges::getFaceSize(const int iC) const {
  check(valid_iC);
  
  int j = iC + 1;
  while (_coordIndex[j] >= 0) j++;
  return _twin[j];
}
  
int HalfEdges::getNumberOfFacesEdge(const int iE) const {
  check(valid_iE);

  int j = _firstCornerEdge[iE];
  if (iE + 1 < _firstCornerEdge.size())
    return _firstCornerEdge[iE + 1] - j;
  else 
    return _cornerEdge.size() - j;
}

bool HalfEdges::isBoundaryEdge(const int iE) const {
  return valid_iE and getNumberOfFacesEdge(iE) == 1;
}

bool HalfEdges::isRegularEdge(const int iE) const {
  return valid_iE and getNumberOfFacesEdge(iE) == 2;
}

bool HalfEdges::isSingularEdge(const int iE) const {
  return valid_iE and getNumberOfFacesEdge(iE) > 2;
}

// HINTS

// - it is best to determine the return values of these methods in the
//   constructor and save them in private variables

bool HalfEdges::hasBoundaryEdges() const {
  return _hasBoundaryEdges;
}

bool HalfEdges::hasRegularEdges() const {
  return _hasRegularEdges;
}

bool HalfEdges::hasSingularEdges() const {
  return _hasSingularEdges;
}

#undef check
#undef nC
#undef iC_in_range
#undef valid_iC
#undef nE 
#undef valid_iE