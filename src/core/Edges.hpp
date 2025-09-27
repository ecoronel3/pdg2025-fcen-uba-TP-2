//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-05 16:33:49 taubin>
//------------------------------------------------------------------------
//
// Edges.hpp
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

#pragma once

#include <vector>

class Edges {

  // - the public interface to the Edges class only allows read-only
  //   access
  // - to get access to the _insertEdge method a subclass has to be
  //   created
  // - the Graph class is identical to the Edges class with the
  //   insertEdge method made public
  
public:

  // create a graph with nV vertices and no edges;
  // the range of valid vertex indices is 0<=iV<nV
  explicit Edges(int nV);

  // returns the number of vertices
  int getNumberOfVertices() const;

  // returns the number of edges nE at the time of the call;
  // at any particular time, the range of valid vertex indices is
  // 0<=iE<getNumberOfEdges()
  int getNumberOfEdges() const;

  // returns -1 if iV0==iV1 or one of the vertex indices is out of
  // range; also returns -1 if the edge (iV0,iV1) has not been
  // inserted into the Edges yet; otherwise it returns the edge index iE
  // assigned to the edge when inserted
  int getEdge(int iV0, int iV1) const;

  // an edge is stored internally as a pair of vertex indices
  // (iV0,iV1) so that iV0<iV1; getVertex0(iE) returns iV0, and
  // getVertex1(iE) returns iV1.
  int getVertex0(int iE) const;
  int getVertex1(int iE) const;

  // Edges Traversal sample code
  //
  // int nE = edges.getNumberOfEdges();
  // int iE,iV0,iV1;
  // for(iE=0;iE<nE;iE++) {
  //   iV0 = edges.getVertex0(iE);
  //   iV1 = edges.getVertex0(iE);
  //   // ...
  // }

protected:

  // remove all the edges, and change the number of vertices
  void reset(int nV);

  // - if iV0==iV1 or one of the two vertex indices is out of range,
  //   _insertEdge() returns -1 ;
  // - if iE=getEdge(iV0,iV1) is a valid edge index,
  //   _insertEdge() returns iE;
  // - otherwise a new edge index iE is assigned to the edge, and
  //   insertEdge() returns the new index iE
  int insertEdge(int iV0, int iV1);

private:

  // representation: array of single-linked lists

  // _first[iV0] is the index into _edge array corresponding to first
  // edge (iV0,iV1) so that iV0<iV1; _first[iV0]==-1 if the list is
  // empty
  std::vector<int> _first;
  // stores triples (iV0,iV1,next), where the next value is an index
  // into _edge array corresponding to next edge (iV0,iV1) so that
  // iV0<iV1; next==-1 indicates the end of the list; the order of the
  // triples in each list is not specified
  std::vector<int> _edge;

};
