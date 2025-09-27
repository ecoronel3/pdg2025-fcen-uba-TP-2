//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-05 16:34:29 taubin>
//------------------------------------------------------------------------
//
// PolygonMesh.cpp
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

#include "PolygonMesh.hpp"

#include "Partition.hpp"

PolygonMesh::PolygonMesh(const int nVertices, const std::vector<int>& coordIndex):
  HalfEdges(nVertices,coordIndex),
  _nPartsVertex(),
  _isBoundaryVertex(),
  _numberOfFaces(0)
{

  const int nC = getNumberOfCorners();
  for (int iC = 0; iC < nC; ++iC) {
    if(_coordIndex[iC]>=0)
      continue;
    ++_numberOfFaces;
  }

  int nV = getNumberOfVertices();
  int nE = getNumberOfEdges(); // Edges method
  // int nF = getNumberOfFaces();


  // 1) classify the vertices as boundary or internal
  _isBoundaryVertex.resize(nV, false);
  for (int iC = 0; iC < nC; ++iC) {
    // - for edge boundary iE label its two end vertices as boundary
    if (coordIndex[iC] < 0)
      continue;

    const int iV0 = getSrc(iC);
    const int iV1 = getDst(iC);
    const int iE = getEdge(iV0, iV1);
    if (getNumberOfEdgeHalfEdges(iE) == 1) {
      _isBoundaryVertex[iV0] = true;
      _isBoundaryVertex[iV1] = true;
    }
  }

  // 2) create a partition of the corners in the stack
  Partition partition(nC);
  // 3) for each regular edge
  //    - get the two half edges incident to the edge
  //    - join the two pairs of corresponding corners across the edge
  //    - you need to take into account the relative orientation of the two incident half-edges

  for (int iE = 0; iE < nE; ++iE) {
    if (getNumberOfEdgeHalfEdges(iE) == 2) {
      const int iC1 = getEdgeHalfEdge(iE, 0);
      const int iC2 = getEdgeHalfEdge(iE, 1);
      partition.join(
        std::min(iC1, iC2), std::max(iC1, iC2));
      // ???
    }
  }

  // consistently oriented
  /* \                  / */
  /*  \ iC01 <-- iC00  /  */
  /*   X ---- iE ---- X   */
  /*  / iC10 --> iC11  \  */
  /* /                  \ */

  // opposite orientation
  /* \                  / */
  /*  \ iC01 --> iC00  /  */
  /*   X ---- iE ---- X   */
  /*  / iC10 --> iC11  \  */
  /* /                  \ */

  // a decision has to be made about inconsistently oriented faces
  // incident to the same edge, as well as how to deal with singular
  // edges; for the moment let's assume that the mesh does not have
  // singular edges, and that pairs of corners corresponding to the
  // same vertex across inconsistently oriented faces will be joined

  // note that the partition will end up with the corner separators as
  // singletons, but it doesn't matter for the last step, and
  // the partition will be deleted upon return
  
  // 4) count number of parts per vertex
  //    - initialize _nPartsVertex array to 0's
  _nPartsVertex.resize(nV, 0);

  //    - for each corner iC which is a representative of its subset, 
  //    - get the corresponding vertex index iV and increment _nPartsVertex[iV]
  //    - note that all the corners in each subset share a common
  //      vertex index, but multiple subsets may correspond to the
  //      same vertex index, indicating that the vertex is singular

  for (int iC = 0; iC < nC; ++iC) {
    // - for edge boundary iE label its two end vertices as boundary
    if (coordIndex[iC] < 0)
      continue;

    const int iV = getSrc(iC);
    _nPartsVertex[iV] += partition.getSize(iC);
  }
}

int PolygonMesh::getNumberOfFaces() const {
  return _numberOfFaces;
}

int PolygonMesh::getNumberOfEdgeFaces(const int iE) const {
  return getNumberOfEdgeHalfEdges(iE);
}

int PolygonMesh::getEdgeFace(const int iE, const int j) const {
  if (iE < 0 || iE >= getNumberOfEdges())
    return -1;

  if (j < 0 || j >= getNumberOfEdgeHalfEdges(iE))
    return -1;

  const int iC = getEdgeHalfEdge(iE, j);
  return getFace(iC);
}

bool PolygonMesh::isEdgeFace(const int iE, const int iF) const {
  if (iE < 0 || iE >= getNumberOfEdges())
    return false;

  for (int j = 0; j < getNumberOfEdgeHalfEdges(iE); ++j) {
    const int iC = getEdgeHalfEdge(iE, j);
    if (getFace(iC) == iF)
      return true;
  }
  return false;
}

// classification of edges

bool PolygonMesh::isBoundaryEdge(const int iE) const {

  return getNumberOfEdgeFaces(iE) == 1;
}

bool PolygonMesh::isRegularEdge(const int iE) const {
  return getNumberOfEdgeFaces(iE) == 2;
}

bool PolygonMesh::isSingularEdge(const int iE) const {
  return getNumberOfEdgeFaces(iE) > 1;
}

// classification of vertices

bool PolygonMesh::isBoundaryVertex(const int iV) const {
  const int nV = getNumberOfVertices();
  return (0<=iV && iV<nV)?_isBoundaryVertex[iV]:false;
}

bool PolygonMesh::isSingularVertex(const int iV) const {
  const int nV = getNumberOfVertices();
  return (0<=iV && iV<nV && _nPartsVertex[iV]>1);
}

// properties of the whole mesh

bool PolygonMesh::isRegular() const {

  for (int iE = 0; iE < getNumberOfEdges(); ++iE) {
    if (isSingularEdge(iE))
      return false;
  }

  for (int iV = 0; iV < getNumberOfVertices(); ++iV) {
    if (isSingularVertex(iV))
      return false;
  }

  return true;
}

bool PolygonMesh::hasBoundary() const {
  for (int iE = 0; iE < getNumberOfEdges(); ++iE) {
    if (isBoundaryEdge(iE))
      return true;
  }
  return false;
}
