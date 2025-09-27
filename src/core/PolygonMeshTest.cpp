//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-05 16:34:30 taubin>
//------------------------------------------------------------------------
//
// PolygonMeshTest.cpp
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

#include "PolygonMeshTest.hpp"

#include <cassert>
#include <iostream>
#include <wrl/Shape.hpp>
#include <wrl/Appearance.hpp>
#include <wrl/Material.hpp>
#include <wrl/IndexedFaceSet.hpp>
#include <wrl/SceneGraphTraversal.hpp>

PolygonMeshTest::PolygonMeshTest(SceneGraph& sceneGraph, const std::string& indent, std::ostream& ostr):_ostr(ostr) {
  _ostr << indent << "PolygonMeshTest {" << endl;

  int nIndexedFaceSet = 0;

  SceneGraphTraversal traversal(sceneGraph);
  Node* node = nullptr;
  while((node=traversal.next())!=nullptr) {
    if(node->isShape()) {
      _ostr << indent << "  Shape {" << endl;
      auto shape = dynamic_cast<Shape*>(node);

      _ostr << indent << "    name = \"" << shape->getName() << "\"" << endl;    

      node = shape->getAppearance();
      if(node==nullptr) {
        _ostr << indent << "    appearance = NULL" << endl;    
      } else if(node->isAppearance()) {
        auto appearance = dynamic_cast<Appearance*>(node);
        node = appearance->getMaterial();
        if(node==nullptr) {
          _ostr << indent << "    appearance->material = NULL" << endl;    
        } else if(node->isMaterial()) {
          auto* material = dynamic_cast<Material*>(node);
          const Color& dc = material->getDiffuseColor();
          _ostr << indent
                << "    diffuseColor = [ "
                << dc.r << " " << dc.g << " " << dc.b
                << " ]" << endl;    
        }
      }

      node = shape->getGeometry();
      if(node==nullptr) {
        _ostr << indent << "    geometry = NULL" << endl;    
      } else if(node->isIndexedFaceSet()) {
        _ostr << indent
              << "    geometry IndexedFaceSet[" << nIndexedFaceSet << "] {" << endl;
        auto* ifs = dynamic_cast<IndexedFaceSet*>(node);

        int nVifs = ifs->getNumberOfCoord();
        vector<int>& coordIndex = ifs->getCoordIndex();

        _ostr << indent << "      nV(ifs) = " << nVifs << endl;

        _ostr << indent << "      PolygonMesh(nV,coordIndex) {" << endl;

        PolygonMesh pMesh(nVifs,coordIndex);
        assert(pMesh.getDst(2) == 2);
        assert(pMesh.getDst(6) == 3);
        assert(pMesh.getDst(5) == 0);

        assert(pMesh.getNext(0) == 1);
        assert(pMesh.getNext(1) == 2);
        assert(pMesh.getNext(2) == 0);

        assert(pMesh.getNext(8) == 9);
        assert(pMesh.getNext(9) == 10);
        assert(pMesh.getNext(10) == 8);

        assert(pMesh.getNext(12) == 13);
        assert(pMesh.getNext(13) == 14);
        assert(pMesh.getNext(14) == 12);

        assert(pMesh.getPrev(10) == 9);
        assert(pMesh.getPrev(9) == 8);
        assert(pMesh.getPrev(8) == 10);

        assert(pMesh.getPrev(2) == 1);
        assert(pMesh.getPrev(1) == 0);
        assert(pMesh.getPrev(0) == 2);

        assert(pMesh.getPrev(14) == 13);
        assert(pMesh.getPrev(13) == 12);
        assert(pMesh.getPrev(12) == 14);

        assert(pMesh.getTwin(0) == 14);
        assert(pMesh.getTwin(8) == 13);
        assert(pMesh.getTwin(6) == 9);

        assert(pMesh.getNumberOfEdgeHalfEdges(0) == 2);
        assert(pMesh.getNumberOfEdgeHalfEdges(1) == 2);
        assert(pMesh.getNumberOfEdgeHalfEdges(2) == 2);
        assert(pMesh.getNumberOfEdgeHalfEdges(3) == 2);
        assert(pMesh.getNumberOfEdgeHalfEdges(4) == 2);
        assert(pMesh.getNumberOfEdgeHalfEdges(5) == 2);

        assert(pMesh.getEdgeHalfEdge(0, 0) == 0);
        assert(pMesh.getEdgeHalfEdge(0, 1) == 14);
        assert(pMesh.getEdgeHalfEdge(0, 2) == -1);

        int nV = pMesh.getNumberOfVertices();
        int nE = pMesh.getNumberOfEdges();
        int nF = pMesh.getNumberOfFaces();
        int nC = pMesh.getNumberOfCorners();

        _ostr << indent << "        nV          = " << nV << endl;
        _ostr << indent << "        nE          = " << nE << endl;
        _ostr << indent << "        nF          = " << nF << endl;
        _ostr << indent << "        nC          = " << nC << endl;

        // print info about the polygon mesh

        int nV_boundary  = 0;
        int nV_internal  = 0;
        int nV_singular  = 0;
        int nV_regular   = 0;
        int nE_boundary  = 0;
        int nE_regular   = 0;
        int nE_singular  = 0;
        int nE_other     = 0;

        int iE,iV;

        for(iE=0;iE<nE;iE++) {
          if(pMesh.isBoundaryEdge(iE)) {
            nE_boundary++;
          } else if(pMesh.isRegularEdge(iE)) {
            nE_regular++;
          } else if(pMesh.isSingularEdge(iE)) {
            nE_singular++;
          } else {
            nE_other++;
          }
        }

        for(iV=0;iV<nV;iV++) {
          if(pMesh.isBoundaryVertex(iV))
            nV_boundary++;
          if(pMesh.isSingularVertex(iV))
            nV_singular++;
        }

        nV_internal = nV-nV_boundary;
        nV_regular  = nV-nV_singular;

        _ostr << indent << "        nV_boundary = " << nV_boundary << endl;
        _ostr << indent << "        nV_internal = " << nV_internal << endl;
        _ostr << indent << "        nV_regular  = " << nV_regular  << endl;
        _ostr << indent << "        nV_singular = " << nV_singular << endl;
        _ostr << indent << "        nE_boundary = " << nE_boundary << endl;
        _ostr << indent << "        nE_regular  = " << nE_regular  << endl;
        _ostr << indent << "        nE_singular = " << nE_singular << endl;
        _ostr << indent << "        nE_other    = " << nE_other    << endl;
        _ostr << indent << "        isRegular   = " << pMesh.isRegular() << endl;
        _ostr << indent << "        hasBoundary = " << pMesh.hasBoundary() << endl;

        _ostr << indent << "      } PolygonMesh" << endl;
        _ostr << indent << "    } IndexedFaceSet" << endl;
        nIndexedFaceSet++;
      } else {
        _ostr << indent << "    geometry " << node->getType() << endl;    
      }

      _ostr << indent << "  } Shape" << endl;
    }
  }

  _ostr << indent << "} PolygonMeshTest" << endl;
}
