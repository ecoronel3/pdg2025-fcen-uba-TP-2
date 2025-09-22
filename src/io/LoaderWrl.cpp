//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-05 16:36:18 taubin>
//------------------------------------------------------------------------
//
// LoaderWrl.cpp
//
// Software developed for the course
// Digital Geometry Processing
// Copyright (c) 2025, Gabriel Taubin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Brown University nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL GABRIEL TAUBIN BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#include "LoaderWrl.hpp"

#include <cstdio>
#include <stdexcept>

#include "TokenizerFile.hpp"

#define VRML_HEADER "#VRML V2.0 utf8"

const char* LoaderWrl::_ext = "wrl";

bool LoaderWrl::loadSceneGraph(TokenizerFile& tkn, SceneGraph& wrl) {

  string name    = "";
  bool   success = false;
  while(tkn.get()) {
    if(tkn.equals("DEF")) {
      tkn.get("missing token after DEF");
      name = tkn;
    } else if(tkn.equals("Group")) {
      Group* g = new Group();
      wrl.addChild(g);
      loadGroup(tkn,*g);
      g->setName(name);
      name = "";
    } else if(tkn.equals("Transform")) {
      Transform* t = new Transform();
      wrl.addChild(t);
      loadTransform(tkn,*t);
      t->setName(name);
      name = "";
    } else if(tkn.equals("Shape")) {
      Shape* s = new Shape();
      wrl.addChild(s);
      loadShape(tkn,*s);
      s->setName(name);
      name = "";
    } else if(tkn.equals("")) {
      break;
    } else {
      fprintf(stderr,"tkn=\"%s\"\n",tkn.c_str());
      throw std::runtime_error("unexpected token while parsing Group");
    }
  }
  success = true;
  return success;
}

bool LoaderWrl::loadGroup(TokenizerFile& tkn, Group& group) {

  // Group {
  //   MFNode children    []
  //   SFVec3f bboxCenter  0 0 0
  //   SFVec3f bboxSize   -1 -1 -1
  // }

  bool success = false;
  if(tkn.expecting("{")==false) throw std::runtime_error("expecting \"{\"");
  while(success==false && tkn.get()) {
    if(tkn.equals("children")) {
      loadChildren(tkn,group);
    } else if(tkn.equals("bboxCenter")) {
      Vec3f v;
      if(tkn.getVec3f(v)==false)
        throw std::runtime_error("expecting Vec3f");
      group.setBBoxCenter(v);
    } else if(tkn.equals("bboxSize")) {
      Vec3f v;
      if(tkn.getVec3f(v)==false)
        throw std::runtime_error("expecting Vec3f");
      group.setBBoxSize(v);
    } else if(tkn.equals("}")) {
      success = true;
    } else {
      throw std::runtime_error("unexpected token while parsing Group");
    }
  }
  return success;
}

bool LoaderWrl::loadTransform(TokenizerFile& tkn, Transform& transform) {

  // Transform {
  //   MFNode     children          []
  //   SFVec3f    bboxCenter        0 0 0
  //   SFVec3f    bboxSize          -1 -1 -1
  //   SFVec3f    center            0 0 0
  //   SFRotation rotation          0 0 1 0
  //   SFVec3f    scale             1 1 1
  //   SFRotation scaleOrientation  0 0 1 0
  //   SFVec3f    translation       0 0 0
  // }

  bool success = false;
  if(tkn.expecting("{")==false) throw std::runtime_error("expecting \"{\"");
  while(success==false && tkn.get()) {
    if(tkn.equals("children")) {
      loadChildren(tkn,transform);
    } else if(tkn.equals("bboxCenter")) {
      Vec3f v;
      if(tkn.getVec3f(v)==false)
        throw std::runtime_error("expecting Vec3f");
      transform.setBBoxCenter(v);
    } else if(tkn.equals("bboxSize")) {
      Vec3f v;
      if(tkn.getVec3f(v)==false)
        throw std::runtime_error("expecting Vec3f");
      transform.setBBoxCenter(v);
    } else if(tkn.equals("center")) {
      Vec3f v;
      if(tkn.getVec3f(v)==false)
        throw std::runtime_error("expecting Vec3f");
      transform.setCenter(v);
    } else if(tkn.equals("rotation")) {
      Vec4f	v;
      if(tkn.getVec4f(v)==false)
        throw std::runtime_error("expecting Vec4f");
      transform.setRotation(v);
    } else if(tkn.equals("scale")) {
      Vec3f v;
      if(tkn.getVec3f(v)==false)
        throw std::runtime_error("expecting Vec3f");
      transform.setScale(v);
    } else if(tkn.equals("scaleOrientation")) {
      // expecting 4 floats
      Vec4f v;
      if(tkn.getVec4f(v)==false)
        throw std::runtime_error("expecting Vec4f");
      transform.setScaleOrientation(v);
    } else if(tkn.equals("translation")) {
      // expecting 3 floats
      Vec3f v;
      if(tkn.getVec3f(v)==false)
        throw std::runtime_error("expecting Vec3f");
      transform.setTranslation(v);
    } else if(tkn.equals("}")) {
      success = true;
    } else {
      throw std::runtime_error("unexpected token while parsing Group");
    }
  }
  return success;
}

bool LoaderWrl::loadChildren(TokenizerFile& tkn, Group& group) {
  string name    = "";
  bool   success = false;
  if(tkn.expecting("[")==false) throw std::runtime_error("expecting \"[\"");
  while(success==false && tkn.get()) {
    if(tkn.equals("DEF")) {
      tkn.get("missing token after DEF");
      name = tkn;
    } else if(tkn.equals("Group")) {
      Group* g = new Group();
      group.addChild(g);
      loadGroup(tkn,*g);
      g->setName(name);
      name = "";
    } else if(tkn.equals("Transform")) {
      Transform* t = new Transform();
      group.addChild(t);
      loadTransform(tkn,*t); 
      t->setName(name);
      name = "";
   } else if(tkn.equals("Shape")) {
      Shape* s = new Shape();
      group.addChild(s);
      loadShape(tkn,*s);
      s->setName(name);
      name = "";
    } else if(tkn.equals("]")) {
      success = true;
    } else {
      throw std::runtime_error("unexpected token while parsing Group");
    }
  }
  return success;
}

bool LoaderWrl::loadShape(TokenizerFile& tkn, Shape& shape) {

  // Shape {
  //   SFNode appearance NULL
  //   SFNode geometry   NULL
  // }

  // TODO Mon Nov 12 10:22:28 2012
  // implement DEF name for appearance and geometry nodes

  string name    = "";
  bool   success = false;
  if(tkn.expecting("{")==false) throw std::runtime_error("expecting \"{\"");
  while(success==false && tkn.get()) {
    if(tkn.equals("appearance")) {
      tkn.get("expecting appearance node");
      if(tkn.equals("DEF")) {
        tkn.get("missing token after DEF");
        name = tkn;
        tkn.get("missing Appearance token");
      }
      if(tkn.equals("Appearance")==false)
        throw std::runtime_error("expecting Appearance");
      Appearance* a = new Appearance();
      a->setName(name);
      name = "";
      shape.setAppearance(a);
      loadAppearance(tkn,*a);
    } else if(tkn.equals("geometry")) {
      tkn.get("expecting geometry node");
      if(tkn.equals("DEF")) {
        tkn.get("missing token after DEF");
        name = tkn;
        tkn.get("missing Appearance token");
      }
      if(tkn.equals("IndexedFaceSet")) {
        IndexedFaceSet* ifs = new IndexedFaceSet();
        ifs->setName(name);
        name = "";
        shape.setGeometry(ifs);
        loadIndexedFaceSet(tkn,*ifs);
      } else if(tkn.equals("IndexedLineSet")) {
        IndexedLineSet* ils = new IndexedLineSet();
        ils->setName(name);
        name = "";
        shape.setGeometry(ils);
        loadIndexedLineSet(tkn,*ils);
      } else {
        throw std::runtime_error("found unexpected geometry node");
      }
    } else if(tkn.equals("}")) {
      success = true;
    } else {
      throw std::runtime_error("found Appearance field");
    }
  }
  return success;
}

bool LoaderWrl::loadAppearance(TokenizerFile& tkn, Appearance& appearance) {

  // Appearance {
  //   SFNode material NULL
  //   SFNode texture NULL
  //   // SFNode textureTransform NULL
  // }

  // TODO Mon Nov 12 10:22:28 2012
  // implement DEF name for matrial and texture nodes

  string name    = "";
  bool   success = false;
  if(tkn.expecting("{")==false) throw std::runtime_error("expecting \"[\"");
  while(success==false && tkn.get()) {
    if(tkn.equals("material")) {
      tkn.get("expecting material node");
      if(tkn.equals("DEF")) {
        tkn.get("missing token after DEF");
        name = tkn;
        tkn.get("missing Appearance token");
      }
      if(tkn.equals("Material")==false)
        throw std::runtime_error("expecting Material");
      Material* m = new Material();
      m->setName(name);
      name = "";
      appearance.setMaterial(m);
      loadMaterial(tkn,*m);
    } else if(tkn.equals("texture")) {
      tkn.get("expecting Texture node");
      if(tkn.equals("DEF")) {
        tkn.get("missing token after DEF");
        name = tkn;
        tkn.get("missing Appearance token");
      }
      if(tkn.equals("ImageTexture")) {
        ImageTexture* it = new ImageTexture();
        it->setName(name);
        name = "";
        appearance.setTexture(it);
        loadImageTexture(tkn,*it);
   // } else if(tkn.equals("PixelTexture")) {
   //   ...
   // } else if(tkn.equals("MovieTexture")) {
   //   ...
      } else {
        throw std::runtime_error("found unexpected Texture node");
      }
    } else if(tkn.equals("}")) {
      success = true;
    } else {
      throw std::runtime_error("found unexpected Group field");
    }
  }
  return success;
}

bool LoaderWrl::loadMaterial(TokenizerFile& tkn, Material& material) {

  // Material {
  //   SFFloat ambientIntensity 0.2
  //   SFColor diffuseColor     0.8 0.8 0.8
  //   SFColor emissiveColor    0 0 0
  //   SFFloat shininess        0.2
  //   SFColor specularColor    0 0 0
  //   SFFloat transparency     0
  // }

  bool success = false;
  if(tkn.expecting("{")==false) throw std::runtime_error("expecting \"{\"");
  while(success==false && tkn.get()) {
    if(tkn.equals("ambientIntensity")) {
      float f;
      if(tkn.getFloat(f)==false)
        throw std::runtime_error("expecting float");
      material.setAmbientIntensity(f);
    } else if(tkn.equals("diffuseColor")) {
      Color c;
      if(tkn.getColor(c)==false)
        throw std::runtime_error("expecting Color");
      material.setDiffuseColor(c);
    } else if(tkn.equals("emissiveColor")) {
      Color c;
      if(tkn.getColor(c)==false)
        throw std::runtime_error("expecting Color");
      material.setEmissiveColor(c);
    } else if(tkn.equals("shininess")) {
      float f;
      if(tkn.getFloat(f)==false)
        throw std::runtime_error("expecting float");
      material.setShininess(f);
    } else if(tkn.equals("specularColor")) {
      Color c;
      if(tkn.getColor(c)==false)
        throw std::runtime_error("expecting Color");
      material.setSpecularColor(c);
    } else if(tkn.equals("transparency")) {
      float f;
      if(tkn.getFloat(f)==false)
        throw std::runtime_error("expecting float");
      material.setTransparency(f);
    } else if(tkn.equals("}")) {
      success = true;
    } else {
      throw std::runtime_error("found unexpected Appearance field");
    }
  }
  return success;

}

bool LoaderWrl::loadImageTexture(TokenizerFile& tkn, ImageTexture& imageTexture) {

  // ImageTexture {
  //   MFString url []
  //   SFBool repeatS TRUE
  //   SFBool repeatT TRUE
  // }

  bool success = false;
  if(tkn.expecting("{")==false) throw std::runtime_error("expecting \"{\"");
  while(success==false && tkn.get()) {
    if(tkn.equals("url")) {
      vector<string>& _url = imageTexture.getUrl();
      if(loadVecString(tkn,_url)==false)
        throw std::runtime_error("loading vector<string>");
    } else if(tkn.equals("repeatS")) {
      bool b;
      if(tkn.getBool(b)==false)
        throw std::runtime_error("expecting boolean value");
      imageTexture.setRepeatS(b);
    } else if(tkn.equals("repeatT")) {
      bool b;
      if(tkn.getBool(b)==false)
        throw std::runtime_error("expecting boolean value");
      imageTexture.setRepeatT(b);
    } else if(tkn.equals("}")) {
      success = true;
    } else {
      throw std::runtime_error("found unexpected Appearance field");
    }
  }
  return success;
}

bool LoaderWrl::loadIndexedFaceSet(TokenizerFile& tkn, IndexedFaceSet& ifs) {

  // IndexedFaceSet {
  //   SFNode  color             NULL
  //   SFNode  coord             NULL
  //   SFNode  normal            NULL
  //   SFNode  texCoord          NULL
  //   SFBool  ccw               TRUE
  //   MFInt32 colorIndex        []        # [-1,)
  //   SFBool  colorPerVertex    TRUE
  //   SFBool  convex            TRUE
  //   MFInt32 coordIndex        []        # [-1,)
  //   SFFloat creaseAngle       0         # [ 0,)
  //   MFInt32 normalIndex       []        # [-1,)
  //   SFBool  normalPerVertex   TRUE
  //   SFBool  solid             TRUE
  //   MFInt32 texCoordIndex     []        # [-1,)
  // }

  bool success = false;
  if(tkn.expecting("{")==false) throw std::runtime_error("expecting \"{\"");
  while(success==false && tkn.get()) {
    if(tkn.equals("color")) {
      //   SFNode  
      vector<float>& _color = ifs.getColor();

      // if DEF name found, skip and ignore ??

      if(tkn.expecting("Color")==false)
        throw std::runtime_error("expecting Color");
      if(tkn.expecting("{")==false)
        throw std::runtime_error("expecting \"{\"");
      if(tkn.expecting("color")==false)
        throw std::runtime_error("expecting color");
      if(loadVecFloat(tkn,_color)==false)
        throw std::runtime_error("loading IndexedFaceSet color field");
      if(tkn.expecting("}")==false)
        throw std::runtime_error("expecting \"}\"");

    } else if(tkn.equals("coord")) {
      //   SFNode  
      vector<float>& _coord = ifs.getCoord();

      // if DEF name found, skip and ignore ??

      if(tkn.expecting("Coordinate")==false)
        throw std::runtime_error("expecting Coordinate");
      if(tkn.expecting("{")==false)
        throw std::runtime_error("expecting \"{\"");
      if(tkn.expecting("point")==false)
        throw std::runtime_error("expecting point");
      if(loadVecFloat(tkn,_coord)==false)
        throw std::runtime_error("loading IndexedFaceSet coord field");
      if(tkn.expecting("}")==false)
        throw std::runtime_error("expecting \"}\"");

    } else if(tkn.equals("normal")) {
      //   SFNode  
      vector<float>& _normal = ifs.getNormal();

      // if DEF name found, skip and ignore ??

      if(tkn.expecting("Normal")==false)
        throw std::runtime_error("expecting Normal");
      if(tkn.expecting("{")==false)
        throw std::runtime_error("expecting \"{\"");
      if(tkn.expecting("vector")==false)
        throw std::runtime_error("expecting vector");
      if(loadVecFloat(tkn,_normal)==false)
        throw std::runtime_error("loading IndexedFaceSet normal field");
      if(tkn.expecting("}")==false)
        throw std::runtime_error("expecting \"}\"");

    } else if(tkn.equals("texCoord")) {
      //   SFNode  
      vector<float>& _texCoord = ifs.getTexCoord();

      // if DEF name found, skip and ignore ??

      if(tkn.expecting("TextureCoordinate")==false)
        throw std::runtime_error("expecting TextureCoordinate");
      if(tkn.expecting("{")==false)
        throw std::runtime_error("expecting \"{\"");
      if(tkn.expecting("point")==false)
        throw std::runtime_error("expecting point");
      if(loadVecFloat(tkn,_texCoord)==false)
        throw std::runtime_error("loading IndexedFaceSet texCoord field");
      if(tkn.expecting("}")==false)
        throw std::runtime_error("expecting \"}\"");

    } else if(tkn.equals("ccw")) {
      //   SFBool
      bool& _ccw = ifs.getCcw();
      if(tkn.getBool(_ccw)==false)
        throw std::runtime_error("loading IndexedFaceSet ccw field");  
    } else if(tkn.equals("colorIndex")) {
      //   MFInt32 
      vector<int>& _colorIndex = ifs.getColorIndex();
      if(loadVecInt(tkn,_colorIndex)==false)
        throw std::runtime_error("loading IndexedFaceSet colorIndex field");
    } else if(tkn.equals("colorPerVertex")) {
      //   SFBool
      bool& _colorPerVertex = ifs.getColorPerVertex();
        if(tkn.getBool(_colorPerVertex)==false)
        throw std::runtime_error("loading IndexedFaceSet colorPerVertex field");
    } else if(tkn.equals("convex")) {
      //   SFBool
      bool& _convex = ifs.getConvex();  
      if(tkn.getBool(_convex)==false)
        throw std::runtime_error("loading IndexedFaceSet convex field");
    } else if(tkn.equals("coordIndex")) {
      //   MFInt32 
      vector<int>& _coordIndex = ifs.getCoordIndex();
      if(loadVecInt(tkn,_coordIndex)==false)
        throw std::runtime_error("loading IndexedFaceSet coordIndex field");
    } else if(tkn.equals("creaseAngle")) {
      //   SFFloat
      float& _creaseAngle = ifs.getCreaseangle();
      if(tkn.getFloat(_creaseAngle)==false)
        throw std::runtime_error("loading IndexedFaceSet creaseAngle value");
    } else if(tkn.equals("normalIndex")) {
      //   MFInt32 
      vector<int>& _normalIndex = ifs.getNormalIndex();
      if(loadVecInt(tkn,_normalIndex)==false)
        throw std::runtime_error("loading IndexedFaceSet normalIndex field");
    } else if(tkn.equals("normalPerVertex")) {
      //   SFBool
      bool& _normalPerVertex = ifs.getNormalPerVertex();
      if(tkn.getBool(_normalPerVertex)==false)
        throw std::runtime_error("loading IndexedFaceSet normalPerVertex field");
    } else if(tkn.equals("solid")) {
      //   SFBool
      bool& _solid = ifs.getSolid();  
      if(tkn.getBool(_solid)==false)
        throw std::runtime_error("loading IndexedFaceSet solid field");
    } else if(tkn.equals("texCoordIndex")) {
      //   MFInt32 
      vector<int>& _texCoordIndex = ifs.getTexCoordIndex();
      if(loadVecInt(tkn,_texCoordIndex)==false)
        throw std::runtime_error("loading IndexedFaceSet texCoordIndex field");
    } else if(tkn.equals("}")) {
      success = true;
    } else {
      throw std::runtime_error("found unexpected IndexedFaceSet field");
    }
  }
  return success;
}

bool LoaderWrl::loadIndexedLineSet(TokenizerFile& tkn, IndexedLineSet& ifs) {

  // IndexedFaceSet {
  //   SFNode  coord             NULL
  //   MFInt32 coordIndex        []        # [-1,)
  //   SFNode  color             NULL
  //   MFInt32 colorIndex        []        # [-1,)
  //   SFBool  colorPerVertex    TRUE
  // }

  bool success = false;
  if(tkn.expecting("{")==false) throw std::runtime_error("expecting \"{\"");
  while(success==false && tkn.get()) {
    if(tkn.equals("color")) {
      //   SFNode  
      vector<float>& _color = ifs.getColor();

      // if DEF name found, skip and ignore ??

      if(tkn.expecting("Color")==false)
        throw std::runtime_error("expecting Color");
      if(tkn.expecting("{")==false)
        throw std::runtime_error("expecting \"{\"");
      if(tkn.expecting("color")==false)
        throw std::runtime_error("expecting color");
      if(loadVecFloat(tkn,_color)==false)
        throw std::runtime_error("loading IndexedLineSet color field");
      if(tkn.expecting("}")==false)
        throw std::runtime_error("expecting \"}\"");

    } else if(tkn.equals("coord")) {
      //   SFNode  
      vector<float>& _coord = ifs.getCoord();

      // if DEF name found, skip and ignore ??

      if(tkn.expecting("Coordinate")==false)
        throw std::runtime_error("expecting Coordinate");
      if(tkn.expecting("{")==false)
        throw std::runtime_error("expecting \"{\"");
      if(tkn.expecting("point")==false)
        throw std::runtime_error("expecting point");
      if(loadVecFloat(tkn,_coord)==false)
        throw std::runtime_error("loading IndexedLineSet coord field");
      if(tkn.expecting("}")==false)
        throw std::runtime_error("expecting \"}\"");

    } else if(tkn.equals("colorIndex")) {
      //   MFInt32 
      vector<int>& _colorIndex = ifs.getColorIndex();
      if(loadVecInt(tkn,_colorIndex)==false)
        throw std::runtime_error("loading IndexedLineSet colorIndex field");
    } else if(tkn.equals("colorPerVertex")) {
      //   SFBool
      bool& _colorPerVertex = ifs.getColorPerVertex();
        if(tkn.getBool(_colorPerVertex)==false)
        throw std::runtime_error("loading IndexedLineSet colorPerVertex field");
    } else if(tkn.equals("coordIndex")) {
      //   MFInt32 
      vector<int>& _coordIndex = ifs.getCoordIndex();
      if(loadVecInt(tkn,_coordIndex)==false)
        throw std::runtime_error("loading IndexedLineSet coordIndex field");
    } else if(tkn.equals("}")) {
      success = true;
    } else {
      throw std::runtime_error("found unexpected IndexedLineSet field");
    }
  }
  return success;
}

bool LoaderWrl::loadVecFloat(TokenizerFile&tkn,vector<float>& vec) {
  bool success = false;
  if(tkn.expecting("[")==false) throw std::runtime_error("expecting \"[\"");
  float value;
  while(success==false && tkn.get()) {
    if(tkn.equals("]")) {
      success = true; // done
    } else if(sscanf(tkn.c_str(),"%f",&value)==1) {
      vec.push_back(value);
    } else {
      throw std::runtime_error("expecting int value");
    }
  }

  return success;
}

bool LoaderWrl::loadVecInt(TokenizerFile&tkn,vector<int>& vec) {
  bool success = false;
  if(tkn.expecting("[")==false) throw std::runtime_error("expecting \"[\"");
  int value;
  while(success==false && tkn.get()) {
    if(tkn.equals("]")) {
      success = true; // done
    } else if(sscanf(tkn.c_str(),"%d",&value)==1) {
      vec.push_back(value);
    } else {
      throw std::runtime_error("expecting int value");
    }
  }
  return success;
}

bool LoaderWrl::loadVecString(TokenizerFile&tkn,vector<string>& vec) {
  bool success = false;
  tkn.get("expecting a token");
  if(tkn.equals("[")) {
    // expecting 0 or more strings followed by "]"
    while(tkn.get()) {
      if(tkn.equals("]"))
        break;
      else
        vec.push_back(tkn);
    }
    success = true;
  } else {
    // expecting a single string
    tkn.get("expecting a token");
    vec.push_back(tkn);
    success = true;
  }
  return success;
}

bool LoaderWrl::load(const char* filename, SceneGraph& wrl) {
  bool success = false;

  FILE* fp = (FILE*)0;
  try {

    // open the file
    if(filename==(char*)0) throw std::runtime_error("filename==null");
    fp = fopen(filename,"r");
    if(fp==(FILE*)0) throw std::runtime_error("fp==(FILE*)0");

    // clear the container
    wrl.clear();
    wrl.setUrl(filename);

    // read and check header line
    char header[16];
    // memset(header,'\0',16);
    for(int i=0;i<16;i++) header[i] = '\0';
    fscanf(fp,"%15c",header);
    if(string(header)!=VRML_HEADER) throw std::runtime_error("header!=VRM_HEADER");

    // create a TokenizerFile and start parsing
    TokenizerFile tkn(fp);
    loadSceneGraph(tkn,wrl);

    // will be done later
    // wrl.updateBBox();
    
    // if we have reached this point we have succeeded
    fclose(fp);
    success = true;

  } catch(const std::exception& e) {

    if(fp!=(FILE*)0) fclose(fp);
    fprintf(stderr,"ERROR | %s\n",e.what());
    wrl.clear();
    wrl.setUrl("");

  }

  return success;
}

