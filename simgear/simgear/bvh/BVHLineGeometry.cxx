// Copyright (C) 2008 - 2009  Mathias Froehlich - Mathias.Froehlich@web.de
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//

#include "BVHLineGeometry.hxx"

#include "BVHVisitor.hxx"

namespace simgear {

BVHLineGeometry::BVHLineGeometry(const SGLineSegmentf& lineSegment, Type type) :
    _lineSegment(lineSegment),
    _type(type)
{
}

BVHLineGeometry::~BVHLineGeometry()
{
}

void
BVHLineGeometry::accept(BVHVisitor& visitor)
{
    visitor.apply(*this);
}

SGSphered
BVHLineGeometry::computeBoundingSphere() const
{
    SGSphered sphere;
    sphere.expandBy(SGVec3d(_lineSegment.getStart()));
    sphere.expandBy(SGVec3d(_lineSegment.getEnd()));
    return sphere;
}

}
