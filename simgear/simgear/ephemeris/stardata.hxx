// stardata.hxx -- manage star data
//
// Written by Curtis Olson, started March 2000.
//
// Copyright (C) 2000  Curtis L. Olson - http://www.flightgear.org/~curt
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
// $Id$


#ifndef _SG_STARDATA_HXX
#define _SG_STARDATA_HXX

#include <vector>
#include <simgear/math/SGMath.hxx>

class SGPath;

class SGStarData {
public:
    // Constructor
    SGStarData( const SGPath& path );

    // Destructor
    ~SGStarData();

    // load the stars database
    bool load( const SGPath& path );

    // stars
    inline int getNumStars() const { return static_cast<int>(_stars.size()); }
    inline SGVec3d *getStars() { return &(_stars[0]); }

private:
    std::vector<SGVec3d> _stars;
};


#endif // _SG_STARDATA_HXX
