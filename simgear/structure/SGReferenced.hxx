/* -*-c++-*-
 *
 * Copyright (C) 2005-2006 Mathias Froehlich
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 */

#ifndef SGReferenced_HXX
#define SGReferenced_HXX

#include "SGAtomic.hxx"

/// Base class for all reference counted SimGear objects
/// Classes derived from this one are meant to be managed with
/// the SGSharedPtr class.
///
/// For more info see SGSharedPtr. For using weak references see
/// SGWeakReferenced.

class SGReferenced {
public:
  SGReferenced(void) : _refcount(0u)
  {}
  /// Do not copy reference counts. Each new object has it's own counter
  SGReferenced(const SGReferenced&) : _refcount(0u)
  {}
  /// Do not copy reference counts. Each object has it's own counter
  SGReferenced& operator=(const SGReferenced&)
  { return *this; }

  static unsigned get(const SGReferenced* ref)
  { if (ref) return ++(ref->_refcount); else return 0; }
  static unsigned put(const SGReferenced* ref) noexcept
  { if (ref) return --(ref->_refcount); else return 0; }
  static unsigned count(const SGReferenced* ref)
  { if (ref) return ref->_refcount; else return 0; }
  static bool shared(const SGReferenced* ref)
  { if (ref) return 1u < ref->_refcount; else return false; }

private:
  mutable SGAtomic _refcount;
};

#endif
