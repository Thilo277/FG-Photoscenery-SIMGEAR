///@file
/// Window for placing a Canvas onto it (for dialogs, menus, etc.)
//
// Copyright (C) 2012  Thomas Geymayer <tomgey@gmail.com>
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
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA

#include <simgear_config.h>
#include "CanvasMgr.hxx"
#include "CanvasSystemAdapter.hxx"
#include "CanvasWindow.hxx"

#include <simgear/canvas/Canvas.hxx>
#include <simgear/misc/strutils.hxx>
#include <simgear/scene/util/OsgMath.hxx>

#include <osgGA/GUIEventHandler>

namespace simgear
{
namespace canvas
{

  //----------------------------------------------------------------------------
  const std::string Window::TYPE_NAME = "window";

  //----------------------------------------------------------------------------
  Window::Window( const CanvasWeakPtr& canvas,
                  const SGPropertyNode_ptr& node,
                  const Style& parent_style,
                  Element* parent ):
    Image(canvas, node, parent_style, parent),
    _resize_top(node, "resize-top"),
    _resize_right(node, "resize-right"),
    _resize_bottom(node, "resize-bottom"),
    _resize_left(node, "resize-left"),
    _resize_status(node, "resize-status")
  {
    node->setFloatValue("source/right", 1);
    node->setFloatValue("source/bottom", 1);
    node->setBoolValue("source/normalized", true);
  }

  //----------------------------------------------------------------------------
  Window::~Window()
  {
    if( _canvas_decoration )
      _canvas_decoration->destroy();
  }

  //----------------------------------------------------------------------------
  void Window::update(double delta_time_sec)
  {
    if( _attributes_dirty & DECORATION )
    {
      updateDecoration();
      _attributes_dirty &= ~DECORATION;
    }

    Image::update(delta_time_sec);
  }

  //----------------------------------------------------------------------------
  void Window::valueChanged(SGPropertyNode * node)
  {
    bool handled = false;
    if( node->getParent() == _node )
    {
      handled = true;
      const std::string& name = node->getNameString();
      if( name  == "resize" )
        _resizable = node->getBoolValue();
      else if( name == "update" )
        update(0);
      else if( name == "capture-events" )
        _capture_events = node->getBoolValue();
      else if( name == "decoration-border" )
        parseDecorationBorder(node->getStringValue());
      else if(    strutils::starts_with(name, "shadow-")
               || name == "content-size" )
        _attributes_dirty |= DECORATION;
      else
        handled = false;
    }

    if( !handled )
      Image::valueChanged(node);
  }

  //----------------------------------------------------------------------------
  const SGVec2<float> Window::getPosition() const
  {
    auto const& m = getMatrix();
    return SGVec2<float>( m(3, 0), m(3, 1) );
  }

  //----------------------------------------------------------------------------
  const SGRect<float> Window::getScreenRegion() const
  {
    return getPosition() + getRegion();
  }

  //----------------------------------------------------------------------------
  void Window::setCanvasContent(CanvasPtr canvas)
  {
    _canvas_content = canvas;
    if( _layout )
    {
      canvas->setLayout(_layout);
      _layout.clear();
    }

    if( _image_content )
      // Placement within decoration canvas
      _image_content->setSrcCanvas(canvas);
    else
      setSrcCanvas(canvas);
  }

  //----------------------------------------------------------------------------
  CanvasWeakPtr Window::getCanvasContent() const
  {
    return _canvas_content;
  }

  //----------------------------------------------------------------------------
  void Window::setLayout(const LayoutRef& layout)
  {
    CanvasPtr canvas = _canvas_content.lock();
    if( canvas )
      canvas->setLayout(layout);
    else
      _layout = layout; // keep layout until content canvas is set
  }

  //----------------------------------------------------------------------------
  CanvasPtr Window::getCanvasDecoration() const
  {
    return _canvas_decoration;
  }

  //----------------------------------------------------------------------------
  bool Window::isResizable() const
  {
    return _resizable;
  }

  //----------------------------------------------------------------------------
  bool Window::isCapturingEvents() const
  {
    return _capture_events;
  }

  //----------------------------------------------------------------------------
  void Window::setVisible(bool visible)
  {
    LayoutItem::setVisible(visible);
    Element::setVisible(LayoutItem::isVisible());
  }

  //----------------------------------------------------------------------------
  bool Window::isVisible() const
  {
    return Element::isVisible();
  }

  //----------------------------------------------------------------------------
  void Window::raise()
  {
    // on writing the z-index the window always is moved to the top of all other
    // windows with the same z-index.
    set<int>("z-index", get<int>("z-index", 0));
  }

  //----------------------------------------------------------------------------
  void Window::handleResize( uint8_t mode,
                             const osg::Vec2f& offset )
  {
    if( mode == NONE )
    {
      _resize_status = 0;
      return;
    }
    else if( mode & INIT )
    {
      _resize_top    = getRegion().t();
      _resize_right  = getRegion().r();
      _resize_bottom = getRegion().b();
      _resize_left   = getRegion().l();
      _resize_status = 1;
    }

    if( mode & BOTTOM )
      _resize_bottom = getRegion().b() + offset.y();
    else if( mode & TOP )
      _resize_top = getRegion().t() + offset.y();

    if( mode & canvas::Window::RIGHT )
      _resize_right = getRegion().r() + offset.x();
    else if( mode & canvas::Window::LEFT )
      _resize_left = getRegion().l() + offset.x();
  }

  //----------------------------------------------------------------------------
  bool Window::handleEvent(const EventPtr& event)
  {
    if( auto mouse_event = dynamic_cast<MouseEvent*>(event.get()) )
    {
      mouse_event->local_pos =
      mouse_event->client_pos =
        mouse_event->screen_pos - toOsg(getPosition());
    }

    return Image::handleEvent(event);
  }

  //----------------------------------------------------------------------------
  void Window::parseDecorationBorder(const std::string& str)
  {
    _decoration_border = simgear::CSSBorder::parse(str);
    _attributes_dirty |= DECORATION;
  }

  //----------------------------------------------------------------------------
  void Window::updateDecoration()
  {
    int shadow_radius = get<float>("shadow-radius") + 0.5;
    if( shadow_radius < 2 )
      shadow_radius = 0;

    CanvasPtr content = _canvas_content.lock();
    SGRect<int> content_view
    (
      0,
      0,
      get<int>("content-size[0]", content ? content->getViewWidth()  : 400),
      get<int>("content-size[1]", content ? content->getViewHeight() : 300)
    );

    if( _decoration_border.isNone() && !shadow_radius )
    {
      setSrcCanvas(content);
      set<int>("size[0]", content_view.width());
      set<int>("size[1]", content_view.height());

      _image_content.reset();
      _image_shadow.reset();
      if( _canvas_decoration )
        _canvas_decoration->destroy();
      _canvas_decoration.reset();
      return;
    }

    if( !_canvas_decoration )
    {
      CanvasMgr* mgr = dynamic_cast<CanvasMgr*>(
        Canvas::getSystemAdapter()->getSubsystem("Canvas")
      );

      if( !mgr )
      {
        SG_LOG(SG_GENERAL, SG_WARN, "canvas::Window: no canvas manager!");
        return;
      }

      _canvas_decoration = mgr->createCanvas("window-decoration");
      _canvas_decoration->getProps()
                        ->setStringValue("background", "rgba(0,0,0,0)");
      setSrcCanvas(_canvas_decoration);

      _image_content = _canvas_decoration->getRootGroup()
                                         ->createChild<Image>("content");
      _image_content->setSrcCanvas(content);

      // Forward keyboard events to content
      _image_content->setFocus();

      // Draw content on top of decoration
      _image_content->set<int>("z-index", 1);
    }

    GroupPtr group_decoration =
      _canvas_decoration->getOrCreateGroup("decoration");
    group_decoration->set<int>("tf/t[0]", shadow_radius);
    group_decoration->set<int>("tf/t[1]", shadow_radius);
    // TODO do we need clipping or shall we trust the decorator not to draw over
    //      the shadow?

    CSSBorder::Offsets const border =
      _decoration_border.getAbsOffsets(content_view);

    int shad2 = 2 * shadow_radius,
        outer_width  = border.l + content_view.width()  + border.r + shad2,
        outer_height = border.t + content_view.height() + border.b + shad2;

    _canvas_decoration->setSizeX( outer_width );
    _canvas_decoration->setSizeY( outer_height );
    _canvas_decoration->setViewWidth( outer_width );
    _canvas_decoration->setViewHeight( outer_height );

    set<int>("size[0]", outer_width - shad2);
    set<int>("size[1]", outer_height - shad2);
    set<int>("outset", shadow_radius);

    assert(_image_content);
    _image_content->set<int>("x", shadow_radius + border.l);
    _image_content->set<int>("y", shadow_radius + border.t);
    _image_content->set<int>("size[0]", content_view.width());
    _image_content->set<int>("size[1]", content_view.height());

    if( !shadow_radius )
    {
      if( _image_shadow )
      {
        _image_shadow->destroy();
        _image_shadow.reset();
      }
      return;
    }

    int shadow_inset = std::max<int>(get<float>("shadow-inset") + 0.5, 0),
        slice_width  = shadow_radius + shadow_inset;

    _image_shadow = _canvas_decoration->getRootGroup()
                                      ->getOrCreateChild<Image>("shadow");
    _image_shadow->set<std::string>("src", "gui/images/shadow.png");
    _image_shadow->set<float>("slice", 7);
    _image_shadow->set<std::string>("fill", "#000000");
    _image_shadow->set<float>("slice-width", slice_width);
    _image_shadow->set<int>("size[0]", outer_width);
    _image_shadow->set<int>("size[1]", outer_height);

    // Draw shadow below decoration
    _image_shadow->set<int>("z-index", -1);
  }

} // namespace canvas
} // namespace simgear
