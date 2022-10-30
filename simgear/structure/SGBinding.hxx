/**
 * \file commands.hxx
 * Interface definition for encapsulated commands.
 * Started Spring 2001 by David Megginson, david@megginson.com
 * This code is released into the Public Domain.
 *
 * $Id$
 */

#ifndef __SGBINDING_HXX
#define __SGBINDING_HXX


#include <simgear/compiler.h>

#include <string>
#include <map>
#include <vector>

#include <simgear/props/props.hxx>
#include <simgear/props/condition.hxx>

#include "commands.hxx"

/**
 * An input binding of some sort.
 *
 * <p>This class represents a binding that can be assigned to a
 * keyboard key, a joystick button or axis, or even a panel
 * instrument.</p>
 */
class SGBinding : public SGConditional
{
public:

  /**
   * Default constructor.
   */
  SGBinding ();

    /**
     * Convenience constructor.
     *
     * @param commandName TODO
     */
    SGBinding(const std::string& commandName);

  /**
   * Convenience constructor.
   *
   * @param node The binding will be built from this node.
   * @param root Property root used while building binding.
   */
  SGBinding( const SGPropertyNode *node,
             SGPropertyNode *root );


  /**
   * Destructor.
   */
  virtual ~SGBinding () = default;

  
  /**
   * Clear internal state of the binding back to empty.
   *
   * This was particularly useful when SGBinding's destructor had its 'remove
   * on delete' behaviour, however this is not the case anymore.
   */
  void clear();


  /**
   * Get the command name.
   *
   * @return The string name of the command for this binding.
   */
  const std::string &getCommandName () const { return _command_name; }

  /**
   * Get the argument that will be passed to the command.
   *
   * @return A property node that will be passed to the command as its
   * argument, or 0 if none was supplied.
   */
  const SGPropertyNode * getArg () { return _arg; }
  

  /**
   * Read a binding from a property node.
   *
   * @param node The property node containing the binding.
   * @param root The property root node used while building the binding from
   *             \a node.
   */
  void read (const SGPropertyNode * node, SGPropertyNode* root);


  /**
   * Fire a binding.
   */
  void fire () const;


  /**
   * Fire a binding with a scaled movement (rather than absolute position).
   */
  void fire (double offset, double max) const;


  /**
   * Fire a binding with a setting (i.e. joystick axis).
   *
   * A double 'setting' property will be added to the arguments.
   *
   * @param setting The input setting, usually between -1.0 and 1.0.
   */
  void fire (double setting) const;

  /**
   * Fire a binding with a number of additional parameters
   * 
   * The children of params will be merged with the fixed arguments.
   */
  void fire (SGPropertyNode* params) const;
  
private:
  void innerFire() const;
                                // just to be safe.
  SGBinding (const SGBinding &binding);

  std::string _command_name;
  mutable SGPropertyNode_ptr _arg;
  mutable SGPropertyNode_ptr _setting;
  mutable SGPropertyNode_ptr _root;
};

typedef SGSharedPtr<SGBinding> SGBinding_ptr;

typedef std::vector<SGBinding_ptr > SGBindingList;
typedef std::map<unsigned,SGBindingList> SGBindingMap;

/**
 * fire every binding in a list, in sequence
 
 */
void fireBindingList(const SGBindingList& aBindings, SGPropertyNode* params = NULL);

/**
 * fire every binding in a list with a setting value
 
 */
void fireBindingListWithOffset(const SGBindingList& aBindings, double offset, double max);

/**
 * read multiple bindings from property-list format
 */
SGBindingList readBindingList(const simgear::PropertyList& aNodes, SGPropertyNode* aRoot);

/**
 * call clear() on every binding in a list
 */
void clearBindingList(const SGBindingList& aBindings);

/**
 * check if at least one binding in the list is enabled. Returns false if bindings
 * list is empty, or all bindings are conditinally disabled.
 */
bool anyBindingEnabled(const SGBindingList& bindings);

#endif
