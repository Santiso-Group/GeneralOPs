/*
** Copyright 2012 Erik Santiso.
** This file is part of mymath.
** mymath is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
**
** mymath is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with mymath. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef H_TREE
#define H_TREE

#include <iostream>
#include <algorithm>
#include <utility>
#include <vector>
#include <deque>
#include <string>
#include "common/include/types.h"
#include "common/include/assert.h"

/**
 * Tree class
 *
 * This implements a very simple templated indexed rooted tree.
 *
 * TODO: Figure out how to group the index-returning methods' documentation
 *       in Doxygen.
 *
 * Notes:
 *
 * - Would be much more elegant to define iterators for different ways to
 *   traverse the tree rather than have all the functions returning
 *   std::vector<size_t>'s. Since this is intended for rather small trees
 *   atm will leave the iterators for a later implementation.
 * 
 * - Also, would be much more efficient to handle this using references, but
 *   again, this is intended to handle small trees.
 *
 * - Might be good to have mutable members for the depth(s).
 */

static size_t const notANode = (size_t)-1;

template<typename Type>
class Tree
{
public:

  // Constructors

  /// Direct constructor
  /**
  * Defines a tree with a (optionally) given root node
  */
  Tree(Type const &root = Type());

  // Destructor

  /// Default destructor (virtual to allow for derivation)
  virtual ~Tree() = default;

// Interface

  /// Add a child to the i-th node, returns the index of the child
  size_t const addChild(size_t const i, Type const &child);
  /// Add a subtree as a child of the i-th node
  void const addSubtree(size_t const i, Tree<Type> const &subtree);
  /// Clear the tree - resulting tree has only the (optionally) given root node
  void clear(Type const &root = Type());
  /// Replace the content of the i-th node
  void setNode(size_t const i, Type const &node);
  /// Removes the i-th node and all its descendents from the tree
  void const prune(size_t const i);

// Accessors

  /// Get the depth of the i-th node
  size_t const getDepth(size_t const i) const;
  /// Get the depth of the tree
  size_t const getDepth() const;
  /// Get the content of the i-th node
  Type const &getNode(size_t const i) const;
  /// Get the number of children of the i-th node
  size_t const getNumChildren(size_t const i) const;
  /// Get the number of nodes of the tree
  size_t const getNumNodes() const;
  /// Get the parent of the i-th node
  size_t const &getParent(size_t const i) const;
  /// Get the subtree that has node i as its root node
  Tree<Type> const getSubtree(size_t const i) const;
  /// Test whether node is a leaf node
  bool const isLeaf(size_t const i) const;

// The following return std::vectors<size_t> containing the node indices
// that satisfy a particular property. Note that, if nodes are deleted, it
// is unsafe to traverse the tree by using raw index values - use one of
// these lists instead.

  /// Get a std::vector containing the ancestors of the given node
  std::vector<size_t> const getAncestors(size_t const i) const;
  /// Get a std::vector containing the children of the given node
  std::vector<size_t> const &getChildren(size_t const i) const;
  /// Get a std::vector containing the descendants of the given node
  std::vector<size_t> const getDescendants(size_t const i) const;
  /// Get a std::vector containing the leaf nodes
  std::vector<size_t> const getLeaves() const;
  /// Get a std::vector containing the nodes at the given depth
  std::vector<size_t> const getLevel(size_t const depth) const;
  /// Get a std::vector containing all the nodes
  std::vector<size_t> const getNodes() const;
  /// Get a std::vector containing all the nodes, in breadth-first order
  std::vector<size_t> const getNodesBreadthFirst() const;
  /// Get a std::vector containing all the nodes, in depth-first order
  std::vector<size_t> const getNodesDepthFirst() const;
  /// Get a std::vector containing the siblings of the given node
  std::vector<size_t> const getSiblings(size_t const i) const;

// Output

  /// Write tree topology in ASCII (mostly for debugging/visualization)
  template<typename NodeType>
  friend std::ostream &operator<<(std::ostream &outStream, 
                                  Tree<NodeType> const &tree);

protected:

// Members

  std::vector<Type> nodes_;                    // The actual nodes
  std::vector<size_t> parents_;                // Parent node of each node
  std::vector<std::vector<size_t> > children_; // Children of each node
  std::vector<size_t> free_;                   // Free space (handles erasing)

// Methods

  bool const isFree(size_t const i) const ; // Checks whether the i-th node has
                                            // been erased.
};

/*
** End of class Tree
*/

// Inlines

template<typename Type>
inline Tree<Type>::Tree(Type const &root)
:
nodes_(1, root), parents_(1, notANode), 
children_(1, std::vector<size_t>()), free_()
{}

template<typename Type>
inline size_t const Tree<Type>::addChild(size_t const i, Type const &child)
{
  assert(i < nodes_.size() && !isFree(i));
  if(free_.size() != 0)
  {
    size_t const indx = free_[free_.size() - 1];
    free_.pop_back();
    nodes_[indx] = child;
    parents_[indx] = i;
    children_[indx] = std::vector<size_t>();
    children_[i].push_back(indx);
    return indx;
  }
  else
  {
    size_t const indx = nodes_.size();
    nodes_.push_back(child);
    parents_.push_back(i);
    children_.push_back(std::vector<size_t>());
    children_[i].push_back(indx);
    return indx;
  }
}

template<typename Type>
void const Tree<Type>::addSubtree(size_t const i, Tree<Type> const &subtree)
{
  assert(i < nodes_.size() && !isFree(i));
  std::deque<size_t> nodesToVisit(1, 0);
  std::deque<size_t> parentNodes(1, i);

  while(nodesToVisit.size() > 0)
  {
    size_t const currentNode = nodesToVisit[0];
    nodesToVisit.pop_front();
    size_t const currentParent = parentNodes[0];
    parentNodes.pop_front();
    size_t newNode = addChild(currentParent, subtree.getNode(currentNode));
    std::vector<size_t> const currentChildren = 
      subtree.getChildren(currentNode);
    for(std::vector<size_t>::const_iterator it = currentChildren.begin();
        it != currentChildren.end(); ++it)
    {
      nodesToVisit.push_back(*it);
      parentNodes.push_back(newNode);
    }
  }
}

template<typename Type>
inline void Tree<Type>::clear(Type const &root)
{
  nodes_ = std::vector<Type>(1, root);
  parents_ = std::vector<size_t>(1, notANode);
  children_ = std::vector<std::vector<size_t> >(1, std::vector<size_t>());
  free_.clear();
}

template<typename Type>
inline void Tree<Type>::setNode(size_t const i, Type const &node)
{
  assert(i < nodes_.size() && !isFree(i));
  nodes_[i] = node;
}

template<typename Type>
void const Tree<Type>::prune(size_t const i)
{
  assert(i < nodes_.size() && !isFree(i));
  assert(i > 0); // Can't prune root node
  std::vector<size_t> const descendants = getDescendants(i);
  // Clear node and update parent's children list
  size_t const parent = parents_[i];
  std::vector<size_t>::iterator it = 
    std::find(children_[parent].begin(), children_[parent].end(), i);
  assert(it != children_[parent].end()); // This would be a bug
  children_[parent].erase(it);
  parents_[i] = notANode;
  children_[i].clear();
  free_.push_back(i);
  // Clear all the descendants
  for(std::vector<size_t>::const_iterator cit = descendants.begin();
      cit != descendants.end(); ++cit)
  {
    parents_[*cit] = notANode;
    children_[*cit].clear();
    free_.push_back(*cit);
  }
}

// None of the methods below modify the tree

template<typename Type>
inline size_t const Tree<Type>::getDepth(size_t const i) const
{
  assert(i < nodes_.size() && !isFree(i));
  size_t depth = 0;
  size_t current = i;
  while(current != 0)
  {
    current = parents_[current];
    ++depth;
  }
  return depth; 
}

template<typename Type>
inline size_t const Tree<Type>::getDepth() const
{
  size_t depth = 0;
  std::vector<size_t> const leaves = getLeaves();
  for(size_t i = 0; i < leaves.size(); ++i)
  {
    size_t const leafDepth = getDepth(leaves[i]);
    if(leafDepth > depth) depth = leafDepth;
  }
  return depth;
}

template<typename Type>
inline Type const &Tree<Type>::getNode(size_t const i) const
{
  assert(i < nodes_.size() && !isFree(i));
  return nodes_[i];
}

template<typename Type>
inline size_t const Tree<Type>::getNumChildren(size_t const i) const
{
  assert(i < nodes_.size() && !isFree(i));
  return children_[i].size();
}

template<typename Type>
inline size_t const Tree<Type>::getNumNodes() const
{
  return nodes_.size() - free_.size();
}

template<typename Type>
inline size_t const &Tree<Type>::getParent(size_t const i) const
{
  if(i == 0) return notANode; // Root has no parent
  assert(i < nodes_.size() && !isFree(i));
  return parents_[i];
}

template<typename Type>
inline Tree<Type> const Tree<Type>::getSubtree(size_t const i) const
{
  assert(i < nodes_.size() && !isFree(i));
  Tree<Type> subtree(getNode(i));
  std::deque<size_t> nodesToVisit(1, i);
  std::deque<size_t> nodesToVisitSub(1, 0);
  while(nodesToVisit.size() > 0)
  {
    size_t const currentNode = nodesToVisit[0];
    nodesToVisit.pop_front();
    size_t const currentNodeSub = nodesToVisitSub[0];
    nodesToVisitSub.pop_front();
    std::vector<size_t> const currentChildren = getChildren(currentNode);
    for(std::vector<size_t>::const_iterator it = currentChildren.begin();
        it != currentChildren.end(); ++it)
    {
      nodesToVisit.push_back(*it);
      size_t indxSub = subtree.addChild(currentNodeSub, getNode(*it));
      nodesToVisitSub.push_back(indxSub);
    }
  }
  return subtree;
}

template<typename Type>
inline bool const Tree<Type>::isLeaf(size_t const i) const
{
  assert(i < nodes_.size() && !isFree(i));
  return (children_[i].size() == 0);
}

template<typename Type>
inline std::vector<size_t> const Tree<Type>::getAncestors(size_t const i) const
{
  if(i == 0) return std::vector<size_t>(); // Root has no ancestors
  assert(i < nodes_.size() && !isFree(i));
  std::vector<size_t> ancestors;
  size_t currentNode = i;
  while(currentNode != 0)
  {
    currentNode = parents_[currentNode];
    ancestors.push_back(currentNode);
  }
  return ancestors;
}

template<typename Type>
inline std::vector<size_t> const &Tree<Type>::getChildren(size_t const i) const
{
  assert(i < nodes_.size() && !isFree(i));
  return children_[i];
}

template<typename Type>
inline std::vector<size_t> const 
  Tree<Type>::getDescendants(size_t const i) const
{
  assert(i < nodes_.size() && !isFree(i));
  std::vector<size_t> descendants;
  std::deque<size_t> nodesToVisit(1, i);
  while(nodesToVisit.size() > 0)
  {
    size_t const currentNode = nodesToVisit[0];
    std::vector<size_t> const currentChildren = getChildren(currentNode);
    if(currentNode != i) descendants.push_back(currentNode);
    nodesToVisit.pop_front();
    for(std::vector<size_t>::const_reverse_iterator rit 
          = currentChildren.rbegin();
        rit != currentChildren.rend(); ++rit)
      nodesToVisit.push_front(*rit);
  }
  return descendants;
}

template<typename Type>
inline std::vector<size_t> const Tree<Type>::getLeaves() const
{
  std::vector<size_t> leaves;
  for(size_t i = 0; i < nodes_.size(); ++i)
  {
    if(isFree(i)) continue;
    if(children_[i].size() == 0) leaves.push_back(i);
  }
  return leaves;
}

template<typename Type>
inline std::vector<size_t> const Tree<Type>::getLevel(size_t const depth) const
{
  assert(depth <= getDepth());
  std::vector<size_t> level;
  for(size_t i = 0; i < nodes_.size(); ++i)
  {
    if(isFree(i)) continue;
    if(getDepth(i) == depth) level.push_back(i);
  }
  return level;
}

template<typename Type>
inline std::vector<size_t> const Tree<Type>::getNodes() const
{
  std::vector<size_t> nodes;
  for(size_t i = 0; i < nodes_.size(); ++i)
    if(!isFree(i)) nodes.push_back(i);
  return nodes;
}

template<typename Type>
inline std::vector<size_t> const Tree<Type>::getNodesBreadthFirst() const
{
  std::vector<size_t> nodes;
  std::deque<size_t> nodesToVisit(1, 0);
  while(nodesToVisit.size() > 0)
  {
    size_t const currentNode = nodesToVisit[0];
    nodes.push_back(currentNode);
    nodesToVisit.pop_front();
    std::vector<size_t> const currentChildren = getChildren(currentNode);
    for(std::vector<size_t>::const_iterator it = currentChildren.begin();
        it != currentChildren.end(); ++it)
      nodesToVisit.push_back(*it);
  }
  return nodes;
}

template<typename Type>
inline std::vector<size_t> const Tree<Type>::getNodesDepthFirst() const
{
  std::vector<size_t> nodes;
  std::deque<size_t> nodesToVisit(1, 0);
  while(nodesToVisit.size() > 0)
  {
    size_t const currentNode = nodesToVisit[0];
    nodes.push_back(currentNode);
    nodesToVisit.pop_front();
    std::vector<size_t> const currentChildren = getChildren(currentNode);
    for(std::vector<size_t>::const_reverse_iterator rit 
          = currentChildren.rbegin();
        rit != currentChildren.rend(); ++rit)
      nodesToVisit.push_front(*rit);
  }
  return nodes;
}

template<typename Type>
inline std::vector<size_t> const Tree<Type>::getSiblings(size_t const i) const
{
  if(i == 0) return std::vector<size_t>(); // Root has no siblings
  assert(i < nodes_.size() && !isFree(i));
  std::vector<size_t> siblings = getChildren(parents_[i]);
  std::vector<size_t>::iterator it = 
    std::find(siblings.begin(), siblings.end(), i);
  siblings.erase(it);
  return siblings;
}

template<typename Type>
inline std::ostream &operator<<(std::ostream &outStream, Tree<Type> const &tree)
{
  // The idea here is similar to the getNodesDepthFirst() method, only we want 
  // to keep track of formatting and padding for special cases (last siblings,
  // nodes with/without children)

  outStream << std::endl;

  // Stacks to keep track of tree topology and formatting
  std::deque<size_t> nodesToVisit(1, 0);
  std::deque<std::string> indents(1, "");
  std::deque<bool> lastSiblings(1, true);

  while(nodesToVisit.size() > 0)
  {
    size_t const currentNode = nodesToVisit[0];
    nodesToVisit.pop_front();
    std::string indent = indents[0];
    indents.pop_front();
    bool const isLast = lastSiblings[0];
    lastSiblings.pop_front();

    std::vector<size_t> const currentChildren = tree.getChildren(currentNode);

// Uncomment to display tree indices - works regardless of node type
    outStream << indent << "+--" << currentNode << std::endl;

// Uncomment to display node contents - requires operator<< for node object
//    outStream << indent << "+--" << tree.getNode(currentNode) << std::endl;

    // Change indents, add formatting for nodes with children and last siblings
    indent += !isLast?"|  ":"   ";
    if(currentChildren.size() > 0) outStream << indent << "|" << std::endl;
    if(currentChildren.size() == 0 && isLast) outStream << indent << std::endl;

    // Update stacks
    for(std::vector<size_t>::const_reverse_iterator rit 
          = currentChildren.rbegin();
        rit != currentChildren.rend(); ++rit)
    {
      nodesToVisit.push_front(*rit);
      lastSiblings.push_front(rit == currentChildren.rbegin());
      indents.push_front(indent);
    }
  } // End of loop over nodes

  return outStream;
}
/*
template<typename Type>
inline std::pair<Tree<Type>, Tree<Type> > const
  Tree<Type>::crossover(Tree<Type> const &firstTree, size_t const firstNode, 
                        Tree<Type> const &secondTree, size_t const secondNode)
{
  Tree<Type> first = firstTree;
  Tree<Type> second = secondTree;
  Tree<Type> firstSub = first.getSubtree(firstNode);
  Tree<Type> secondSub = second.getSubtree(secondNode);
  size_t const firstParent = first.getParent(firstNode);
  size_t const secondParent = second.getParent(secondNode);
  first.prune(firstNode);
  second.prune(secondNode);
  first.addSubtree(firstParent, secondSub);
  second.addSubtree(secondParent, firstSub);
  return std::pair<Tree<Type>, Tree<Type> >(first, second);
}
*/
template<typename Type>
inline bool const Tree<Type>::isFree(size_t const i) const
{
  return std::find(free_.begin(), free_.end(), i) != free_.end();
}

#endif

