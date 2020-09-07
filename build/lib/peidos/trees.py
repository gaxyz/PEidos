#!/usr/bin/env python3
import dendropy as dp

# Tree functions
def compute_generation(node,firstSplitTime):
    if node.parent_node:
        node.generation = int(node.parent_node.generation) + int(node.edge.length)
    else:
        node.generation = firstSplitTime

def prepare_tree(tree, firstSplitTime):
    
    """
    Assign labels to internal nodes to properly write an eidos script.
    """
    tree.leaves = []
    for node in tree.preorder_node_iter():
        
        if node.label is None:
            new_label = node.taxon.__str__().replace("'","")
            node.label = new_label
        if node.is_leaf():
            node.slabel = node.label
            tree.leaves.append( node.slabel )
            
        compute_generation(node,firstSplitTime)
    
    for node in tree.postorder_node_iter():
        
        sister_slabels = []
        if node.parent_node: 
            for i in node.parent_node.child_nodes(): 
                sister_slabels.append( i.slabel )
            
            sister_slabels.sort()
            node.parent_node.slabel = sister_slabels[0]
            
        

def read_tree( treeFile, firstSplitTime ):
    tree = dp.Tree.get(path = treeFile, 
                   schema = "newick")
    prepare_tree(tree, firstSplitTime)
    
    return tree






