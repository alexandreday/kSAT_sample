import numpy as np

class TREENODE:

    def __init__(self, id_ = -1, parent = None, child = None, value = 0):
        if child is None:
            self.child = [] # has to be list of TreeNode
        else:
            self.child = child
        self.value = value
        self.parent = parent
        self.id_ = id_

    def __repr__(self):
        return ("Node: [%s] @ s = %.3f" % (self.id_,self.scale)) 

    def is_leaf(self):
        return len(self.child) == 0

    def get_child(self, id_ = None):
        if id_ is None:
            return self.child
        else:
            for c in self.child:
                if c.get_id() == id_:
                    return c

    def get_child_id(self):
        if len(self.child) == 0:
            return []
        else:
            return [c.id_ for c in self.child]

    def get_id(self):
        return self.id_

    def add_child(self, treenode):
        self.child.append(treenode)
    
    def remove_child(self, treenode):
        self.child.remove(treenode)
        
    def get_rev_child(self):
        child = self.child[:]
        child.reverse()
        return child 
    
class DPLL:
    def __init__(self):
        return 
