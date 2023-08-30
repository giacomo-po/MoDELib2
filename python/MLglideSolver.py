import os
os.environ["DGLBACKEND"] = "pytorch"
import math
import numpy as np
import sys
import torch
import torch.nn as nn
import torchvision
import dgl
import dgl.data
import torch.nn.functional as F
from dgl.nn import GraphConv
import vtk

# dglContainer Class
# Linking a Machine Learning Model to the c++
class GCN(nn.Module):
    def __init__(self, in_feats, h_feats, num_classes):
        super(GCN, self).__init__()
        self.conv1 = GraphConv(in_feats, h_feats)
        self.conv2 = GraphConv(h_feats, num_classes)
        
    def forward(self, g, in_feat):
        h = self.conv1(g, in_feat)
        h = F.relu(h)
        h = self.conv2(g, h)
        return h
        
    def passing(self):
        src_ids = torch.tensor([1, 2, 3, 4])     # Source nodes for edges
        dst_ids = torch.tensor([0, 1, 2, 3])     # Destination nodes for edges
        g2 = dgl.graph((src_ids, dst_ids))
        l = torch.randn(5,5)
        g2 = dgl.add_self_loop(g2)
        return self.forward(g2,l)
    
def dglgraph():
    graph = GCN(5, 16, 2)  #Create the model with given dimensions
    return graph


# MultContainer Class
def mult(A,B):
    result = np.matmul(A,B)
    return result

