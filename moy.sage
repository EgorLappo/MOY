class X():
    
    def __init__(self, elist):
        # read in a crossing as a list of 4 edges
        assert elist[0][1] == elist[2][1] and elist[1][1] == elist[3][1]
        self.crossing = list(elist)
    
    @property
    def i(self):
        return self.crossing[0]
    
    @property
    def j(self):
        return self.crossing[1]
    
    @property
    def k(self):
        return self.crossing[2]
    
    @property
    def l(self):
        return self.crossing[3]
    
    @property
    def type(self):
        return 'x'
    
    @property
    def seg_ids(self):
        return [s[0] for s in self.crossing]
    
    
    @property
    def sign(self):
        if self.crossing[1][0] - self.crossing[3][0] == 1 or self.crossing[3][0] - self.crossing[1][0] > 1:
            return 1
        elif self.crossing[3][0] - self.crossing[1][0] == 1 or self.crossing[1][0] - self.crossing[3][0] > 1:
            return -1
        
    def __repr__(self):
        return 'Crossing '+str([i[0] for i in self.crossing])

class V():
    
    def __init__(self, inc, out):
        self.incoming_edges = inc
        self.outgoing_edges = out
        self.inc_len = len(inc)
        self.out_len = len(out)
        
    @property
    def incoming(self):
        return self.incoming_edges
    
    @property
    def outgoing(self):
        return self.outgoing_edges
    
    @property
    def n_incoming(self):
        return len(self.incoming_edges)
    
    @property
    def n_outgoing(self):
        return len(self.outgoing_edges)
    
    @property
    def in_seg_ids(self):
        return [s[0] for s in self.incoming_edges]
    
    @property
    def out_seg_ids(self):
        return [s[0] for s in self.outgoing_edges]
    
    @property
    def type(self):
        return 'v'
    
    def __repr__(self):
        return 'Vertex with ' + str([i[0] for i in self.incoming_edges]) \
             + ' incoming and '+ str([i[0] for i in self.outgoing_edges]) + ' outgoing.'


class MOY():
    def __init__(self, xvs, ):
        self.t = var('t')
        self.abs_delta = 0
        
        # store the initial list of crossings and vertices
        self.xvs = xvs  
        
        #separate it into crossings and vertices
        self.xs = [i for i in xvs if i.type == 'x']
        self.vs = [i for i in xvs if i.type == 'v']
        
        # we get segments as inputs, so we want to go from segments to actual edges and vice versa
        max_seg_id = -1
        for x in self.xs:
            max_seg_id = max(max_seg_id, max(x.seg_ids))
        for v in self.vs:
            max_seg_id = max(max_seg_id, max(v.in_seg_ids), max(v.out_seg_ids))
        self.n_seg = max_seg_id + 1
        
        # this is a standard union-find implementation to find the maps from segments to edges and back
        self.segment_parents = [i for i in range(self.n_seg)]
        
        for x in self.xs:
            self._union(self.segment_parents, x.i[0], x.k[0])
            self._union(self.segment_parents, x.j[0], x.l[0])
        
        # this list gives an (arbitrary) index of an edge for each segment
        self.segment_parents = [self._find(self.segment_parents, i) for i in range(self.n_seg)]
        
        # self.edges stores the indices of segments that represent edges, as chosed by union-find
        self.n_edges = len(set(self.segment_parents))
        self.edges = list(set(self.segment_parents))
        
        self.edge_to_seg = {i:[j for j in range(self.n_seg) if self._find(self.segment_parents, j) == i] for i in self.edges}
        self.seg_to_edge = self.segment_parents
        
        
        # this is a straightforward way to record a color of each segment
        self.seg_to_col = {}
        for x in self.xs:
            for s in x.crossing:
                if s[0] not in self.seg_to_col:
                    self.seg_to_col[s[0]] = s[1]
        
        for v in self.vs:
            for s in v.incoming:
                if s[0] not in self.seg_to_col:
                    self.seg_to_col[s[0]] = s[1]
            for s in v.outgoing:
                if s[0] not in self.seg_to_col:
                    self.seg_to_col[s[0]] = s[1]
        
        # in case we ever need this, we also record a color for each edge in a separate list
        self.edge_to_col = {e:self.seg_to_col[self.edge_to_seg[e][0]] for e in self.edges}
                    
        # |Cr(D)| in the paper
        self.n_crossings = len(self.xs) + sum([v.n_incoming for v in self.vs])
        
        # under the assumption of a connected graph diagram, we get |Re(D)| as follows
        self.n_regions = self.n_crossings + 2
        
        # now, we need to somehow store which segments are adjacent to which regions
        # so, we first have 2*n_segments 'candidate regions' -- on the left and right of each segment
        # then, for each crossing we unite all we can:
        # for example, in positive X(i,j,k,l) we unite 
        #    - 'right of i' and 'right of j'
        #    - 'left of j'  and 'right of k'
        #    - 'left of k'  and 'left of l'
        #    - 'right of l' and 'left of i'
        # in negative X(i,j,k,l) we unite 
        #    - 'right of i' and 'left of j'
        #    - 'right of j' and 'right of k'
        #    - 'left of k'  and 'right of l'
        #    - 'left of l'  and 'left of i'
        
        # for each vertex, we first do not want to introduce the circle regions
        # we take care about them later, since we would know which other two regions are adjacent to them!
        # so, we just unite left of first incoming with left of first outgoing, etc.
        
        reg_c = [i for i in range(2*self.n_seg)]
        for x in self.xs:
            if x.sign == 1:
                # 2*x.i -- left; 2*x.i + 1 -- right
                self._union(reg_c, 2*x.i[0] + 1, 2*x.j[0] + 1)
                self._union(reg_c, 2*x.j[0],     2*x.k[0] + 1)
                self._union(reg_c, 2*x.k[0],     2*x.l[0])
                self._union(reg_c, 2*x.l[0] + 1, 2*x.i[0])
            elif x.sign == -1:
                self._union(reg_c, 2*x.i[0] + 1, 2*x.j[0])
                self._union(reg_c, 2*x.j[0] + 1, 2*x.k[0] + 1)
                self._union(reg_c, 2*x.k[0],     2*x.l[0] + 1)
                self._union(reg_c, 2*x.l[0],     2*x.i[0])
            else: 
                print('something is really wrong!')
                
        for v in self.vs:
            self._union(reg_c, 2*v.incoming[0][0], 2*v.outgoing[-1][0])
            self._union(reg_c, 2*v.incoming[-1][0]+1, 2*v.outgoing[0][0]+1)
            if v.n_incoming > 1:
                for i in range(v.n_incoming - 1):
                    self._union(reg_c, 2*v.incoming[i][0]+1, 2*v.incoming[i+1][0])
            if v.n_outgoing > 1:
                for i in range(v.n_outgoing - 1):
                    self._union(reg_c, 2*v.outgoing[i][0], 2*v.outgoing[i+1][0]+ 1)
                    
        self.reg_parents = [self._find(reg_c, i) for i in range(2*self.n_seg)]
        
        # seg to region in the form (r_left, r_right), i. e. for each segment we record 
        # a region on the left and on the right
        self.seg_to_reg = {i:[self.reg_parents[2*i], self.reg_parents[2*i+1]] for i in range(self.n_seg)}
        
        # this list stores representatives of noncircular regions as with edges above
        self.n_noncirc_regions = len(set(self.reg_parents))
        self.noncirc_regions = list(set(self.reg_parents))
        
        # for each region we store the segments that bound it
        self.reg_to_seg = {i:list(set([j//2 for j, r in enumerate(self.reg_parents) if r == i])) for i in self.noncirc_regions}
        
        # for each region we store the crossings and vertices adjacent to it
        # we encode it by putting a tuple: an object itself plus some additional info.
        # for crossings:  position of a region
        # for vertices: whether the region is adjacent to incoming edges in case we need it

        # also add if's to avoid duplicates (for example, planar handcuffs would give duplicates without it)

        self.reg_to_xvs = {i:[] for i in self.noncirc_regions}
        for x in self.xs:
            if x.sign == 1:
                if (x, 'S') not in self.reg_to_xvs[self.seg_to_reg[x.i[0]][0]]:
                    self.reg_to_xvs[self.seg_to_reg[x.i[0]][0]].append((x,'S'))
                if (x, 'E') not in self.reg_to_xvs[self.seg_to_reg[x.i[0]][1]]:
                    self.reg_to_xvs[self.seg_to_reg[x.i[0]][1]].append((x,'E'))
                if (x, 'W') not in self.reg_to_xvs[self.seg_to_reg[x.k[0]][0]]:
                    self.reg_to_xvs[self.seg_to_reg[x.k[0]][0]].append((x,'W'))
                if (x, 'N') not in self.reg_to_xvs[self.seg_to_reg[x.k[0]][1]]:
                    self.reg_to_xvs[self.seg_to_reg[x.k[0]][1]].append((x,'N'))
            elif x.sign == -1:
                if (x, 'W') not in self.reg_to_xvs[self.seg_to_reg[x.i[0]][0]]:
                    self.reg_to_xvs[self.seg_to_reg[x.i[0]][0]].append((x,'W'))
                if (x, 'S') not in self.reg_to_xvs[self.seg_to_reg[x.i[0]][1]]:
                    self.reg_to_xvs[self.seg_to_reg[x.i[0]][1]].append((x,'S'))
                if (x, 'N') not in self.reg_to_xvs[self.seg_to_reg[x.k[0]][0]]:
                    self.reg_to_xvs[self.seg_to_reg[x.k[0]][0]].append((x,'N'))
                if (x, 'E') not in self.reg_to_xvs[self.seg_to_reg[x.k[0]][1]]:
                    self.reg_to_xvs[self.seg_to_reg[x.k[0]][1]].append((x,'E'))
        for  v in self.vs:
            for s, _ in v.incoming:
                if (v, 'inc') not in self.reg_to_xvs[self.seg_to_reg[s][0]]:
                    self.reg_to_xvs[self.seg_to_reg[s][0]].append((v,'inc'))
            if (v, 'inc') not in self.reg_to_xvs[self.seg_to_reg[v.incoming[-1][0]][1]]:
                self.reg_to_xvs[self.seg_to_reg[v.incoming[-1][0]][1]].append((v,'inc'))
            
            if len(v.outgoing) > 1:
                for i in range(len(v.outgoing)-1):
                    s, _ = v.outgoing[i]
                    if (v, 'inc') not in self.reg_to_xvs[self.seg_to_reg[s][0]]:
                        self.reg_to_xvs[self.seg_to_reg[s][0]].append((v,'not inc'))
        
        self.A = self.get_alexander_matrix()
        
    # just a check
    def is_MOY(self):
        return all([sum([c for s, c in v.outgoing]) == sum([c for s, c in v.outgoing]) for v in self.vs])
    
    def get_alexander_matrix(self):
        # get the dimensions of a matrix
        nrows = self.n_crossings
        ncols = self.n_regions
        
        # use sage to define a symbolic matrix of zeros
        A = matrix(SR, nrows, ncols)
        
        # first we iterate over the crossings, filling out all of the related entries of the matrix
        # (it is more convenient to make the outer loop go over regions)
        for j, (reg, adj_xvs) in enumerate(self.reg_to_xvs.items()):
            # this is a helper array to check if we have an adjacency
            adj_xvs_elems_only = [e[0] for e in adj_xvs]
            for i, x in enumerate(self.xs):
                if x in adj_xvs_elems_only:
                    contrib = 0
                    # this covers the possibility of the same crossing being present multiple times 
                    # for the same region 
                    x_present = [e for e in adj_xvs if e[0] == x]
                    for x, pos in x_present:
                        x_type = '+' if x.sign == 1 else '-'
                        # extract the color of a proper edge
                        col = self.seg_to_col[x.j[0]]
                        contrib += self.m_Cp(x_type, pos)*self.A_Cp(x_type, pos, col)
                    A[i,j] += contrib
        # ^ this way we have filled out the first len(xs) rows, corresponding to the crossings 
        # that do not arise from vertices. Next, i switch and iterate over the remaining crossings, adding contributions to each row
    
        # is still the row index that i will increment
        i = len(self.xs)
        # k serves to keep track of the circle region
        for k, v in enumerate(self.vs):
            # for each incoming edge we get a crossing
            for s, col in v.incoming:
                left_region = list(self.reg_to_xvs.keys()).index(self.seg_to_reg[s][0])
                circ_region = k + len(self.reg_to_xvs.keys())
                right_region = list(self.reg_to_xvs.keys()).index(self.seg_to_reg[s][1])

                A[i,left_region] += self.m_Cp('circ', 'left')*self.A_Cp('circ', 'left', col)
                A[i,circ_region] += self.m_Cp('circ', 'top')*self.A_Cp('circ', 'top', col)
                A[i,right_region] += self.m_Cp('circ', 'right')*self.A_Cp('circ', 'right', col)
                i += 1
        return A
    
    def get_state_sum(self, delta_segment = 0, left_region_index = 0):
        # first create a delta factor, as described in Definition 2.10
        n = left_region_index
        self.abs_delta = self.t**(n - self.seg_to_col[delta_segment]) - self.t**(n)  
        
        # now delete the two columns and compute the deretminant
        R_u, R_v = self.seg_to_reg[delta_segment]
        square_A = self.A[:,[i for i in range(self.n_regions) if i != R_u and i != R_v]]
        
        return square_A.det()/self.abs_delta
    
    def get_alexander_invariant(self, delta_segment = 0, left_region_index = 0):
        ss = self.get_state_sum(delta_segment, left_region_index)
        denum = (self.t**(-1/2)- self.t**(1/2))**(len(self.vs) - 1)
        return ss/denum

    def get_square_A(self, delta_segment = 0, left_region_index = 0):
        n = left_region_index
        self.abs_delta = self.t**(n - self.seg_to_col[delta_segment]) - self.t**(n)  
        
        # now delete the two columns and compute the deretminant
        R_u, R_v = self.seg_to_reg[delta_segment]
        square_A = self.A[:,[i for i in range(self.n_regions) if i != R_u and i != R_v]]
        return square_A
                
    # get the m_Cp contribution of a vertex to a region
    # given the type of vertex and the relative location of the region
    def m_Cp(self, xv_type, pos):   
        if xv_type == 'circ':
            if pos == 'top':
                return 1
            elif pos == 'left':
                return -1
            elif pos == 'right':
                return 1
        elif xv_type == '+':
            if pos == 'N':
                return 1
            elif pos == 'S':
                return 1
            elif pos == 'W':
                return -1
            elif pos == 'E':
                return -1
        elif xv_type == '-':
            if pos == 'N':
                return -1
            elif pos == 'S':
                return -1
            elif pos == 'W':
                return 1
            elif pos == 'E':
                return 1
            
    # same for the other contribution A_Cp
    def A_Cp(self, xv_type, pos, color):
        if xv_type == 'circ':
            if pos == 'top':
                return self.t**(-color/2) - self.t**(color/2)
            elif pos == 'left':
                return self.t**(-color/2)
            elif pos == 'right':
                return self.t**(color/2)
        elif xv_type == '+':
            if pos == 'N':
                return self.t**(-color)
            elif pos == 'S':
                return 1
            elif pos == 'W':
                return self.t**(-color)
            elif pos == 'E':
                return 1
        elif xv_type == '-':
            if pos == 'N':
                return self.t**(color)
            elif pos == 'S':
                return 1
            elif pos == 'W':
                return 1
            elif pos == 'E':
                return self.t**(color)
        
    def _find(self, l, seg):
        if seg == l[seg]:
            return seg
        l[seg] = self._find(l, l[seg])
        return l[seg]
        
    def _union(self, l, seg_i, seg_j):
        iloc = self._find(l, seg_i)
        jloc = self._find(l, seg_j)
        
        if iloc <= jloc:
            l[jloc] = iloc
        else:
            l[iloc] = jloc
        return 0



# some debug code

# def find(l, seg):
#     if seg == l[seg]:
#         return seg
#     l[seg] = find(l, l[seg])
#     return l[seg]

# def union(l, seg_i, seg_j):
#     iloc = find(l, seg_i)
#     jloc = find(l, seg_j)

#     if iloc <= jloc:
#         l[jloc] = iloc
#     else:
#         l[iloc] = jloc
#     return 0

# reg_c = [i for i in range(2*hc1.n_seg)] 
# print(reg_c)
# for x in hc1.xs:
#     if x.sign == 1:
#         # 2*x.i -- left; 2*x.i + 1 -- right
#         union(reg_c, 2*x.i[0] + 1, 2*x.j[0] + 1)
#         union(reg_c, 2*x.j[0],     2*x.k[0] + 1)
#         union(reg_c, 2*x.k[0],     2*x.l[0])
#         union(reg_c, 2*x.l[0] + 1, 2*x.i[0])
        
#         print("Joining left of " + str(x.j[0])+ " and right of " + str(x.k[0]))
#         print("Joining left of " + str(x.i[0])+ " and right of " + str(x.l[0]))
#         print("Joining left of " + str(x.k[0])+ " and left of " + str(x.l[0]))
#         print("Joining right of " + str(x.i[0])+ " and right of " + str(x.j[0]))
#     elif x.sign == -1:
#         union(reg_c, 2*x.i[0] + 1, 2*x.j[0])
#         union(reg_c, 2*x.j[0] + 1, 2*x.k[0] + 1)
#         union(reg_c, 2*x.k[0],     2*x.l[0] + 1)
#         union(reg_c, 2*x.l[0],     2*x.i[0])
        
#         print("Joining left of " + str(x.j[0])+ " and right of " + str(x.i[0]))
#         print("Joining left of " + str(x.k[0])+ " and right of " + str(x.l[0]))
#         print("Joining left of " + str(x.i[0])+ " and left of " + str(x.l[0]))
#         print("Joining right of " + str(x.j[0])+ " and right of " + str(x.k[0]))
#     else: 
#         print('something is really wrong!')
        
#     print("-"*80)

# for v in hc1.vs:
#     union(reg_c, 2*v.incoming[0][0], 2*v.outgoing[-1][0])
#     union(reg_c, 2*v.incoming[-1][0]+1, 2*v.outgoing[0][0]+1)
    
#     print("Joining left of " + str(v.incoming[0][0])+ " and left of " + str(v.outgoing[-1][0]))
#     print("Joining right of " + str(v.incoming[-1][0])+ " and right of " + str(v.outgoing[0][0]))
    
#     if v.n_incoming > 1:
#         for i in range(v.n_incoming - 1):
#             union(reg_c, 2*v.incoming[i][0]+1, 2*v.incoming[i+1][0])
#             print("Joining left of " + str(v.incoming[i+1][0])+ " and right of " + str(v.incoming[i][0]))
#     if v.n_outgoing > 1:
#         for i in range(v.n_outgoing - 1):
#             union(reg_c, 2*v.outgoing[i][0], 2*v.outgoing[i+1][0]+ 1)
#             print("Joining left of " + str(v.outgoing[i][0])+ " and right of " + str(v.outgoing[i+1][0]))
#     print("-"*80)

# [find(reg_c, i) for i in range(2*hc1.n_seg)]

# nrows = hc1.n_crossings
# ncols = hc1.n_regions

# # use sage to define a symbolic matrix of zeros
# A = matrix(SR, nrows, ncols)

# # first we iterate over the crossings, filling out all of the related entries of the matrix
# # (it is more convenient to make the outer loop go over regions)
# for j, (reg, adj_xvs) in enumerate(hc1.reg_to_xvs.items()):
#     # this is a helper array to check if we have an adjacency
#     adj_xvs_elems_only = [e[0] for e in adj_xvs]
#     for i, x in enumerate(hc1.xs):
#         if x in adj_xvs_elems_only:
#             print(i,j, reg, x)
#             contrib = 0
#             # this covers the possibility of the same crossing being present multiple times 
#             # for the same region 
#             x_present = [e for e in adj_xvs if e[0] == x]
#             for x, pos in x_present:
#                 x_type = '+' if x.sign == 1 else '-'
#                 # extract the color of a proper edge
#                 col = hc1.seg_to_col[x.j[0]]
#                 contrib += hc1.m_Cp(x_type, pos)*hc1.A_Cp(x_type, pos, col)
#             print(contrib)
#             A[i,j] += contrib
#     print("-"*80)
