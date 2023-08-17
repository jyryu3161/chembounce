# source code from https://github.com/iwatobipen/CATS2D/blob/master/cats2d/rd_cats2d.py

from rdkit import Chem
import numpy as np


cats_smarts = {
    'D' : ['[!$([#6,H0,-,-2,-3])]'],
    'A' : ['[!$([#6,F,Cl,Br,I,o,s,nX3,#7v5,#15v5,#16v4,#16v6,*+1,*+2,*+3])]'],
    'P' : ['[*+]', '[#7H2]'],
    'N' : ['[*-]', '[C&$(C(=O)O),P&$(P(=O)),S&$(S(=O)O)]'],
    'L' : ['[Cl,Br,I]', '[S;D2;$(S(C)(C))]',
       '[C;D2;$(C(=C)(=C))]', '[C;D3;$(C(=C)(C)(C))]', '[C;D4;$(C(C)(C)(C)(C))]',
       '[C;D3;H1;$(C(C)(C)(C))]', '[C;D2;H2;$(C(C)(C))]'],
    'R' : ['[a]']
}
    
cats_desc = ['DD','AD', 'DP', 'DN', 'DL', 'DR',
         'AA', 'AP', 'AN', 'AL', 'AR',
         'PP', 'NP', 'LP', 'PR',
         'NN', 'LN', 'NR',
         'LL', 'LR',
         'RR']

class CATS2D():

    def __init__(self, max_path_len=9, verbose=False, scale_type='raw', cats_smarts=cats_smarts, cats_desc=cats_desc):
        self.max_path_len = max_path_len
        self.verbose = verbose
        self.scale_type = scale_type
        self.cats_smarts = cats_smarts
        self.cats_desc = cats_desc
    
    def getPcoreGroups(self, mol):
        """
        Given a molecule it assigns PPP's to individual atoms
    
        The return value is a list of length number of atoms in
        the input molecule. The i'th element of the list contains
        a list of PPP labels that were identified for the i'th atom
        """
    
        ret = ['' for x in range(0,mol.GetNumAtoms())]
    
        labels = self.cats_smarts.keys()
        for label in labels:
            patterns = self.cats_smarts[label]
            for pattern in patterns:
             
                patt = Chem.MolFromSmarts(pattern)
                matched = False
                for matchbase in mol.GetSubstructMatches(patt, uniquify=True):
                    for idx in matchbase:
                        if ret[idx] == '': ret[idx] = [label]
                        else: 
                            tmp = ret[idx]
                            tmp.append(label)
                            ret[idx] = tmp
                    matched = True
                if matched: break
        return ret
    
    def _getZeroMatrix(self, r,c):
        return [[0 for x in range(0,c)] for x in range(0,r)]
    
    def getAdjacencyMatrix(self, mol):
        """
        Generates an adjacency matrix for a molecule. Note that this
        will be a matrix with 0 along the diagonals
        """
        n = mol.GetNumAtoms()
        admat = self._getZeroMatrix(n,n)
        for bond in mol.GetBonds():
            bgn_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            admat[bgn_idx][end_idx] = 1
            admat[end_idx][bgn_idx] = 1
        return admat
    
    def getTopologicalDistanceMatrix(self, admat):
        """
        Generates the topological distance matrix given
        an adjacency matrix
        """
        n = len(admat)
        d = self._getZeroMatrix(n,n)
        for i in range(0,n):
            for j in range(0,n):
                if admat[i][j] == 0: d[i][j] = 99999999;
                else: d[i][j] = 1
        for i in range(0,n): d[i][i] = 0
        
        for k in range(0,n):
            for i in range(0,n):
                for j in range(0,n):
                    if d[i][k]+d[k][j] < d[i][j]: d[i][j] = d[i][k]+d[k][j]
        return d
    
    def getPPPMatrix(self, admat, ppp_labels):
        pppm = {}
        n = len(admat)
        for i in range(0,n):
            ppp_i = ppp_labels[i]
            if ppp_i == '': continue
            for j in range(0,n):
                ppp_j = ppp_labels[j]
                if ppp_j == '': continue
                pairs = []
                for x in ppp_i:
                    for y in ppp_j:
                        if (x,y) not in pairs and (y,x) not in pairs:
                            ## make sure to add the labels in increasing
                            ## lexicographical order
                            if x < y: tmp = (x,y)
                            else: tmp = (y,x)
                            pairs.append(tmp)
                pppm[i,j] = pairs
        return pppm
    
    def getCATs2D(self, mol):
        natom = mol.GetNumAtoms()
    
        ppp_labels = self.getPcoreGroups(mol)
        admat = self.getAdjacencyMatrix(mol)
        tdistmat = self.getTopologicalDistanceMatrix(admat)
        pppmat = self.getPPPMatrix(tdistmat, ppp_labels)
    
        # get the occurence of each of the PPP's
        ppp_count = dict(zip(['D','N', 'A', 'P', 'L', 'R'], [0]*6))
        for label in ppp_labels:
            for ppp in label:
                ppp_count[ppp] = ppp_count[ppp] + 1
    
        # lets calculate the CATS2D raw descriptor
        desc = [ [0 for x in range(0, self.max_path_len+1)] for x in range(0,21)]
        for x,y in pppmat.keys():
            labels = pppmat[x,y]
            dist = tdistmat[x][y]
            if dist > self.max_path_len: continue
            for pair in labels:
                id = '%s%s' % (pair[0],pair[1])
                idx = cats_desc.index(id)
                vals = desc[idx]
                vals[dist] += 1
                desc[idx] = vals
    
        if self.scale_type == 'num':
            for row in range(0, len(desc)):
                for col in range(0, len(desc[0])):
                    desc[row][col] = float(desc[row][col]) / natom
        elif self.scale_type == 'occ':
            #  get the scaling factors
            facs = [0]*len(cats_desc)
            count = 0
            for ppp in cats_desc:
                facs[count] = ppp_count[ppp[0]] + ppp_count[ppp[1]]
                count += 1
    
            # each row in desc corresponds to a PPP pair
            # so the scale factor is constant over cols of a row
            count = 0
            for i in range(0, len(desc)):
                if facs[i] == 0: continue
                for j in range(0, len(desc[0])):
                    desc[i][j] = desc[i][j] / float(facs[i])
    
        res = []
        for row in desc:
            for col in row:
                res.append(col)
        return res

    def getCATS2Ddesc(self):
        res = []
        for label in self.cats_desc:
            for i in range(0, self.max_path_len+1):
                res.append('%s.%d ' % (label,i))
        return res