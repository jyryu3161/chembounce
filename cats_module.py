# source code from https://github.com/iwatobipen/CATS2D/blob/master/cats2d/rd_cats2d.py

from rdkit import Chem
import numpy as np
from scipy.sparse.csgraph import shortest_path


cats_smarts = {
    "D": ["[!$([#6,H0,-,-2,-3])]"],
    "A": ["[!$([#6,F,Cl,Br,I,o,s,nX3,#7v5,#15v5,#16v4,#16v6,*+1,*+2,*+3])]"],
    "P": ["[*+]", "[#7H2]"],
    "N": ["[*-]", "[C&$(C(=O)O),P&$(P(=O)),S&$(S(=O)O)]"],
    "L": [
        "[Cl,Br,I]",
        "[S;D2;$(S(C)(C))]",
        "[C;D2;$(C(=C)(=C))]",
        "[C;D3;$(C(=C)(C)(C))]",
        "[C;D4;$(C(C)(C)(C)(C))]",
        "[C;D3;H1;$(C(C)(C)(C))]",
        "[C;D2;H2;$(C(C)(C))]",
    ],
    "R": ["[a]"],
}

cats_desc = [
    "DD",
    "AD",
    "DP",
    "DN",
    "DL",
    "DR",
    "AA",
    "AP",
    "AN",
    "AL",
    "AR",
    "PP",
    "NP",
    "LP",
    "PR",
    "NN",
    "LN",
    "NR",
    "LL",
    "LR",
    "RR",
]


class CATS2D:
    def __init__(
        self, max_path_len=9, verbose=False, scale_type="raw", cats_smarts=cats_smarts, cats_desc=cats_desc
    ):
        self.max_path_len = max_path_len
        self.verbose = verbose
        self.scale_type = scale_type
        self.cats_smarts = cats_smarts
        self.cats_desc = cats_desc
        if scale_type not in ("raw", "num", "occ") or max_path_len < 0:
            raise ValueError("Invalid CATS scaling or path length")
        self._patterns = {k: [Chem.MolFromSmarts(p) for p in v] for k, v in cats_smarts.items()}
        if any(pattern is None for patterns in self._patterns.values() for pattern in patterns):
            raise ValueError("Invalid CATS SMARTS")
        if any(len(pair) != 2 or any(label not in cats_smarts for label in pair) for pair in cats_desc):
            raise ValueError("Invalid CATS pair definitions")
        self._pair_index = {"".join(sorted(k)): i for i, k in enumerate(cats_desc)}

    def getPcoreGroups(self, mol):
        """
        Given a molecule it assigns PPP's to individual atoms

        The return value is a list of length number of atoms in
        the input molecule. The i'th element of the list contains
        a list of PPP labels that were identified for the i'th atom
        """

        labels = [set() for _ in mol.GetAtoms()]
        for label, patterns in self._patterns.items():
            for pattern in patterns:
                for match in mol.GetSubstructMatches(pattern):
                    for idx in match:
                        labels[idx].add(label)
        return [sorted(group) for group in labels]

    def _getZeroMatrix(self, r, c):
        return [[0 for x in range(0, c)] for x in range(0, r)]

    def getAdjacencyMatrix(self, mol):
        """
        Generates an adjacency matrix for a molecule. Note that this
        will be a matrix with 0 along the diagonals
        """
        n = mol.GetNumAtoms()
        admat = self._getZeroMatrix(n, n)
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
        distances = shortest_path(np.asarray(admat), directed=False, unweighted=True)
        return np.where(np.isfinite(distances), distances, 99999999).astype(int).tolist()

    def getPPPMatrix(self, admat, ppp_labels):
        pppm = {}
        n = len(admat)
        for i in range(0, n):
            ppp_i = ppp_labels[i]
            if ppp_i == "":
                continue
            for j in range(0, n):
                ppp_j = ppp_labels[j]
                if ppp_j == "":
                    continue
                pairs = []
                for x in ppp_i:
                    for y in ppp_j:
                        if (x, y) not in pairs and (y, x) not in pairs:
                            ## make sure to add the labels in increasing
                            ## lexicographical order
                            if x < y:
                                tmp = (x, y)
                            else:
                                tmp = (y, x)
                            pairs.append(tmp)
                pppm[i, j] = pairs
        return pppm

    def getCATs2D(self, mol):
        natom = mol.GetNumAtoms()
        if not natom:
            raise ValueError("CATS requires a nonempty molecule")
        labels = self.getPcoreGroups(mol)
        distances = Chem.GetDistanceMatrix(mol).astype(int)
        within = distances <= self.max_path_len
        features = {label: np.array([label in group for group in labels]) for label in self.cats_smarts}
        descriptor = []
        for left, right in self.cats_desc:
            a, b = features[left], features[right]
            pairs = np.logical_and.outer(a, b)
            if left != right:
                pairs |= np.logical_and.outer(b, a)
            counts = np.bincount(distances[pairs & within], minlength=self.max_path_len + 1).astype(float)
            if self.scale_type == "num":
                counts /= natom
            elif self.scale_type == "occ":
                divisor = a.sum() + b.sum()
                if divisor:
                    counts /= divisor
            descriptor.extend(counts.tolist())
        return descriptor

    def getCATS2Ddesc(self):
        res = []
        for label in self.cats_desc:
            for i in range(0, self.max_path_len + 1):
                res.append("%s.%d " % (label, i))
        return res
