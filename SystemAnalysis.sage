class SystemAnalysis: 
    def __init__(self,F):
        self.F = Sequence(F)
        self.debug = False
    
    def polynomial_ring(self):
        return self.F.ring()
    
    def coefficient_ring(self):
        return self.F.ring().base_ring()
    
    def variables(self):
        return self.F.variables()[::-1]
    
    def nvariables(self):
        return self.F.nvariables()
    
    def nequations(self):
        return len(self.F)
    
    def nmonomials(self):
        return [len(f.monomials()) for f in self.F]
    
    def lm_degrees(self):
        return [f.lm().degree() for f in self.F]
    
    def partitions(self, K=None):
        '''
        Returns the set of all possible partitions of the variable set.
        
        You may specify the argument K. 
        - If K is an integer, "partitions" returns all partitions of the variable set into k parts.
        - If K is a list of integers (in descending order), "partitions" returns all partitions of 
          the variable set into len(k) groups of sizes given in k.
        
        Returns:
            (SetPartitions): (All) partitions of the variable set.
        '''  
        return SetPartitions(self.F.variables(), K)
        #return SetPartitions(self.polynomial_ring().gens(), K)
    
    def npartitions(self, K=None):
        '''
        Returns number of possible partitions in set of variables of polynomial system self.F.
        
        Returns:
            (int): Number of partitions on variables set. See Bell number.
        '''  
        return self.partitions(K).cardinality()
    
    def is_homogeneous(self):
        return all(f.is_homogeneous() for f in self.F)
    
    def homogeneous_components(self):
        return [f.homogeneous_components() for f in self.F]
    
    def _is_weighted_homogeneous(self, f, normalize=True):
        A = matrix(QQ,f.exponents())  # exponent matrix
        b = column_matrix([QQ.one()]*f.number_of_terms())  #free vector
        Ab = A.augment(b, subdivide=true);  #augmented matrix
        Ab.echelon_form() #echelon form

        # Ab has n+1 columns (for w1,...,wn,w), where n=f.nvariables()
        # Hyperplane is (n+1)-1=n dim subspace of K^(n+1)
        if (Ab.rank() != f.nvariables()):
            return False, None, None

        weights = Ab.echelon_form()[:-1,-1].list()  
        weighted_deg = 1
        # Since w1*a1 + ... + wn*an = 1, for (a1,---,an) being the exponent of a term, all w1,...,wn will be in form 1/...

        # TODO how to know if all weights are positive?
        assert(all([w > 0 for w in weights]))

        # (w1,...,wn,1) -> Make sure w1,...,wn are integers
        if normalize:
            for i in range(len(weights)):
                factor = weights[i].denominator()
                weighted_deg *= factor
                weights = [w * factor  for w in weights]

        return True, weights, weighted_deg
    
    def is_weighted_homogeneous(self, idx=None, normalize=True):
        if idx != None:
            return self._is_weighted_homogeneous(self.F[idx])
        
        return [self._is_weighted_homogeneous(f, normalize) for f in self.F]
    
    def is_m_homogeneous(self, partition=None):
        '''
        Checks if polynomial system self.F is (multi-)homogeneous (w.r.t. given partiton).
        If partition is None, the trivial partition into 1 large group is assumed.
        
        Definition:
            If a polynomial f is multi-homogeneous w.r.t. a partition consisting of l groups, 
            f is called l-homogeneous. In the trivial case, one simply writes "homogeneous"
            instead of 1-homogeneous.

        Parameters:
            partition (list of lists, list of sets, SetPartition): The string w.r.t. which the multidegrees are computed.

        Returns:
            (bool): True if each polynomial in self.F is l-homogeneous, False otherwise.
        '''   
        if partition is None:
            return self.is_homogeneous()
        
        partition = OrderedSetPartition(partition)
        assert(partition in self.partitions())
        
        for f in self.F:
            for j, group in enumerate(partition):
                substitution = {var:f.base_ring().one() for var in f.variables() if var not in group}
                f_ = f.substitute(substitution)
                if not f_.is_homogeneous():
                    return False
        return True
    
    def homogenize(self):
        # return PolynomialSequence
        return Sequence([f.homogenize() for f in self.F])
    
    def m_homogenize(self, partition=None):
        # 1 additional homogeneous coordinate per group in partition
        # If partition is None, 1-homogenization
        # Otherwise: m-homogenization (wrt given partition)
        if partition is None:
            return self.homogenize()
        
        if self.debug:
            print(f"Partition: {partition}")
        
        D, _ = self.m_degrees(partition)
        #print(D)
        (n,l) = D.dimensions()
        
        G = []
        h_vars = list(self.polynomial_ring().variable_names()) + [f"h{j+1}" for j in range(l)]  # One homogenization variable per group.
        
        if self.debug:
            print(f"Variables for {l}-homogenized system: {h_vars}")

        P = PolynomialRing(self.coefficient_ring(), h_vars, order=self.polynomial_ring().term_order())
        #print(P)
        h_vars = [P(v) for v in h_vars]
        
        print('\n'.join([f"{h} -> {set(p)}" for p, h in zip(partition,h_vars[-l:])]))
        
        for i, f in enumerate(self.F):
            g = P(f)
            #print(g)
            for j, group in enumerate(partition):
                h = h_vars[-l+j]
                #print(g.variables())
                substitution = {var:(var/h) for var in g.variables() if var in [P(x) for x in group]}
                #print(substitution)
                g = g.substitute(substitution) # yields FractionFieldElement
                #print(g)
                g = g * h**D[i][j] # multiply out divisors
                #print(h**D[i][j])
                #print(g)
                g = P(g) # convert back to correct ring
                #print("---")
            G.append(P(g))

        if self.debug:
            print(f"{l}-homogenized system wrt partition {partition}:")
            print(Sequence(G,cr=True))
        
        return Sequence(G)
    
    def max_degree(self):
        return self.F.maximal_degree()
    
    def degrees(self):
        '''
        Returns list containing degree of each polynomial in the polynomial system self.F.

        Returns:
            (list): List of degrees, where value at position i is the (maximal) degree of the i-th polynomial in self.F.
        '''   
        return [f.degree() for f in self.F]
    
    def m_degrees(self, partition=None):
        '''
        Returns multidegree matrix and list of group sizes w.r.t. given partiton.
        If partition is None, the trivial partition into 1 large group is assumed.

        Parameters:
            partition (list of lists, list of sets, SetPartition): The string w.r.t. which the multidegrees are computed.

        Returns:
            D (nxl matrix): The multidegree matrix, where D[i,j] is the degree of the i-th polynomial homogenized w.r.t. the j-th set in the given partition.
            K (1xl vector): Group sizes, where K[j] is the size of the j-th group in the partition.
        '''        
        if partition is None:
            D = matrix(ZZ, self.degrees()).transpose()
            K = vector(ZZ, [self.F.nvariables()])
            return D, K 
        
        partition = OrderedSetPartition(partition)
        #assert(partition in self.partitions())

        # Group sizes in partition
        K = vector(ZZ, [len(group) for group in partition])

        # multidegrees |F| x |K| (number of polynomials x number of groups in partition)
        D = matrix(ZZ, len(self.F), len(partition))
        for i, f in enumerate(self.F):
            for j, group in enumerate(partition):
                substitution = {var:f.base_ring().one() for var in f.variables() if var not in group}
                #print(substitution)
                f_ = f.substitute(substitution)
                # when std_grading is not set, degree of zero polynomial is -1 -> does not work for algebraic closure
                #D[i,j] = f_.degree(std_grading=True) # TODO: check if this has downsides
                if f_ == f.parent().zero():
                    D[i,j] = 0
                else: 
                    D[i,j] = f_.degree() # TODO: check if this has downsides
        
        return D, K
    
    def bezout_bound(self):
        '''
        Returns the Bézout bound of the polynomial system self.F.
        
        Assumption: 
            Number of sulutions is finite.
        
        Bézout theorem:
            n homogeneous polynomials of degree d1, ..., dn in n + 1 indeterminates define either an 
            algebraic set of positive dimension, or a zero-dimensional algebraic set consisting of 
            d1 * ... * dn points counted with their multiplicities. 
        
        Returns:
            b (int): The Bézout bound of the system self.F.
        '''
        return prod(self.degrees())
    
    def mh_bezout_test_monomial(self, P, K):
        assert(P.ngens() == len(K))
        return prod(list(map(lambda tj, nj : tj**nj, P.gens(), K)))
        
    def mh_bezout_linear_form(self, P, D):
        assert(P.ngens() == D.ncols()) # One indeterminate for each group
        lf = P.one()
        for mdf in D.rows(): # mdf = multidegree of function f
            lf = lf * sum(list(map(lambda dj, tj : dj * tj, mdf, P.gens())))
        return lf
        
    def mh_bezout_coeff_linear_form(self, D, K):
        (n,l) = D.dimensions()
        P = PolynomialRing(ZZ, [f't{i}' for i in range(l)])
        
        testMonomial = self.mh_bezout_test_monomial(P, K)
        if self.debug:
            print(f"Test monomial: {testMonomial}")
        
        testPolynomial = self.mh_bezout_linear_form(P, D)
        if self.debug:
            print(f"Product of linear forms: {testPolynomial}")
        
        return testPolynomial.monomial_coefficient(testMonomial) # 0 if no intersections (or only 1 poly given), Bezout bound otherwise
    
    def row_expansion(self, D, K, i):
        # TODO memory efficient version?
        
        (n,l) = D.dimensions()
        assert(len(K) == l)
        
        if i == n:
            return 1

        b = 0
        for j in range(l):
            if K[j] != 0 and D[i,j] != 0:
                K_ = copy(K)
                K_[j] -= 1
                b = b + D[i,j] * self.row_expansion(D, K_, i + 1)     
        return b
    
    def mh_bezout_bound(self, partition=None, algorithm='rowExpansion'):
        '''
        Returns the multihomogeneous Bézout bound of the polynomial system self.F w.r.t. given partiton.
        If partition is None, the trivial partition into 1 large group is assumed implicitely 
        (by the called algorithms).
        
        Assumption: 
            Number of sulutions is finite.
        
        Multi-homogeneous Bézout theorem:
            Let n = n1 + ... + nl. n multi-homogeneous polynomials (in n+l variables) of multi-degrees 
            d1, ..., dn (rows of D) define either a multi-projective algebraic set of positive dimension, 
            or a zero-dimensional algebraic set consisting of b points, counted with multiplicities, 
            where b is the coefficient of t1^n1 * ... * tl^nl in the product of linear forms d1 * ... * dn.

        Parameters:
            partition (list of lists, list of sets, SetPartition): The string w.r.t. which the multi-homogeneous Bézout bound is computed.
            algorithm ('rowExpansion', 'bezoutTheorem'): Algorithm that is used for the computation of the multi-homogeneous Bézout bound.

        Returns:
            b (int): The multihomogeneous Bézout bound of the system self.F.
        '''
        if algorithm == 'rowExpansion':
            D, K = self.m_degrees(partition)
            b = self.row_expansion(D, K, 0) # TODO maybe switch rows to start expansion with row with max number of zeros
        elif algorithm == 'bezoutTheorem':
            D, K = self.m_degrees(partition)
            b = self.mh_bezout_coeff_linear_form(D, K)
        else:
            print("Undefined algorithm")
            b = None
            
        return b
    
    def min_mh_bezout_bound(self, max_npartitions=50000, K=None):
        # Partitions to consider
        partitions = self.partitions(K)
        npartitons = partitions.cardinality()
        
        min_mhb = +Infinity
        min_partitions = []
        mhb_dict = {} # for plot
        
        if self.debug:
            print(f"Considered partitions: {partitions}")
            print(f"Number of partitions: {npartitons}")
        if npartitons > max_npartitions:
            print(f"Number of partions is larger than proposed max ({max_npartitions}). Consider increasing the default maximum or analyzing only a subset of all possible partitions.")
            return min_mhb, min_partitions, mhb_dict
        
        for partition in partitions:
            mhb = self.mh_bezout_bound(partition)
            # Count for plot
            if mhb in mhb_dict:
                mhb_dict[mhb] += 1
            else:
                mhb_dict[mhb] = 1
            # Comparison
            if mhb < min_mhb:
                min_mhb = mhb
                min_partitions = [partition]
            elif mhb == min_mhb:
                min_partitions.append(partition)
        
        if self.debug:
            print(f"Minimal multi-homogeneous Bézout bound: {min_mhb}")
            print("Minimal partitions:")
            for i, partition in enumerate(min_partitions):
                print(i, partition)
        
        # plot
        #import matplotlib.pyplot as plt
        #import numpy as np
        
        #plt.bar(list(mhb_dict.keys()), mhb_dict.values(), color='g')
        #plt.bar([ str(i) for i in mhb_dict.keys()], mhb_dict.values(), color='g')
        
        #t = [min(mhb_dict.keys()), max(mhb_dict.keys())]
        #plt.xticks(t,t)
        #plt.xticks(fontsize=45, rotation=90)
        
        #plt.title("MH Bézout bounds")
        #plt.xlabel("Bound")
        #plt.ylabel("Frequency")
        #plt.show()
        
        #print(f"Max bound: {max(mhb_dict.keys())}")
        
        return min_mhb, min_partitions, mhb_dict