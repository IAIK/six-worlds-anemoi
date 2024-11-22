load('anemoi.sage')
load('models.sage')

import sys
from time import time


#=============================================================================================
def main(P, model, final_ll, varord, steps, rmax, tmax):
    """    
    Description
    
    Inputs:
    P            Specific instance of Anemoi permutation.
    model        Function used to model specific instance.
    final_ll     Use final linear layer in model.
    varord       Variable ordering (as defined in models.sage or None for automatic selection)
    steps        1: GB, 2: GB + FGLM, 3: GB + FGLM + ELIM. 4: GB + FGLM + ELIM + multiplicities.
    rmax         Maximum number of rounds.
    tmax         Maximum time [s] for computations.
    
    MAGMA algorithms:
    
    GroebnerBasis():
        Explicitly force a Gröbner basis (GB) for the ideal I to be constructed. 
        The internally saved Gröbner basis in I is updated.
        Returns GröbnerBasis wrt ordering of I and degrees from during computation.
        
        The parameter "algorithm" may be set to one of: "Default", "Direct", "FGLM" or "Walk".
        
        The value "Direct" specifies that Magma should compute the GB of I (with respect
        to the order of I) by a direct algorithm alone, so that an order-conversion algorithm
        is not used (the parameter Faugere below controls which direct algorithm is used).
        
        Two direct algorithms are available:
        (1) The Faugère F4 algorithm (faugere='true'), which works by specialized sparse linear al-
        gebra and is applicable to ideals defined over a finite field or the rational field;
        (2) The Buchberger algorithm (faugere='false') for ideals defined over any field.
    
    ChangeOrder():
        Given an ideal I of the polynomial ring P = R[x1, . . . , xn], together with a monomial
        order order, construct the polynomial ring Q = R[x1, . . . , xn] with order order, and 
        then return the ideal J of Q corresponding to I and the isomorphism f from P to Q. 
         
        The point of the function is that one can change the order on monomials of I to
        be that of Q. When a Gr ̈obner basis of J is needed to be calculated, Magma uses
        a conversion algorithm starting from a Gröbner basis of I if possible—this usually
        makes order conversion much more efficient than by computing a Gr ̈obner basis of
        J from scratch.

    Variety():
        Given a zero-dimensional ideal I of a polynomial ring P , return the variety of I over
        its coefficient field K as a sequence of tuples. Each tuple is of length n, where n is
        the rank of P , and corresponds to an assignment of the n variables of P (in order)
        such that all polynomials in I vanish with this assignment.
        
        The function works in the zero-dimensional case by first computing a triangular
        or radical decomposition of I. This reduces the problem to successively computing roots 
        of univariate polynomials.
    """
    
    M = Magma()
    
    for nr in range(1,rmax+1):
        print("="*100, flush=True)
        
        # Iterate over number of rounds
        system, _ = model(P, n_rounds=nr, final_ll=final_ll, ordering=varord, debug=False)

        # Print system information
        R = system[0].parent() # the polynomial ring with designated variables
        print(f"n_r: {nr}, n_v: {len(R.variable_names())}, n_e: {len(system)}, mdeg: {max([f.degree() for f in system])}", flush=True)
        print(f"degs = {[f.degree() for f in system]}", flush=True)
        print(R, flush=True)
        
        # Define polynomial ring (implicitely defines base-ring/-field)
        r = M(R)
        M.eval(f"r := {r.name()};")
            
        # Define system
        s = M(system)
        M.eval(f"s := {s.name()};")
        
        # Define ideal
        M.eval(f"i := ideal<r|s>;")
        
        # Magma GröbnerBasis computation
        M.eval("SetNthreads(5);")
        M.eval("SetVerbose(\"Groebner\", 1);")
        M.eval("SetVerbose(\"FGLM\", 3);")
        
        # Step 1: degrevlex Gröbner basis computation
        if steps >= 1:
            print("-"*100, flush=True)
            print(f"Starting GB computation... (r={nr})", flush=True)
            
            time_gb = time()
            out = M.eval(f"time gb, dregs := GroebnerBasis(i : Al := \"Direct\", Faugere := true);") # F4 algorithm
            time_gb = time() - time_gb
            print(out, flush=True) # this includes log of F4 algorithm and timing information
            
            out = M.eval("print dregs")
            print(f"dregs = {out}", flush=True)

            if time_gb > tmax:
                break
        
        # Step 2: FGLM change of order algorithm: degrevlex -> lex
        if steps >= 2:
            print("-"*100, flush=True)
            print(f"Starting FGLM computation... (r={nr})", flush=True)

            M.eval(f"i := ChangeOrder(i,\"lex\");")
            time_fglm = time()
            out = M.eval(f"time gb_lex := GroebnerBasis(i: Al := \"FGLM\");")
            time_fglm = time() - time_fglm
            print(out, flush=True) # this includes log of FGLM algorithm and timing information
            
            out = M.eval(f"print Degree(gb_lex[#gb_lex]);")
            print(f"unideg = {out}", flush=True)
            
            # If FGLM exceeds a certain time limit and it is the last step in this round
            if (steps == 2) and (time_fglm > tmax):
                # Switch to computing GB only from the next round on
                steps = 1
        
        # Step 3: Compute the affine variety
        if steps >= 3:
            print("-"*100, flush=True)
            print(f"Starting UNIV SOLVING computation... (r={nr})", flush=True)
            
            time_elim = time()
            out = M.eval(f"time v := Variety(i);")
            # out = M.eval(f"time v := Variety(i, AlgebraicClosure());")
            time_elim = time() - time_elim
            print(out, flush=True)
            
            nsols = M("#v;").sage()
            print(f"nsols = {nsols}", flush=True)
                
            # For Anemoi, this step is pretty fast, no need to stop here if lex basis was already computed.
        
        # Compute the sum of multiplicities for the points in the affine variety
        if steps >= 4:
            print("-"*100, flush=True)
            print(f"Starting MULTIPLICITY computation... (r={nr})", flush=True)
            
            M.eval("get_multiplicity := func<pt | Dimension(quo<r | gb cat [r.j - pt[j] : j in [1..#pt]]>)>;")
            time_mult = time()
            out = M.eval(f"time nsols_mult := (#v eq 0) select 0 else &+ [get_multiplicity(pt) : pt in v];")
            time_mult = time() - time_mult
            #print(out, flush=True)  # Triggers GB again (but faster)
            
            nsols_mult = M("nsols_mult;").sage()
            print(f"nsols_mult = {nsols_mult}", flush=True)
            
            # If multiplicity computation exceeds a certain time limit
            if time_mult > tmax:
                # Switch to not doing this step from the next round on
                steps = 3
            if time_fglm > tmax:
                # FGLM was too expensive before
                steps = 1
                        


#=============================================================================================
if __name__ == "__main__":
    """    
    Inputs:
    p            Prime.
    n            Power for q=p^n, where n=1 for p odd.
    l            Input size to permutation is 2l. For our purpose: l=1.
    alpha        Alpha in power function.
    steps        1: GB, 2: GB + FGLM, 3: GB + FGLM + ELIM. 4: GB + FGLM + ELIM + multiplicities.
    rmax         Maximum number of rounds.
    tmax         Maximum time [s] for computations.
    model        Used algebraic model. (FCICO, PCICO)
    final_ll     Include final linear layer.
    varord       Variable ordering (as defined in models.sage) (optional)
    QUAD         Exponent for polynomials Q to be used in Flystel (optional)
    """
    args = sys.argv
    assert(len(args) in [1+9,1+10,1+11])  # varord and QUAD arguments are optional
    
    ANEMOI_PRIMES = {
        'BLS12-381': BLS12_381_SCALARFIELD,
        'BN-254': BN_254_SCALARFIELD
    }
    
    try:
        p = int(args[1])
    except:
        try:
            p = ANEMOI_PRIMES[args[1]]
        except:
            print('Invalid prime')
            exit(-1)
        
    n = int(args[2])
    l = int(args[3])
    alpha = int(args[4])
    
    # Attack steps to be performed. Default: 4
    steps = int(args[5])
    steps = steps if steps in [1,2,3,4] else 4
    
    # Maximum number of rounds. Default: 10
    rmax = int(args[6])
    rmax = rmax if rmax > 0 else 10
    
    # Maximum time [s] for individual computation steps. Default: 3600s (1h)
    tmax = float(args[7])
    tmax = tmax if tmax > 0 else 3600
    
    model = args[8]
    final_ll = args[9].capitalize() == 'True'
    varord = int(args[10]) if len(args) > 10 else None
    QUAD = int(args[11]) if len(args) > 11 else None
    print(f"Model: {model}. Use final LL: {final_ll}.")
    
    
    # Check arguments
    if (p < 2) or (n < 1) or (l < 1) or (alpha < 3) or (steps < 0) or (rmax < 0) or (tmax < 0):
        print("Input arguments not correct.")
        exit(-1)
    
    # Check cases
    if p == 2 and is_odd(n):
        # alpha = 2^i + 1
        i = log(alpha-1,2)
        if not isinstance(i, Integer):
            print(f"alpha = {alpha} invalid for field of size {p}^{n}.")
            exit(-1)
    elif p > 2 and is_prime(p):
        assert(n==1)
        if gcd(alpha, p-1) != 1:
            print(f"alpha = {alpha} invalid for field of size {p}^{n}.")
            exit(-1)
    else:
        print(f"Invalid field size {p}^{n}.")
        exit(-1)
        
    P = AnemoiPermutation(q=p^n, alpha=alpha, n_rounds=rmax, n_cols=l, QUAD=QUAD)
    P.alpha = alpha # alpha=3 hard coded in anemoi code for characteristic 2
    
    # Anemoi information
    print("="*100, flush=True)
    print(f"ANEMOI GF({p}^{n}), alpha = {alpha}, QUAD = {P.QUAD}. {'With ' if final_ll else 'No '}final LL. Model: {model}", flush=True)
    print(f"rmax = {rmax}, tmax = {tmax}.", flush=True)
    print("="*100, flush=True)
    
    # Run attack using specific model
    if model == 'FCICO':
        main(P, model_F_CICO, final_ll, varord, steps, rmax, tmax)
    elif model == 'PCICO':
        main(P, model_P_CICO, final_ll, varord, steps, rmax, tmax)
    else:
        print(f"Unknown model: {model}.")
        exit(-1)

