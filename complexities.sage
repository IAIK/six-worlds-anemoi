OMEGAS = [2,2.37,3]

# -----------------------------------------------------------------------------------------------------
# Step 1: GB complexity (different formulas found in literature)
# -----------------------------------------------------------------------------------------------------
def gb_complexity_1(ne, nv, dreg, w, bitcompl=True, degs=None):  # formula used by Faugere
    c = nv * binomial(ne + dreg, dreg)**w
    return floor(log(c,2)) if bitcompl else c

def gb_complexity_2(ne, nv, dreg, w, bitcompl=True, degs=None):  # Formula used in Anemoi paper
    c = binomial(nv + dreg, nv)**w
    return floor(log(c,2)) if bitcompl else c

def gb_complexity_3(ne, nv, dreg, w, bitcompl=True, degs=None):  # Formula in Anemoi paper, including factor nv
    c = nv * binomial(nv + dreg, nv)**w
    return floor(log(c,2)) if bitcompl else c

def gb_complexity_4(ne, nv, dreg, w, bitcompl=True, degs=None):  # Formula used in Spaenlehauer PhD Thesis
    c = 0
    for i in IntegerRange(dreg+1):
        c_ = 0
        for j in IntegerRange(1,ne+1):
            deg_fj = degs[j-1]
            c_ += binomial(nv + i - deg_fj - 1, i - deg_fj)
        c += binomial(nv + i - 1, i)**(w-1) * c_  # full rank case, HF = 0 in formula
    return floor(log(c,2)) if bitcompl else c

# Used in our paper
def gb_complexity(ne, nv, dreg, w, bitcompl=True, degs=None, fun_gb_complexity=gb_complexity_4):
    """Return complexity of GB computation step.

    ne        Number of equations.
    nv        Number of variables.
    dreg      Degree of regularity (either from experiments or approximation / bound).
    w         Linear algebra 'constant'.
    bitcompl  Return log(c,2) for complexity c calculated from formula if true, else c.
    degs      List of polynomial degrees in system. Used only by gb_complexity_4.
    """
    if dreg is None or dreg == -1:
        return None
    assert(degs is None or len(degs) == ne)
    dreg = ceil(dreg) # extrapolated dreg might be real, not integer
    return fun_gb_complexity(ne, nv, dreg, w, bitcompl, degs)

# -----------------------------------------------------------------------------------------------------
# Step 2: FGLM complexity
# -----------------------------------------------------------------------------------------------------
def fglm_complexity(nv, dI, w, bitcompl=True): 
    """Return complexity of FGLM computation step.

    nv        Number of variables.
    dI        Quotient space dimension (either from experiments or approximation / bound).
    w         Linear algebra 'constant'.
    bitcompl  Return log(c,2) for complexity c calculated from formula if true, else c.
    """
    if dI is None or dI == -1:
        return None
    c = nv * dI**w
    return floor(log(c,2)) if bitcompl else c

# -----------------------------------------------------------------------------------------------------
# Step 3: FAC complexity (for GBlex in shape position, factorization of last univariate polynomial)
# -----------------------------------------------------------------------------------------------------
def fac_complexity(d, bitcompl=True): 
    """Return complexity of FAC computation step.

    nv        Number of variables.
    dI        Quotient space dimension (either from experiments or approximation / bound).
    w         Linear algebra 'constant'.
    bitcompl  Return log(c,2) for complexity c calculated from formula if true, else c.
    """
    if d is None or d == -1:
        return None
    c = d**1.815
    return floor(log(c,2)) if bitcompl else c

# -----------------------------------------------------------------------------------------------------
# Find minimum round number: N* >= min{r : C_alg(N) >= 2^s}
# -----------------------------------------------------------------------------------------------------
def find_opt_num_rounds_gb(fun_ne, fun_nv, fun_dreg, fun_degs, alpha, w, s):
    """Return round number derived from FGLM complexity: N* >= min{r : C_FGLM(N) >= 2^s}
    
    fun_ne      Function to compute number of equations. fun_ne: N -> ne.
    fun_nv      Function to compute number of variables. fun_nv: N -> nv.
    fun_dreg    Function to compute degree of regularity (either approximation or bound). fun_dreg: (N,alpha) -> dreg.
    fun_degs    Function to compute list of degrees of polynomial equation system. fun_degs: (N,alpha) -> [...].
    alpha       Degree of permutation E(x) = x^alpha used in Flystel construction.
    w           Linear algebra 'constant'.
    s           Target security level.
    """
    N_opt = 1
    while(gb_complexity(ne=fun_ne(N_opt), nv=fun_nv(N_opt), dreg=fun_dreg(N_opt,alpha), 
                        w=w, bitcompl=False, degs=fun_degs(N_opt,alpha)) < 2^s):
        N_opt += 1
    return N_opt


def find_opt_num_rounds_fglm(fun_nv, fun_dI, alpha, w, s):
    """Return round number derived from FGLM complexity: N* >= min{r : C_FGLM(N) >= 2^s}
    
    fun_nv    Function to compute number of variables. fun_nv: N -> nv.
    fun_dI    Function to compute quotient space dimension (either approximation or bound). fun_dI: (N,alpha) -> dI.
    alpha     Degree of permutation E(x) = x^alpha used in Flystel construction.
    w         Linear algebra 'constant'.
    s         Target security level.
    """
    N_opt = 1
    while(fglm_complexity(nv=fun_nv(N_opt), dI=fun_dI(N_opt,alpha), w=w, bitcompl=False) < 2^s):
        N_opt += 1
    return N_opt


def find_opt_num_rounds_fac(fun_udeg, alpha, s):
    """Return round number derived from FAC complexity: N* >= min{r : C_FGLM(N) >= 2^s}
    
    fun_udeg  Function to compute univariate degree (either approximation or bound). fun_udeg: (N,alpha) -> d.
    alpha     Degree of permutation E(x) = x^alpha used in Flystel construction.
    s         Target security level.
    """
    N_opt = 1
    while(fac_complexity(d=fun_udeg(N_opt,alpha), bitcompl=False) < 2^s):
        N_opt += 1
    return N_opt
