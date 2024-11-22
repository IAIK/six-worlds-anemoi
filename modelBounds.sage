# TODO those statements only hold for prime field!

# ----------------------------------------------------------------------------------
# Model details
# ----------------------------------------------------------------------------------

fun_nv_FCICO = lambda N: 2*N
fun_nv_PCICO = lambda N: N + 1

fun_ne_FCICO = lambda N: 2*N
fun_ne_PCICO = lambda N: N + 1

fun_degs_FCICO = lambda N, alpha: flatten([[alpha, alpha] for _ in range(1, N+1)])
fun_degs_PCICO = lambda N, alpha: [max(2*r, alpha) for r in range(1, N+1)] + [N+1]

# ----------------------------------------------------------------------------------
# Extrapolations
# ----------------------------------------------------------------------------------

fun_dreg_approx = {'FCICO': [lambda N, alpha: (alpha + 1) / 2 * (N+1),
                            ], 
                   'PCICO': [lambda N, alpha: (alpha + 3)/2 * N + (alpha-1)/2,
                            ]}
title_dreg_approx = {'FCICO': ['Approx', 'dexp AP'], 'PCICO': ['Approx','Interpol']}

fun_dI_FCICO = lambda N, alpha : (alpha + 2)^N
fun_dI_PCICO = fun_dI_FCICO

fun_udeg_FCICO = fun_dI_FCICO
fun_udeg_PCICO = fun_udeg_FCICO

# ----------------------------------------------------------------------------------
# Theoretical bounds
# ----------------------------------------------------------------------------------

# --------------------------------------
# Macaulay bound
# --------------------------------------

fun_mac_FCICO = lambda N, alpha: sum(map(lambda d: d - 1, fun_degs_FCICO(N,alpha))) + 1
fun_mac_PCICO = lambda N, alpha: sum(map(lambda d: d - 1, fun_degs_PCICO(N,alpha))) + 1

# --------------------------------------
# Bézout bound
# --------------------------------------

fun_b_FCICO = lambda N, alpha: prod(fun_degs_FCICO(N,alpha))
fun_b_PCICO = lambda N, alpha: prod(fun_degs_PCICO(N,alpha))

# --------------------------------------
# Multihomogeneous Bézout bound
# --------------------------------------

# For FCICO, "minimal" multihomogenoues Bézout equals Bézout bound

fun_mhb_FCICO = fun_b_FCICO 

# Derived "minimal" multihomogenoues Bézout bound for PCICO

# For fixed values in the case N < r_alpha
# TODO update! first value might be for N = 2
mhb_dict = {
    3: [6, 36,252,1764,12348,86436,605052,4235364],
    5: [9, 75,600,5400,48600,437400,3936600,35429400],
    7: [11, 121,1331,13720,150920,1660120,18261320,200874520],
    11: [15, 225,3375,50625,759375,11390625,170859375,2562890625]
}

# First round r for which 2*r > alpha (poi = point of interest = r_alpha)
def poi(alpha):
    r = 1
    while(2*r < alpha):
        r += 1
    return r

# Given that alpha is always odd, we have the following
r_alpha = lambda alpha: Integer((alpha + 1) / 2)
assert(all([poi(a) == r_alpha(a) for a in [3,5,7,9,11]]))

def t_alpha(a): 
    assert(a in [3,5,7,11])
    ra = r_alpha(a)
    return (a + 4)^(ra) if a == 11 else 2*ra*a^(ra - 1)*(ra + 1)

fun_mhb_PCICO = lambda N, alpha: mhb_dict[alpha][N-1] if N < r_alpha(alpha) else t_alpha(alpha) * (alpha + 4)^(N - r_alpha(alpha))
