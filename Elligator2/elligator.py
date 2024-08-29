from core import *

A = Mt.A
B = GF(1)
# Note that B is not the same coeficient as Montgomery curves.
# Montgomery curves always set this particular coeficient to 1.


############################
# Reference Implementation #
############################
def dir_map_ref(r):
    """Computes a point (u, v) from the representative r in GF(p)

    Always succeds
    """
    w = -A / (GF(1) + Z * r**2)
    e = legendre(w**3 + A*w**2 + w)
    u = e*w - (GF(1)-e)*(A//2)
    v = -e * sqrt(u**3 + A*u**2 + B*u)
    return (u, v)

def rev_map_ref(u, v_is_negative):
    """Computes the representative of the point (u, v), if possible

    Returns None if the point cannot be mapped.
    """
    if u == -A or not is_square(-Z * u * (u+A)):
        return None
    sq1 = sqrt(-u     / (Z * (u+A)))
    sq2 = sqrt(-(u+A) / (Z * u    ))
    rep = sq1
    rep = cmove(rep, sq2, v_is_negative)
    rep = cmove(rep, -rep, rep.is_negative()) # abs(rep)
    return rep


###########################################
# Fast Implementation (explicit formulas) #
###########################################
def dir_map_fast(r):
    """Computes a point (u, v) from the representative r in GF(p)

    Always succeeds
    """
    u  = r**2
    t1 = u * Z
    v  = t1 + GF(1)
    t2 = v**2
    t3 = A**2
    t3 = t3 * t1
    t3 = t3 - t2
    t3 = t3 * A
    t1 = t2 * v
    t1, is_square = inv_sqrt(t3 * t1)
    u  = u * ufactor  # no-op if ufactor == 1
    v  = r * vfactor  # copy  if vfactor == 1
    u  = cmove(u, GF(1), is_square)
    v  = cmove(v, GF(1), is_square)
    v  = v * t3
    v  = v * t1
    t1 = t1**2
    u  = u * -A
    u  = u * t3
    u  = u * t2
    u  = u * t1
    t1 = -v
    v  = cmove(v, t1, is_square != v.is_negative()) # use constant time XOR
    return (u, v)

def rev_map_fast(u, v_is_negative):
    """Computes the representative of the point (u, v), if possible

    Returns None if the point cannot be mapped.
    """
    t = u + A
    r = -Z * u
    r = r * t
    r, is_square = inv_sqrt(r)
    if not is_square:
        return None
    u = cmove(u, t, v_is_negative)
    r = u * r
    t = -r
    r = cmove(r, t, r.is_negative()) # abs(rep)
    return r


################################
# Compare both implementations #
################################
def dir_map(r):
    p_ref  = dir_map_ref(r)
    p_fast = dir_map_fast(r)
    p_neg  = dir_map_fast(-r)
    u, v   = p_ref
    r_back = rev_map_fast(u, v.is_negative())
    if p_ref  != p_fast : raise ValueError('ref/fast map mismatch')
    if p_ref  != p_neg  : raise ValueError('+r/-r map mismatch')
    if r_back != r.abs(): raise ValueError('roundtrip map mismatch')
    return p_ref

def rev_map(u, v_is_negative):
    r_ref  = rev_map_ref (u, v_is_negative)
    r_fast = rev_map_fast(u, v_is_negative)
    if r_ref != r_fast: raise ValueError('r mismatch (ref vs fast rev_map')
    if r_ref:
        u_back, v_back = dir_map_fast(r_ref)
        if u_back != u or v_back.is_negative() != v_is_negative:
            raise ValueError('roundtrip mismatch (rev_map)')
    return r_ref
