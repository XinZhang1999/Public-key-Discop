import sys

# collect arguments
if len(sys.argv) < 3:
    raise ValueError('Usage: gen_vectors.py curve vectors')
curve   = sys.argv[1]
vectors = sys.argv[2]

# Import curve module
if   curve == "curve25519": from curve25519 import *
elif curve == "curve448"  : from curve448   import *
else: raise ValueError('Uknnown curve module')

# remaining imports
from elligator import *
from random    import randrange
from random    import seed


##############
# Direct map #
##############
def direct_map_all_vectors():
    """All test vectors for the direct map"""
    seed(12345)  # cheap determinism for the random test vectors
    vectors = []

    # Representative 0 maps to the point 0, 0
    vectors.append(vectors_to_string([0, 0, 0]))

    # Random representatives map to their respective point
    for _ in range(256):
        r    = GF(randrange(0, GF.p - 1)).abs()
        u, v = dir_map(r)
        vectors.append(vectors_to_string([r, u, v]))

    return "\n\n".join(vectors)


###############
# Reverse map #
###############
def random_curve_point():
    u = GF(randrange(0, GF.p - 1))
    while not is_square(u**3 + A * u**2 + B * u):
        u = GF(randrange(0, GF.p - 1))
    return u

def reverse_map_all_vectors():
    """All test vectors for the reverse map"""
    seed(12345)  # cheap determinism for the random test vectors
    vectors = []

    # point (0, 0) maps to representative 0
    vectors.append(vectors_to_string([0, False, "00:", 0]))

    # some points that do not map
    for i in range(16):
        u = random_curve_point()
        r = rev_map(u, False)
        while not r is None:
            u = random_curve_point()
            r = rev_map(u, False)
        if not rev_map(u, True) is None:
            raise ValueError('Reverse map should fail')
        vectors.append(vectors_to_string([u, False, "ff:", ":"]))
        vectors.append(vectors_to_string([u, True , "ff:", ":"]))

    # lots of points that do map
    for i in range(256):
        u  = random_curve_point()
        rp = rev_map(u, False)
        while rp is None:
            u  = random_curve_point()
            rp = rev_map(u, False)
        rn = rev_map(u, True)
        if rn is None: raise ValueError('Reverse map should succeed')
        vectors.append(vectors_to_string([u, False, "00:", rp]))
        vectors.append(vectors_to_string([u, True , "00:", rn]))

    return "\n\n".join(vectors)


##############
# Scalarmult #
##############
def scalarmult_all_vectors():
    """All test vectors for scalar multiplication"""
    seed(12345)  # cheap determinism for the random test vectors
    vectors = []
    for i in range(64):
        c      = i % Mt.cofactor
        scalar = randrange(2**(GF.nb_bytes * 8))      # lower bits = random
        scalar = scalar // Mt.cofactor * Mt.cofactor  # lower bits = 0
        scalar = scalar + c                           # lower bits = c
        vectors.append(vectors_to_string([
            scalar,
            co_scalarmult(scalar, c)
        ]))
    return "\n\n".join(vectors)


################
# Main program #
################
vectors_map = {"direct"    : direct_map_all_vectors,
               "inverse"   : reverse_map_all_vectors,
               "scalarmult": scalarmult_all_vectors,
               }
print(vectors_map[vectors]())