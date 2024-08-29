from curve25519 import *
from elligator import *
from random    import randrange
from random    import seed
import hashlib

def map_to_curve(random):
    """Maps a uniform random 256 bit number to a curve point

    This is compatible with libsodium.

    Contrary to what the Hash to Curve RFC recommends, we don't use a
    big hash to reduce it modulo p.  Instead we just chop off one
    bit. The resulting field is close enough to a power of two that the
    deviation from perfect randomness is undetectable.
    """
    y_sign  = random // 2**255            # Get sign of Edwards x coordinate
    r       = random %  2**255            # Elligator representative
    u, _    = dir_map(GF(r))              # Ignore Montgomery v coordinate
    x, y, z = to_edwards(u)               # Convert to Edwards
    if x.to_num() % 2 != y_sign: x = -x   # Set sign of Edwards x coordinate
    x, y, z = Ed.scalarmult((x, y, z), 8) # Multiply by cofactor

    # Serialise Edwards point (divide, get sign of x)
    z       = z.invert()
    y       = y * z
    x       = x * z
    x_sign  = x.to_num() % 2              # Negative means odd here.
    point   = y.to_num() + x_sign * 2**255
    return point

# Generate the actual test vectors, print them in stdout.
seed(12345)  # cheap determinism for the random test vectors
vectors = []
for i in range(64):
    r = randrange(2**256)
    p = map_to_curve(r)
    vectors.append(vectors_to_string([r, p]))
print("\n\n".join(vectors))