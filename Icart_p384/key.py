import hashlib
import random
import cryptography
def modinv(a, n):
    """Compute the modular inverse of a modulo n using the extended Euclidean
    Algorithm. See https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Modular_integers.
    """
    # TODO: Change to pow(a, -1, n) available in Python 3.8
    t1, t2 = 0, 1
    r1, r2 = n, a
    while r2 != 0:
        q = r1 // r2
        t1, t2 = t2, t1 - q * t2
        r1, r2 = r2, r1 - q * r2
    if r1 > 1:
        return None
    if t1 < 0:
        t1 += n
    return t1

def jacobi_symbol(n, k):
    """Compute the Jacobi symbol of n modulo k

    See https://en.wikipedia.org/wiki/Jacobi_symbol

    For our application k is always prime, so this is the same as the Legendre symbol."""
    assert k > 0 and k & 1, "jacobi symbol is only defined for positive odd k"
    n %= k
    t = 0
    while n != 0:
        while n & 1 == 0:
            n >>= 1
            r = k & 7
            t ^= (r == 3 or r == 5)
        n, k = k, n
        t ^= (n & k & 3 == 3)
        n = n % k
    if k == 1:
        return -1 if t else 1
    return 0

def modsqrt(a, p):
    """
    Compute the square root of a modulo p when p % 4 = 3.
    """
    if p % 4 != 3:
        raise NotImplementedError("modsqrt only implemented for p % 4 = 3")
    sqrt = pow(a, (p + 1)//4, p)
    if pow(sqrt, 2, p) == a % p:
        return sqrt
    return None

def modcubert(a, p):
    """
    Compute the cube root of a modulo p when p % 3 = 2.
    """
    if p % 3 != 2:
        raise NotImplementedError("modcubert only implemented for p % 3 = 2")
    cubert = pow(a, (2*p-1)//3, p)
    if pow(cubert, 3, p) == a % p:
        return cubert
    return None

def TSmodsqrt(a, p):
    import math
    if (int(math.pow(a, (p - 1)/2)) % p != 1):
        return("No solutions")
    # find max power of 2 dividing p-1
    s = 0
    while((p - 1) % math.pow(2, s) == 0):
        s += 1
    s -= 1
    q = int((p - 1) / math.pow(2, s))# p-1=q*2^s
    # Select a z such that z is a quadratic non-residue modulo p
    z = 1
    res = int(math.pow(z, (p - 1) / 2)) % p
    while (res != p - 1):
        z += 1
        res = math.pow(z, (p - 1) / 2) % p
    c = int(math.pow(z, q)) % p
    r = int(math.pow(a, (q + 1) / 2)) % p
    t = int(math.pow(a, q)) % p
    m = s
    while(t % p != 1):
        i = 0
        div = False
        while (div == False):
            i += 1
            t = int(math.pow(t, 2)) % p
            if (t % p == 1):
                div = True
        b = int(math.pow(c, int(math.pow(2, m - i - 1)))) % p
        r = (r * b) % p
        t = t * (b ** 2) % p
        c = (b ** 2) % p
        m = i
        
    return r

P_384_FIELD_SIZE = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000ffffffff

class fe:
    """Prime field over F"""
    def __init__(self, x):
        if x is None:
            self.val = 0 
        else:
            self.val = x % P_384_FIELD_SIZE

    def __add__     (self, o): return fe(self.val + o.val)
    def __eq__      (self, o): return self.val == o.val
    def __hash__    (self   ): return id(self)
    def __mul__     (self, o): return fe(self.val * o.val)
    def __neg__     (self   ): return fe(-self.val)
    def __pow__     (self, s): return fe(pow(self.val, s, P_384_FIELD_SIZE))
    def __sub__     (self, o): return fe(self.val - o.val)
    def __truediv__ (self, o): return fe(self.val * o.invert().val)
    def __str__     (self): return str(self.val)

    def invert      (self   ):
        return fe(modinv(self.val, P_384_FIELD_SIZE))
    def is_odd(self): return (self.val & 1) != 0
    def is_square(self):
        return jacobi_symbol(self.val, P_384_FIELD_SIZE) >= 0
    def sqrt(self):
        return fe(modsqrt(self.val, P_384_FIELD_SIZE))
    def cubert(self):  # Only For p = 2 mod 3 
        return fe(modcubert(self.val, P_384_FIELD_SIZE))

    @staticmethod
    def from_bytes(b): return fe(int.from_bytes(b, 'big'))
    def to_bytes(self): return self.val.to_bytes(48, 'big')

class EllipticCurve:
    def __init__(self, p, a, b):
        """Initialize elliptic curve y^2 = x^3 + a*x + b over GF(p)."""
        self.p = p
        self.a = a % p
        self.b = b % p

    def affine(self, p1):
        """Convert a Jacobian point tuple p1 to affine form, or None if at infinity.

        An affine point is represented as the Jacobian (x, y, 1)"""
        x1, y1, z1 = p1
        if z1 == 0:
            return None
        inv = modinv(z1, self.p)
        inv_2 = (inv**2) % self.p
        inv_3 = (inv_2 * inv) % self.p
        return ((inv_2 * x1) % self.p, (inv_3 * y1) % self.p, 1)

    def has_even_y(self, p1):
        """Whether the point p1 has an even Y coordinate when expressed in affine coordinates."""
        return not (p1[2] == 0 or self.affine(p1)[1] & 1)

    def negate(self, p1):
        """Negate a Jacobian point tuple p1."""
        x1, y1, z1 = p1
        return (x1, (self.p - y1) % self.p, z1)

    def on_curve(self, p1):
        """Determine whether a Jacobian tuple p is on the curve (and not infinity)"""
        x1, y1, z1 = p1
        z2 = pow(z1, 2, self.p)
        z4 = pow(z2, 2, self.p)
        return z1 != 0 and (pow(x1, 3, self.p) + self.a * x1 * z4 + self.b * z2 * z4 - pow(y1, 2, self.p)) % self.p == 0

    def is_infinity(self, p1):
        """Return true if Jacobian tuple p is at infinity"""
        return p1[2] == 0

    def is_x_coord(self, x):
        """Test whether x is a valid X coordinate on the curve."""
        x_3 = pow(x, 3, self.p)
        return jacobi_symbol(x_3 + self.a * x + self.b, self.p) != -1

    def lift_x(self, x):
        """Given an X coordinate on the curve, return a corresponding affine point for which the Y coordinate is even."""
        x_3 = pow(x, 3, self.p)
        v = x_3 + self.a * x + self.b
        y = modsqrt(v, self.p)
        if y is None:
            return None
        return (x, self.p - y if y & 1 else y, 1)

    def double(self, p1):
        """Double a Jacobian tuple p1

        See https://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Jacobian_Coordinates - Point Doubling"""
        x1, y1, z1 = p1
        if z1 == 0:
            return (0, 1, 0)
        y1_2 = (y1**2) % self.p
        y1_4 = (y1_2**2) % self.p
        x1_2 = (x1**2) % self.p
        s = (4*x1*y1_2) % self.p
        m = 3*x1_2
        if self.a:
            m += self.a * pow(z1, 4, self.p)
        m = m % self.p
        x2 = (m**2 - 2*s) % self.p
        y2 = (m*(s - x2) - 8*y1_4) % self.p
        z2 = (2*y1*z1) % self.p
        return (x2, y2, z2)

    def add_mixed(self, p1, p2):
        """Add a Jacobian tuple p1 and an affine tuple p2

        See https://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Jacobian_Coordinates - Point Addition (with affine point)"""
        x1, y1, z1 = p1
        x2, y2, z2 = p2
        assert(z2 == 1)
        # Adding to the point at infinity is a no-op
        if z1 == 0:
            return p2
        z1_2 = (z1**2) % self.p
        z1_3 = (z1_2 * z1) % self.p
        u2 = (x2 * z1_2) % self.p
        s2 = (y2 * z1_3) % self.p
        if x1 == u2:
            if (y1 != s2):
                # p1 and p2 are inverses. Return the point at infinity.
                return (0, 1, 0)
            # p1 == p2. The formulas below fail when the two points are equal.
            return self.double(p1)
        h = u2 - x1
        r = s2 - y1
        h_2 = (h**2) % self.p
        h_3 = (h_2 * h) % self.p
        u1_h_2 = (x1 * h_2) % self.p
        x3 = (r**2 - h_3 - 2*u1_h_2) % self.p
        y3 = (r*(u1_h_2 - x3) - y1*h_3) % self.p
        z3 = (h*z1) % self.p
        return (x3, y3, z3)

    def add(self, p1, p2):
        """Add two Jacobian tuples p1 and p2

        See https://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Jacobian_Coordinates - Point Addition"""
        x1, y1, z1 = p1
        x2, y2, z2 = p2
        # Adding the point at infinity is a no-op
        if z1 == 0:
            return p2
        if z2 == 0:
            return p1
        # Adding an Affine to a Jacobian is more efficient since we save field multiplications and squarings when z = 1
        if z1 == 1:
            return self.add_mixed(p2, p1)
        if z2 == 1:
            return self.add_mixed(p1, p2)
        z1_2 = (z1**2) % self.p
        z1_3 = (z1_2 * z1) % self.p
        z2_2 = (z2**2) % self.p
        z2_3 = (z2_2 * z2) % self.p
        u1 = (x1 * z2_2) % self.p
        u2 = (x2 * z1_2) % self.p
        s1 = (y1 * z2_3) % self.p
        s2 = (y2 * z1_3) % self.p
        if u1 == u2:
            if (s1 != s2):
                # p1 and p2 are inverses. Return the point at infinity.
                return (0, 1, 0)
            # p1 == p2. The formulas below fail when the two points are equal.
            return self.double(p1)
        h = u2 - u1
        r = s2 - s1
        h_2 = (h**2) % self.p
        h_3 = (h_2 * h) % self.p
        u1_h_2 = (u1 * h_2) % self.p
        x3 = (r**2 - h_3 - 2*u1_h_2) % self.p
        y3 = (r*(u1_h_2 - x3) - s1*h_3) % self.p
        z3 = (h*z1*z2) % self.p
        return (x3, y3, z3)

    def mul(self, ps):
        """Compute a (multi) point multiplication

        ps is a list of (Jacobian tuple, scalar) pairs.
        """
        r = (0, 1, 0)
        for i in range(255, -1, -1):
            r = self.double(r)
            for (p, n) in ps:
                if ((n >> i) & 1):
                    r = self.add(r, p)
        return r

P_384_A = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000fffffffc
P_384_B = 0xb3312fa7e23ee7e4988e056be3f82d19181d9c6efe8141120314088f5013875ac656398d8a2ed19d2a85c8edd3ec2aef
P_384 = EllipticCurve(P_384_FIELD_SIZE, P_384_A, P_384_B)

P_384_G = (0xaa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7, 0x3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f, 1)
P_384_ORDER = 0xffffffffffffffffffffffffffffffffffffffffffffffffc7634d81f4372ddf581a0db248b0a77aecec196accc52973

class ECPubKey():
    """A secp256k1 public key"""

    def __init__(self):
        """Construct an uninitialized public key"""
        self.valid = False

    def set(self, data):
        """Construct a public key from a serialization in compressed or uncompressed format"""
        if (len(data) == 65 and data[0] == 0x04):
            p = (int.from_bytes(data[1:33], 'big'), int.from_bytes(data[33:65], 'big'), 1)
            self.valid = P_384.on_curve(p)
            if self.valid:
                self.p = p
                self.compressed = False
        elif (len(data) == 33 and (data[0] == 0x02 or data[0] == 0x03)):
            x = int.from_bytes(data[1:33], 'big')
            if P_384.is_x_coord(x):
                p = P_384.lift_x(x)
                # Make the Y coordinate odd if required (lift_x always produces
                # a point with an even Y coordinate).
                if data[0] & 1:
                    p = P_384.negate(p)
                self.p = p
                self.valid = True
                self.compressed = True
            else:
                self.valid = False
        else:
            self.valid = False

    def set_from_curve_point(self, curve_point):
        x, y = fe(curve_point[0]), fe(curve_point[1])
        if y.val % 2 == 0:
            compressed_sec = b'\x02' + x.val.to_bytes(48, 'big')
        else:
            compressed_sec = b'\x03' + x.val.to_bytes(48, 'big')
        self.set(compressed_sec)

    @property
    def is_compressed(self):
        return self.compressed

    @property
    def is_valid(self):
        return self.valid

    def get_bytes(self):
        assert(self.valid)
        p = P_384.affine(self.p)
        if p is None:
            return None
        if self.compressed:
            return bytes([0x02 + (p[1] & 1)]) + p[0].to_bytes(48, 'big')
        else:
            return bytes([0x04]) + p[0].to_bytes(48, 'big') + p[1].to_bytes(48, 'big')

    def get_group_element(self):
        assert(self.valid)
        p = P_384.affine(self.p)
        return fe(p[0]), fe(p[1])

    def verify_ecdsa(self, sig, msg, low_s=True):
        """Verify a strictly DER-encoded ECDSA signature against this pubkey.

        See https://en.wikipedia.org/wiki/Elliptic_Curve_Digital_Signature_Algorithm for the
        ECDSA verifier algorithm"""
        assert(self.valid)

        # Extract r and s from the DER formatted signature. Return false for
        # any DER encoding errors.
        if (sig[1] + 2 != len(sig)):
            return False
        if (len(sig) < 4):
            return False
        if (sig[0] != 0x30):
            return False
        if (sig[2] != 0x02):
            return False
        rlen = sig[3]
        if (len(sig) < 6 + rlen):
            return False
        if rlen < 1 or rlen > 33:
            return False
        if sig[4] >= 0x80:
            return False
        if (rlen > 1 and (sig[4] == 0) and not (sig[5] & 0x80)):
            return False
        r = int.from_bytes(sig[4:4+rlen], 'big')
        if (sig[4+rlen] != 0x02):
            return False
        slen = sig[5+rlen]
        if slen < 1 or slen > 33:
            return False
        if (len(sig) != 6 + rlen + slen):
            return False
        if sig[6+rlen] >= 0x80:
            return False
        if (slen > 1 and (sig[6+rlen] == 0) and not (sig[7+rlen] & 0x80)):
            return False
        s = int.from_bytes(sig[6+rlen:6+rlen+slen], 'big')

        # Verify that r and s are within the group order
        if r < 1 or s < 1 or r >= P_384_ORDER or s >= P_384_ORDER:
            return False
        if low_s and s >= P_384_ORDER // 2 :
            return False
        z = int.from_bytes(msg, 'big')

        # Run verifier algorithm on r, s
        w = modinv(s, P_384_ORDER)
        u1 = z*w % P_384_ORDER
        u2 = r*w % P_384_ORDER
        R = P_384.affine(P_384.mul([(P_384_G, u1), (self.p, u2)]))
        if R is None or (R[0] % P_384_ORDER) != r:
            return False
        return True

def generate_privkey():
    """Generate a valid random 48-byte private key."""
    return random.randrange(1, P_384_ORDER).to_bytes(48, 'big')

class ECKey():
    """A secp256k1 private key"""

    def __init__(self):
        self.valid = False

    def set(self, secret, compressed):
        """Construct a private key object with given 48-byte secret and compressed flag."""
        assert(len(secret) == 48)
        secret = int.from_bytes(secret, 'big')
        self.valid = (secret > 0 and secret < P_384_ORDER)
        if self.valid:
            self.secret = secret
            self.compressed = compressed

    def generate(self, compressed=True):
        """Generate a random private key (compressed or uncompressed)."""
        self.set(generate_privkey(), compressed)

    def get_bytes(self):
        """Retrieve the 48-byte representation of this key."""
        assert(self.valid)
        return self.secret.to_bytes(48, 'big')

    @property
    def is_valid(self):
        return self.valid

    @property
    def is_compressed(self):
        return self.compressed

    def get_pubkey(self):
        """Compute an ECPubKey object for this secret key."""
        assert(self.valid)
        ret = ECPubKey()
        p = P_384.mul([(P_384_G, self.secret)])
        ret.p = p
        ret.valid = True
        ret.compressed = self.compressed
        return ret

    def sign_ecdsa(self, msg, low_s=True, rfc6979=False):
        """Construct a DER-encoded ECDSA signature with this key.

        See https://en.wikipedia.org/wiki/Elliptic_Curve_Digital_Signature_Algorithm for the
        ECDSA signer algorithm."""
        assert(self.valid)
        z = int.from_bytes(msg, 'big')
        # Note: no RFC6979 by default, but a simple random nonce (some tests rely on distinct transactions for the same operation)
        if rfc6979:
            k = int.from_bytes(rfc6979_nonce(self.secret.to_bytes(48, 'big') + msg), 'big')
        else:
            k = random.randrange(1, P_384_ORDER)
        R = P_384.affine(P_384.mul([(P_384_G, k)]))
        r = R[0] % P_384_ORDER
        s = (modinv(k, P_384_ORDER) * (z + self.secret * r)) % P_384_ORDER
        if low_s and s > P_384_ORDER //2:
            s = P_384_ORDER - s
        # Represent in DER format. The byte representations of r and s have
        # length rounded up (255 bits becomes 48 bytes and 256 bits becomes 33
        # bytes).
        rb = r.to_bytes((r.bit_length() + 8) // 8, 'big')
        sb = s.to_bytes((s.bit_length() + 8) // 8, 'big')
        return b'\x30' + bytes([4 + len(rb) + len(sb), 2, len(rb)]) + rb + bytes([2, len(sb)]) + sb