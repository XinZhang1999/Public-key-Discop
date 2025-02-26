from Crypto.Cipher import AES
from Crypto.Util.Padding import pad, unpad
import hashlib
import secrets

def modinv(a, n):
    """
    Compute the modular inverse of a modulo n using the extended Euclidean Algorithm
    """
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
    """
    Compute the Jacobi symbol of n modulo k
    For our application k is always prime, so this is the same as the Legendre symbol.
    """
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
    The Tonelli-Shanks algorithm can be used. See https://en.wikipedia.org/wiki/Tonelli-Shanks_algorithm
    Limiting this function to only work for p % 4 = 3 means we don't need to
    iterate through the loop. The highest n such that p - 1 = 2^n Q with Q odd
    is n = 1. Therefore Q = (p-1)/2 and sqrt = a^((Q+1)/2) = a^((p+1)/4)
    secp256k1's is defined over field of size 2**256 - 2**32 - 977, which is 3 mod 4.
    """
    if p % 4 != 3:
        raise NotImplementedError("modsqrt only implemented for p % 4 = 3")
    sqrt = pow(a, (p + 1)//4, p)
    if pow(sqrt, 2, p) == a % p:
        return sqrt
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

class fe:
    """
    Prime field over 2^256 - 2^32 - 977
    """
    def __init__(self, x):
        if x is None:
            self.val = 0 
        else:
            self.val = x % SECP256K1_FIELD_SIZE

    def __add__     (self, o): return fe(self.val + o.val)
    def __eq__      (self, o): return self.val == o.val
    def __hash__    (self   ): return id(self)
    def __mul__     (self, o): return fe(self.val * o.val)
    def __neg__     (self   ): return fe(-self.val)
    def __pow__     (self, s): return fe(pow(self.val, s, SECP256K1_FIELD_SIZE))
    def __sub__     (self, o): return fe(self.val - o.val)
    def __truediv__ (self, o): return fe(self.val * o.invert().val)
    def __str__     (self): return str(self.val)

    def invert      (self   ):
        return fe(modinv(self.val, SECP256K1_FIELD_SIZE))
    def is_odd(self): return (self.val & 1) != 0
    def is_square(self):
        return jacobi_symbol(self.val, SECP256K1_FIELD_SIZE) >= 0
    def sqrt(self):
        return fe(modsqrt(self.val, SECP256K1_FIELD_SIZE))

    @staticmethod
    def from_bytes(b): return fe(int.from_bytes(b, 'big'))
    def to_bytes(self): return self.val.to_bytes(32, 'big')
    
class EllipticCurve:
    def __init__(self, p, a, b):
        """
        Initialize elliptic curve y^2 = x^3 + a*x + b over GF(p).
        """
        self.p = p
        self.a = a % p
        self.b = b % p

    def affine(self, p1):
        """
        Convert a Jacobian point tuple p1 to affine form, or None if at infinity.
        An affine point is represented as the Jacobian (x, y, 1)
        """
        x1, y1, z1 = p1
        if z1 == 0:
            return None
        inv = modinv(z1, self.p)
        inv_2 = (inv**2) % self.p
        inv_3 = (inv_2 * inv) % self.p
        return ((inv_2 * x1) % self.p, (inv_3 * y1) % self.p, 1)

    def has_even_y(self, p1):
        """
        Whether the point p1 has an even Y coordinate when expressed in affine coordinates.
        """
        return not (p1[2] == 0 or self.affine(p1)[1] & 1)

    def negate(self, p1):
        """
        Negate a Jacobian point tuple p1.
        """
        x1, y1, z1 = p1
        return (x1, (self.p - y1) % self.p, z1)

    def on_curve(self, p1):
        """
        Determine whether a Jacobian tuple p is on the curve (and not infinity)
        """
        x1, y1, z1 = p1
        z2 = pow(z1, 2, self.p)
        z4 = pow(z2, 2, self.p)
        return z1 != 0 and (pow(x1, 3, self.p) + self.a * x1 * z4 + self.b * z2 * z4 - pow(y1, 2, self.p)) % self.p == 0

    def is_infinity(self, p1):
        """
        Return true if Jacobian tuple p is at infinity
        """
        return p1[2] == 0

    def is_x_coord(self, x):
        """
        Test whether x is a valid X coordinate on the curve.
        """
        x_3 = pow(x, 3, self.p)
        return jacobi_symbol(x_3 + self.a * x + self.b, self.p) != -1

    def lift_x(self, x):
        """
        Given an X coordinate on the curve, return a corresponding affine point for which the Y coordinate is even.
        """
        x_3 = pow(x, 3, self.p)
        v = x_3 + self.a * x + self.b
        y = modsqrt(v, self.p)
        if y is None:
            return None
        return (x, self.p - y if y & 1 else y, 1)

    def double(self, p1):
        """
        Double a Jacobian tuple p1
        See https://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Jacobian_Coordinates - Point Doubling
        """
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
        """
        Add a Jacobian tuple p1 and an affine tuple p2
        See https://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Jacobian_Coordinates - Point Addition (with affine point)
        """
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
        """
        Add two Jacobian tuples p1 and p2
        See https://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Jacobian_Coordinates - Point Addition
        """
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
        """
        Compute a (multi) point multiplication
        ps is a list of (Jacobian tuple, scalar) pairs.
        """
        r = (0, 1, 0)
        for i in range(255, -1, -1):
            r = self.double(r)
            for (p, n) in ps:
                if ((n >> i) & 1):
                    r = self.add(r, p)
        return r

SECP256K1_FIELD_SIZE = 2**256 - 2**32 - 977
SECP256K1 = EllipticCurve(SECP256K1_FIELD_SIZE, 0, 7)
SECP256K1_G = (0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798, 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8, 1)
SECP256K1_ORDER = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

C1 = fe(-3).sqrt()
C2 = (C1 - fe(1)) / fe(2)
B = fe(7)

def forward_map(u):
    """Forward mapping function

    Parameters:
        u (of type fe) : any field element
    Returns:
        fe, fe : affine X and Y coordinates of a point on the secp256k1 curve
    """
    s = u**2
    x1 = C2 - C1*s / (fe(1)+B+s)
    g1 = x1**3 + B
    if g1.is_square():
        x, g = x1, g1
    else:
        x2 = -x1 - fe(1)
        g2 = x2**3 + B
        if g2.is_square():
            x, g = x2, g2
        else:
            x3 = fe(1) - (fe(1)+B+s)**2 / (fe(3)*s)
            g3 = x3**3 + B
            x, g = x3, g3
    y = g.sqrt()
    if y.is_odd() == u.is_odd():
        return x, y
    else:
        return x, -y

def reverse_map(x, y, i):
    """Reverse mapping function

    Parameters:
        fe, fe : X and Y coordinates of a point on the secp256k1 curve
        i      : integer in range [0,3]
    Returns:
        u (of type fe) : such that forward_map(u) = (x,y), or None.

        - There can be up to 4 such inverses, and i selects which formula to use.
        - Each i can independently from other i values return a value or None.
        - All non-None values returned across all 4 i values are guaranteed to be distinct.
        - Together they will cover all inverses of (x,y) under forward_map.
    """
    if i == 0 or i == 1:
        z = fe(2)*x + fe(1)
        t1 = C1 - z
        t2 = C1 + z
        if not (t1*t2).is_square():
            return None
        if i == 0:
            if t2 == fe(0):
                return None
            if t1 == fe(0) and y.is_odd():
                return None
            u = ((fe(1)+B)*t1/t2).sqrt()
        else:
            x1 = -x-fe(1)
            if (x1**3 + B).is_square():
                return None
            u = ((fe(1)+B)*t2/t1).sqrt()
    else:
        z = fe(2) - fe(4)*B - fe(6)*x
        if not (z**2 - fe(16)*(B+fe(1))**2).is_square():
            return None
        if i == 2:
            s = (z + (z**2 - fe(16)*(B+fe(1))**2).sqrt()) / fe(4)
        else:
            if z**2 == fe(16)*(B+fe(1))**2:
                return None
            s = (z - (z**2 - fe(16)*(B+fe(1))**2).sqrt()) / fe(4)
        if not s.is_square():
            return None
        x1 = C2 - C1*s / (fe(1)+B+s)
        if (x1**3 + B).is_square():
            return None
        u = s.sqrt()
    if y.is_odd() == u.is_odd():
        return u
    else:
        return 
    
def encode(P, randombytes):
    '''
    - P -> u, v:
    - P is a point on an elliptic curve.
    - u and v are two field elements derived from P.
    - The relationship between P, u, and v satisfies the admissible encoding condition:
      forward_map(u) + forward_map(v) = P.
    '''
    count = 0
    while True:
        '''
        Random field element u 
        SHA256("secp256k1_ellsq_encode\x00" + uint32{count} + rnd32 + X + byte{Y & 1})
        '''
        m = hashlib.sha256()
        m.update(b"secp256k1_ellsq_encode\x00")
        m.update(count.to_bytes(4, 'little'))
        m.update(randombytes)
        m.update(P[0].to_bytes(32, byteorder='big'))
        m.update((P[1] & 1).to_bytes(1, 'big'))
        hash = m.digest()
        u = fe(int.from_bytes(hash, 'big'))
        count += 1
        ge = forward_map(u)
        # convert ge to jacobian form for EC operations
        ge = (ge[0].val, ge[1].val, 1)
        T = SECP256K1.negate(ge)
        Q = SECP256K1.add(T, SECP256K1.affine(P))
        if SECP256K1.is_infinity(Q):
            Q = T
        j = secrets.randbelow(4)
        Q = SECP256K1.affine(Q)
        v = reverse_map(fe(Q[0]), fe(Q[1]), j)
        if v is not None:
            return u, v

def decode(u, v):
    '''
    The process u, v -> P is the inverse of forward_map: 
    It reconstructs the elliptic curve point P from two field elements u and v.
    '''
    ge1 = forward_map(u)
    ge2 = forward_map(v)
    # convert ge1 and ge2 to jacobian form for EC operations
    T = ge1[0].val, ge1[1].val, 1
    S = ge2[0].val, ge2[1].val, 1
    P = SECP256K1.add(T, S)
    if SECP256K1.is_infinity(P):
        P = T
    P = SECP256K1.affine(P)
    return fe(P[0]), fe(P[1])

class PRNEncryption_SECP256k1:
    '''
    IND$-CPA secure pseudorandom public-key encryption using admissible encoding.
    Launch on SECP256K1 using SW.
    '''
    def __init__(self):
        self.curve = SECP256K1
        self.SK = secrets.randbelow(SECP256K1_ORDER)
        self.PK = SECP256K1.affine(SECP256K1.mul([(SECP256K1_G, self.SK)]))
        self.LP = SECP256K1.p.bit_length()
        self.LT = 32
        self.bitlength = 256
        self.pointbytes = (self.bitlength + self.LT) // 8
        
    def set_redundancy(self, t):
        self.LT = t
    
    def generate_random_bytes(self, num_bytes):
        return secrets.token_bytes(num_bytes)

    def encode_bytes(self, x, P, LT, K):
        k = secrets.randbelow(((2 ** K) * (2 ** LT) - x) // P)
        return x + k * P

    def decode_bytes(self, enc, P, LT):
        return enc % P
    
    def encrypt(self, plaintext):
        ## Key Deriving
        a = secrets.randbelow(SECP256K1_ORDER)
        P = SECP256K1.affine(SECP256K1.mul([(SECP256K1_G, a)]))
        x, y, z = P
        point_bytes = x.to_bytes(32, 'big') + y.to_bytes(32, 'big') + z.to_bytes(32, 'big')
        self.aes_key = hashlib.sha256(point_bytes).digest()
        
        # Point Hiding
        ge = (P[0], P[1], P[2])
        random_bytes = self.generate_random_bytes(32)
        u, v = encode(ge, random_bytes)
        
        # Bias Eliminating
        u_ = self.encode_bytes(u.val, SECP256K1.p, self.LT, self.bitlength)
        v_ = self.encode_bytes(v.val, SECP256K1.p, self.LT, self.bitlength)
        
        cipher = AES.new(self.aes_key, AES.MODE_CBC)
        self.nonce = cipher.iv
        ciphertext = cipher.encrypt(pad(plaintext.encode('utf-8'), AES.block_size))
        return u_.to_bytes(self.pointbytes,'big') + v_.to_bytes(self.pointbytes,'big') + self.nonce + ciphertext 

    def decrypt(self, ciphertext):
        u_ = int.from_bytes(ciphertext[:self.pointbytes],'big')
        v_ = int.from_bytes(ciphertext[self.pointbytes : 2 * self.pointbytes],'big')
        u = self.decode_bytes(u_, SECP256K1.p, self.LT)
        v = self.decode_bytes(v_, SECP256K1.p, self.LT)
        x, y = decode(fe(u), fe(v))
        z = 1
        point_bytes = x.val.to_bytes(32, 'big') + y.val.to_bytes(32, 'big') + z.to_bytes(32, 'big')
        self.aes_key = hashlib.sha256(point_bytes).digest()
        
        nonce = ciphertext[self.pointbytes * 2 : self.pointbytes * 2 + AES.block_size]
        cipher = AES.new(self.aes_key, AES.MODE_CBC, nonce)
        plaintext = unpad(cipher.decrypt(ciphertext[self.pointbytes * 2 + AES.block_size:]), AES.block_size)
        return plaintext.decode('utf-8')
    
def save_bytes_as_binary_text(ciphertext, filename='message.txt'):

    binary_string = ''.join(format(byte, '08b') for byte in ciphertext)
    with open(filename, 'w') as file:
        for bit in binary_string:
            file.write('0' if bit == '0' else '1')

import argparse
from tqdm import tqdm


generate_bit = False
parser = argparse.ArgumentParser()
parser.add_argument('-generate_bit', action='store_true')
args = parser.parse_args()
if args.generate_bit:
    generate_bit = True
else:
    generate_bit = False

if __name__ == '__main__':
    if not generate_bit:
        encryption_system = PRNEncryption_SECP256k1()
        message = "Attack at 9:00"
        ciphertext = encryption_system.encrypt(message)
        print(f"Ciphertext: {ciphertext.hex()}")
        decrypted_message = encryption_system.decrypt(ciphertext)
        print(f"Decrypted Message: {decrypted_message}")
        assert message == decrypted_message, "Decryption failed!"
        print("Encryption and decryption succeeded.")
    else:
        encryption_system = PRNEncryption_SECP256k1()
        output_file = "prn_test_data"
        num_bits = 100_000_000
        num_bytes = num_bits // 8
        with tqdm(total=num_bytes, unit='B', unit_scale=True, desc="Writing PRN Data") as pbar:
            with open(output_file, "wb") as f:
                bytes_written = 0
                while bytes_written < num_bytes:
                    message = f"Message {bytes_written}"  
                    ciphertext = encryption_system.encrypt(message)
                    prn_data = ciphertext[:min(len(ciphertext), num_bytes - bytes_written)]
                    f.write(prn_data)
                    bytes_written += len(prn_data)
                    pbar.update(len(prn_data))
        print(f"Generated {num_bits} bits of pseudorandom data and saved to {output_file}")

