from Crypto.Cipher import AES
from Crypto.Util.Padding import pad, unpad
import hashlib
import secrets

def modinv(a, n):
    """
    Compute the modular inverse of a modulo n using the extended Euclidean Algorithm.
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
    See https://en.wikipedia.org/wiki/Jacobi_symbol
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
    return cubert

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
        """Compute a (multi) point multiplication
        ps is a list of (Jacobian tuple, scalar) pairs.
        """
        r = (0, 1, 0)
        for i in range(383, -1, -1):
            r = self.double(r)
            for (p, n) in ps:
                if ((n >> i) & 1):
                    r = self.add(r, p)
        return r
    
P_384_FIELD_SIZE = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000ffffffff
P_384_A = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000fffffffc
P_384_B = 0xb3312fa7e23ee7e4988e056be3f82d19181d9c6efe8141120314088f5013875ac656398d8a2ed19d2a85c8edd3ec2aef
P_384 = EllipticCurve(P_384_FIELD_SIZE, P_384_A, P_384_B)
P_384_G = (0xaa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7, 0x3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f, 1)
P_384_ORDER = 0xffffffffffffffffffffffffffffffffffffffffffffffffc7634d81f4372ddf581a0db248b0a77aecec196accc52973

CUBERT2 = fe(2).cubert()

def Solve_Quartic_Equation_Over_Fq(a, b, c, d, e):
    # require p mod 3 = 2 and p mod 4 = 3
    # a,b,c,d,e all field element
    result = []
    delta1 = c**2 - fe(3) * b * d + fe(12) * a * e
    delta2 = fe(2) * c**3 - fe(9) * b * c * d + fe(27) * a * d**2 + fe(27) * b**2 * e - fe(72) * a * c * e
    tmp1 = -fe(4) * delta1**3 + delta2**2
    if not tmp1.is_square():
        return []
    tmp2 = delta2 + tmp1.sqrt()
    tmp3 = tmp2.cubert()
    delta = CUBERT2 * delta1 / (fe(3) * a * tmp3) + tmp3 / (fe(3) * CUBERT2 * a)
    
    L = -b / (fe(4) * a)
    
    tmp4 = b ** 2 / (fe(4) * a**2) - fe(2) * c / (fe(3) * a) + delta
    if not tmp4.is_square() or tmp4 == fe(0):
        return []
    tmp5 = tmp4.sqrt()
    
    M = tmp5 / fe(2)
    
    R1 = b ** 2 / (fe(2) * a**2) - fe(4) * c / (fe(3) * a) - delta
    R2 = (-b ** 3 / (a**3) + fe(4) * b * c / (a**2) - fe(8) * d /a) / (fe(4) * tmp5)
    Tstx12 = R1 - R2
    Tstx34 = R1 + R2
    if Tstx12.is_square():
        R = (Tstx12.sqrt()) / fe(2)
        x1 = L - M - R
        x2 = L - M + R
        result.append(x1)
        result.append(x2)
    if Tstx34.is_square():
        R = (Tstx34.sqrt()) / fe(2)
        x3 = L + M - R
        x4 = L + M + R    
        result.append(x3)
        result.append(x4)
    
    return result

def forward_map(u):
    """Forward mapping function

    Parameters:
        u (of type fe) : any field element
    Returns:
        fe, fe : affine X and Y coordinates of a point on the P-384 curve
    """
    a = fe(P_384.a)
    b = fe(P_384.b)
    v = (fe(3)*a - u**4) / (fe(6) * u)
    t = (v**2 - b - u**6 / fe(27))
    x = t.cubert() + u**2 / fe(3)
    y = u * x + v
    return x, y

def reverse_map(x, y, i):
    """Reverse mapping function

    Parameters:
        fe, fe : X and Y coordinates of a point on the P-384 curve
        i      : integer in range [0,3]
    Returns:
        u (of type fe) : such that forward_map(u) = (x,y), or None.

        - There can be up to 4 such inverses, and i selects which formula to use.
        - Each i can independently from other i values return a value or None.
        - All non-None values returned across all 4 i values are guaranteed to be distinct.
        - Together they will cover all inverses of (x,y) under forward_map.
    """
    a = fe(1)
    b = fe(0)
    c = -fe(6) * x
    d = fe(6) * y
    e = -fe(3) * fe(P_384.a)
    lst = Solve_Quartic_Equation_Over_Fq(a, b, c, d, e)
    if len(lst) > i:
        return lst[i]
    return None

    
def encode(P, randombytes):
    # P -> u, v; forward_map(u)+forward_map(v) = P; 
    count = 0
    while True:
        # Random field element u and random number j is extracted from
        # SHA512("P-384_ellsq_encode\x00" + uint32{count} + rnd32 + X + byte{Y & 1}) 
        m = hashlib.sha512()
        m.update(b"P-384_ellsq_encode\x00")
        m.update(count.to_bytes(4, 'little'))
        m.update(randombytes)
        m.update(P[0].to_bytes(48, byteorder='big'))
        m.update((P[1] & 1).to_bytes(1, 'big'))
        hash512 = m.digest()
        u = fe(int.from_bytes(hash512[:48], 'big'))
        count += 1
        ge = forward_map(u)
        # convert ge to jacobian form for EC operations
        ge = (ge[0].val, ge[1].val, 1)
        T = P_384.negate(ge)
        Q = P_384.add(P_384.affine(P), P_384.affine(T))
        Q = P_384.affine(Q)
        if P_384.is_infinity(Q):
            Q = T
        j = secrets.randbelow(4)
        x, y, z = Q
        v = reverse_map(fe(x), fe(y), j)
        if v is not None:
            x1, y1 = forward_map(u)
            x2, y2 = forward_map(v)
            Sum = P_384.add((x1.val, y1.val, 1),(x2.val, y2.val, 1))
            Sum = P_384.affine(Sum)
            if (P[0] == Sum[0] and P[1] == Sum[1]):
                return u, v

def decode(u, v):
    # u, v -> P
    ge1 = forward_map(u)
    ge2 = forward_map(v)
    # convert ge1 and ge2 to jacobian form for EC operations
    T = ge1[0].val, ge1[1].val, 1
    S = ge2[0].val, ge2[1].val, 1
    P = P_384.add(T, S)
    if P_384.is_infinity(P):
        P = T
    P = P_384.affine(P)
    return fe(P[0]), fe(P[1])

class PRNEncryption_SECP256k1:
    '''
    IND$-CPA secure pseudorandom public-key encryption using admissible encoding.
    Launch on SECP256K1 using SW.
    '''
    def __init__(self):
        self.curve = P_384
        self.SK = secrets.randbelow(P_384_ORDER)
        self.PK = P_384.affine(P_384.mul([(P_384_G, self.SK)]))
        self.LP = P_384.p.bit_length()
        self.LT = 48
        self.bitlength = 384
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
        a = secrets.randbelow(P_384_ORDER)
        P = P_384.affine(P_384.mul([(P_384_G, a)]))
        x, y, z = P
        point_bytes = x.to_bytes(48, 'big') + y.to_bytes(48, 'big') + z.to_bytes(48, 'big')
        self.aes_key = hashlib.sha256(point_bytes).digest()
        
        # Point Hiding
        ge = (P[0], P[1], P[2])
        random_bytes = self.generate_random_bytes(32)
        u, v = encode(ge, random_bytes)
        
        # Bias Eliminating
        u_ = self.encode_bytes(u.val, P_384.p, self.LT, self.bitlength)
        v_ = self.encode_bytes(v.val, P_384.p, self.LT, self.bitlength)
        
        cipher = AES.new(self.aes_key, AES.MODE_CBC)
        self.nonce = cipher.iv
        ciphertext = cipher.encrypt(pad(plaintext.encode('utf-8'), AES.block_size))
        return u_.to_bytes(self.pointbytes,'big') + v_.to_bytes(self.pointbytes,'big') + self.nonce + ciphertext 

    def decrypt(self, ciphertext):
        u_ = int.from_bytes(ciphertext[:self.pointbytes],'big')
        v_ = int.from_bytes(ciphertext[self.pointbytes : 2 * self.pointbytes],'big')
        u = self.decode_bytes(u_, P_384.p, self.LT)
        v = self.decode_bytes(v_, P_384.p, self.LT)
        x, y = decode(fe(u), fe(v))
        z = 1
        point_bytes = x.val.to_bytes(48, 'big') + y.val.to_bytes(48, 'big') + z.to_bytes(48, 'big')
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
                    message = f"Message"  
                    ciphertext = encryption_system.encrypt(message)
                    prn_data = ciphertext[:min(len(ciphertext), num_bytes - bytes_written)]
                    f.write(prn_data)
                    bytes_written += len(prn_data)
                    pbar.update(len(prn_data))
        print(f"Generated {num_bits} bits of pseudorandom data and saved to {output_file}")

