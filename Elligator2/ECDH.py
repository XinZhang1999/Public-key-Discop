from curve25519 import *
from elligator import *
from random    import randrange
from random    import seed
import secrets
import time

def random_point():
    x = GF(randrange(0, Mt.order))
    xP = Mt.scalarmult(Mt.base, x.val)
    r = rev_map(xP, False)
    while r is None:
        x = GF(randrange(0, Mt.order))
        xP = Mt.scalarmult(Mt.base, x.val)
        r = rev_map(xP, False)
    return x, xP
    

def int_254_tobytes_256(inp):
    # add 2 random bits on last bytes
    # inp: 256 bits int;  return 32bytes;
    # 生成一个2位的随机数
    random_2bits = secrets.randbelow(4)  # 4表示2的2次方，即2位的可能值数量
    bytelist = inp.to_bytes(32, byteorder='little')
    handle = bytelist[31]
    result_8bits = (random_2bits << 6) | handle
    return bytelist[0:31]+result_8bits.to_bytes(1, byteorder='little')


def bytes_256_to_int_254(byteslist):
    handle = byteslist[31]
    handle %=64
    result_254bytes = byteslist[0:31]+handle.to_bytes(1, byteorder='little')
    return int.from_bytes(result_254bytes, byteorder = 'little')


x_a, xP_a = random_point()
r_a = rev_map(xP_a, False)

x_b, xP_b = random_point()
r_b = rev_map(xP_b, False)

receive_a = int_254_tobytes_256(r_b.val)
receive_b = int_254_tobytes_256(r_a.val)

PK_a = dir_map_ref(GF(bytes_256_to_int_254(receive_b)))
PK_b = dir_map_ref(GF(bytes_256_to_int_254(receive_a)))


K1 = Mt.scalarmult(PK_a[0], x_b.val)
K2 = Mt.scalarmult(PK_b[0], x_a.val)

print(K1)
print(K2)