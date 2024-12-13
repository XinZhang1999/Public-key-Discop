# Public-key Discop

Provably Secure Public-Key Steganography Based on Admissible Encoding

Xin Zhang, [Kejiang Chen](http://home.ustc.edu.cn/~chenkj/),  Na Zhao, [Weiming Zhang](http://staff.ustc.edu.cn/~zhangwm/), and [Nenghai Yu](http://staff.ustc.edu.cn/~ynh/)

SUBMITTED TO IEEE TRANSACTIONS ON INFORMATION FORENSICS AND SECURITY (TIFS)


## Usage

### Preparation

First, please ensure that you have installed all the required libraries for this repository.

We recommend using [Anaconda](https://anaconda.org/anaconda/conda) and execute the following commands.

```shell
conda create -n pkdiscop python=3.8.12
conda activate pkdiscop

# Visit the PyTorch website (https://pytorch.org/get-started/locally/) for installation commands tailored to your environment
# We have not tested PyTorch versions other than v1.12.0.
conda install pytorch==1.12.0 torchvision==0.13.0 torchaudio==0.12.0 cudatoolkit=11.3 -c pytorch

# Build the Cython files
python src/setup.py build_ext --build-lib=src/

# Make sure that you have downloaded models such as 'gpt2-medium', 'gpt2-large', or 'Llama-7b-hf'.
# You can modify the model path in the `src/model.py` file if needed to point to your local model files.
# These models can be downloaded from Hugging Face or other model repositories.

# Note: The `random_sample_cy.cpython-312-x86_64-linux-gnu.so` file we provided may not work directly in your environment. 
# If you encounter issues with it, please delete this file and rebuild it by running:
python src/setup.py build_ext --build-lib=src/
```

### Run Single Example

```shell
python src/run_single_example.py
```


The example is based on admissible encoding constructed by SW encoding on curve `SECP256K1`.

You can modify the `curve` parameter in the `run_single_example.py` file to select the elliptic curve. The available options for `curve` are:

- `p256`: The admissible encoding is constructed by SWU encoding.
- `p384`: The admissible encoding is constructed by SWU Icart.
- `secp256k1`: The admissible encoding is constructed by SW encoding.

By default, the curve is set to `secp256k1`, but you can change it to either `p256` or `p384` depending on your use case.

#### Elliptic Curve Equation of SECP256k1

$$
y^2 \equiv x^3 + ax + b
$$

#### Parameters

| Name | Value |
| --- | --- |
| p | 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f |
| a | 0x0000000000000000000000000000000000000000000000000000000000000000 |
| b | 0x0000000000000000000000000000000000000000000000000000000000000007 |
| G | (0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798, 0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8) |
| n | 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141 |
| h | 0x01 |

#### Elliptic Curve Equation of P256

$$
y^2 \equiv x^3 + ax + b
$$

#### Parameters

| Name | Value |
| --- | --- |
| p | 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff |
| a | 0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc |
| b | 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b |
| G | (0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296, 0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5) |
| n | 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551 |
| h | 0x01 |

#### Elliptic Curve Equation of P384

$$
y^2 \equiv x^3 + ax + b
$$

#### Parameters

| Name | Value |
| --- | --- |
| p | 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000ffffffff |
| a | 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000fffffffc |
| b | 0xb3312fa7e23ee7e4988e056be3f82d19181d9c6efe8141120314088f5013875ac656398d8a2ed19d2a85c8edd3ec2aef |
| G | (0xaa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7, 0x3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f) |
| n | 0xffffffffffffffffffffffffffffffffffffffffffffffffc7634d81f4372ddf581a0db248b0a77aecec196accc52973 |
| h | 0x01 |



### IND$-CPA Encryption and Pseudorandomness

You can directly run the following commands to generate the corresponding 100 million bits of binary encryption results:

- `python src/PRNEncryption_P256.py -generate`
- `python src/PRNEncryption_P384.py -generate`
- `python src/PRNEncryption_SECP256k1.py -generate`

These outputs can be tested for pseudorandomness using the NIST Statistical Test Suite, a widely used standard for evaluating the randomness of binary sequences. The test suite provides a comprehensive set of tests to determine whether the generated sequences exhibit statistical properties expected of truly random sequences.





## Steganalysis 

The complete code will be released later.
