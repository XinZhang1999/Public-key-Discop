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
```

### Run Single Example
The example is based on admissible encoding constructed by SW encoding on curve SECP256K1.

#### Elliptic Curve Equation

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



```shell
python src/run_single_example.py
```





### Change Admissible Encoding

If you need to deploy on different types of curves, we provide the core code for Icart, SW, and SWU, which are based on admissible encoding. These can be found in their respective named folders. Additionally, we offer Elligator2, which is based on point compression. You can modify the provided examples according to your needs. The code interfaces are quite similar and it won't take much time.

### Change Model

To change the model to llama2-7B, modify the code in `src/run_single_example.py` by setting the model_name as follows:

```python
settings.model_name = "llama2"
