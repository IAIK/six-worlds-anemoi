# Exploring the Six Worlds of Gröbner Basis Cryptanalysis: Application to Anemoi

This repository accompanies the paper **Exploring the Six Worlds of Gröbner Basis Cryptanalysis: Application to Anemoi**. It contains the code used to generate the models, create Magma scripts, and evaluate results.


## Files

Overview and evaluation of results: [`results.ipynb`](results.ipynb)

Other important files:

- [`anemoi.sage`](anemoi.sage): Implementation of Anemoi, following the one provided by the authors [here](https://github.com/anemoi-hash/anemoi-hash).
- [`attack.sage`](attack.sage): Attack script calling Magma functions.
    - Parsing log output from attack: [`parseResults.sage`](parseResults.sage) (saves to `results_Fp.sobj` and `results_F2n.sobj`)
- [`constants.py`](constants.py): Relevant constants, such as large named primes.
- [`models.sage`](models.sage): Implementation of the algebraic models $F_{CICO}$ and $P_{CICO}$ and 3 variable orderings for every model.
    - File for testing model implementation: [`models_poc.ipynb`](models_poc.ipynb)
- [`SystemAnalysis.sage`](SystemAnalysis.sage): Bézout bound, multihomogeneous Bézout bound, etc.
- [`modelBounds.sage`](modelBounds.sage): All derived formulas for theoretical bounds and experimental conjectures.


## Possible reduced instances for Anemoi for different field sizes

For specifications and constraints, see [Anemoi paper](https://eprint.iacr.org/2022/840.pdf), page 10.


### Prime fields $\mathbb{F}_p$

| Prime                                     | $\lceil \log_2(p) \rceil $ | $3$ | $5$ | $7$ | $9$ | $11$ |
|-------------------------------------------|----------------------------|-----|-----|-----|-----|------|
| 65537 = 0x10001                           | $16$                       | yes | yes | yes | yes | yes  |
| 4294967087 = 0xffffff2f                   | $32$                       | yes | yes | yes | yes | yes  |
| 18446744073709551263 = 0xfffffffffffffe9f | $64$                       | yes | yes | yes | yes | yes  |
| PALLAS_BASEFIELD                          | $255$                      | no  | yes | yes | no  | yes  |
| VESTA_BASEFIELD                           | $255$                      | no  | yes | yes | no  | yes  |
| BLS12_377_SCALARFIELD                     | $253$                      | no  | no  | no  | no  | yes  |
| BLS12_381_SCALARFIELD                     | $255$                      | no  | yes | yes | no  | no   |
| BN_254_SCALARFIELD                        | $254$                      | no  | yes | yes | no  | yes  |



### Binary fields $\mathbb{F}_{2^n}$

| $n$   | 3    |  5   |   9  |
|-------|------|------|------|
| $15$  | yes  | yes  | no   |
| $17$  | yes  | yes  | yes  |
| $31$  | yes  | yes  | yes  |
| $63$  | yes  | yes  | no   |
| $65$  | yes  | yes  | yes  |
| $127$ | yes  | yes  | yes  |
| $255$ | yes  | yes  | no   |
| $257$ | yes  | yes  | yes  |

Note: $\alpha = 7$ and $\alpha=11$ not possible due to construction of $\alpha$.
