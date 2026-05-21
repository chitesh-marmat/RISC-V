#  RISC-V Matrix-Vector Implementation

**Course:** ELL305, IIT Delhi  
**ISA:** RV32IMF (RARS simulator)

---

## Problem Statement

Given a 6×6 floating-point matrix **A** and a 6×1 vector **X**:

1. Compute **Y = A × X**
2. Compute **det(A)** via LU decomposition with partial pivoting
3. Print the boolean flag: `sum(Y) < det(A)` (0 or 1)

---

## Running

Requires [RARS](https://github.com/TheThirdOne/rars) (RISC-V Assembler and Runtime Simulator).

```bash
java -jar rars.jar final.asm
```

Or via the RARS GUI: File → Open `final.asm` → Assemble (F3) → Run (F5).

**Expected output (three lines):**
```
<6 space-separated floats>   ← Y vector
<1 float>                    ← det(A)
<0 or 1>                     ← comparison flag
```

---

## Implementation

The program is split into three subroutines, all following the RISC-V calling convention (callee-saved `s0`/`s1` saved on stack).

#### `matvec_mul` — Matrix-vector multiply
Computes `Y[row] = Σ A[row][col] * X[col]` with a nested row/column loop.  
Row-major addressing: `&A[i][j] = base + (i*N + j)*4`.

#### `det_approx` — Determinant via LU decomposition
Performs Gaussian elimination with partial pivoting **in-place** on matrix `A`.  
Per pivot column `i`:
- Find the row with the largest `|A[k][i]|` (abs via `fmv.x.w` + bitmask `0x7FFFFFFF`)
- Swap rows and flip the sign accumulator
- Eliminate all rows below the pivot

The running product of all pivots (with sign) gives `det(A)`.

> **Note:** `det_approx` overwrites `A` — `matvec_mul` must be called first.

#### `compare_sum_det` — Comparison
Sums `Y[0..N-1]` and compares against the stored `det_result` using `flt.s`.

---

## Key Notes

| Detail | Explanation |
|--------|-------------|
| Float abs | No `fabs.s` in RARS; done with `fmv.x.w` → `and 0x7FFFFFFF` → `fmv.w.x` |
| Sign tracking | Row swaps multiply a sign factor by −1, applied to `ft9` via `fmul.s ft9, ft9, ft8` |
| Stack frames | `det_approx`: 32 bytes (4 saved regs); others: 8 bytes (2 saved regs) |
| Rounding mode | Some drafts use `, dyn` suffix — drop it for stricter assemblers |

---

## File Structure

```
.
├── final.asm          # submitted solution
├── docs/
│   ├── Assignment_RV32IMF_Matrix_Determinant_6x6.pdf   # problem statement
│   └── report.pdf                                       # submission report
└── drafts/            # earlier iterations
    ├── round1.asm
    ├── change1.asm
    ├── change2.asm
    ├── riscv3.asm
    └── working
```
