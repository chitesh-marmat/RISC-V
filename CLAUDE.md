# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

ELL305 assignment at IIT Delhi — RISC-V RV32IMF assembly program that:
1. Computes `Y = A * X` (6×6 matrix times 6×1 vector)
2. Computes `det(A)` via LU decomposition with partial pivoting
3. Prints `sum(Y) < det(A)` as a boolean flag (0 or 1)

The submitted file is `final.asm`. Other `.asm` files (`round1.asm`, `change1.asm`, `change2.asm`, `riscv3.asm`, `working`) are earlier drafts and alternatives.

## Running

Load and run `.asm` files in **RARS** (RISC-V Assembler and Runtime Simulator). RARS is a JAR-based GUI tool:

```
java -jar rars.jar final.asm
```

Or open RARS GUI → File → Open → Assemble (F3) → Run (F5).

Expected output (three lines):
1. Six space-separated floats — the Y vector
2. One float — the determinant
3. One integer (0 or 1) — the comparison flag

## Architecture

All three subroutines follow the RISC-V calling convention with explicit stack save/restore of `s0`/`s1` (and `s2`/`s3` in `det_approx`). Callee-saved registers hold base pointers; temporaries `t0`–`t6` hold loop indices and scratch values.

### `matvec_mul` (a0=A, a1=X, a2=Y, a3=N)
Nested loop over rows then columns. Row-major address: `A[row][col] = base + (row*N + col)*4`. Accumulates dot product in `ft3`, stores to `Y[row]`.

### `det_approx`
Modifies matrix `A` in-place (address hardcoded via `la s0, A`). Three nested phases per pivot column `i`:
- **Pivot search** — find row with max `|A[k][i]|` below row `i` using bit-mask `0x7FFFFFFF` for float abs
- **Row swap** — swap rows `i` and `max_row`; multiply `ft9` by −1 to track sign
- **Elimination** — for each row `k > i`, subtract `(A[k][i]/A[i][i]) * row_i` from `row_k`

Running product `ft9` (initialized to 1.0) accumulates all pivots. Stored to `det_result`.

### `compare_sum_det`
Sums `Y[0..N-1]` into `fa0`, loads `det_result` into `ft1`, uses `flt.s a0, fa0, ft1` to set `comparison_flag`.

## Key implementation notes

- **Float abs via bit manipulation**: RISC-V RV32F has no `fabs.s` in all simulators; abs is done by moving to integer register with `fmv.x.w`, masking with `0x7FFFFFFF`, then back with `fmv.w.x`.
- **`det_approx` modifies A**: The matrix stored at label `A` is overwritten during LU factorization. `matvec_mul` must be called first.
- **Rounding mode suffix**: Some drafts use `fmul.s ft6, ft5, ft7, dyn` — the `, dyn` rounding mode is RARS-specific syntax; drop it if targeting a stricter assembler.
- **Stack frame size**: `det_approx` saves four `s` registers so its frame is 32 bytes; the other two subroutines use 8-byte frames.
