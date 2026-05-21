# final_submission_constraint_ok.S
# All integer loop/address registers restricted to x5..x15.
# Floating ops use RV32F single-precision instructions.

.data
N: .word 6

A: .float
  1.25, -2.40,  3.75,  0.50, -1.10,  4.20,
 -0.80,  1.10,  0.95, -3.30,  2.75,  1.60,
  0.60, -1.70,  2.20,  1.30, -0.90,  0.40,
  4.70, -2.20,  0.00,  1.20,  3.30,  0.45,
 -1.50,  2.80,  1.75, -0.60,  2.25, -3.10,
  0.95, -1.40,  4.60,  3.20, -0.80,  1.75

X: .float 1.0, 2.0, 3.0, 4.0, 5.0, 6.0
Y: .space 24

det_result: .float 0.0
comparison_flag: .word 0

one:     .float 1.0
neg_one: .float -1.0

.text
.globl main
main:
    # prepare pointers and N
    la   a0, A
    la   a1, X
    la   a2, Y
    la   t2, N
    lw   t2, 0(t2)      # t2 = N
    mv   a3, t2         # a3 = N for calls

    # compute Y = A * X
    jal  ra, matvec_mul

    # print Y (reload Y base to avoid caller-saved clobber)
    la   s1, Y          # s1 = Y base (x9)
    mv   t2, a3         # t2 = N
    li   t0, 0

print_loop:
    bge  t0, t2, end_print

    slli a4, t0, 2
    add  a4, s1, a4
    flw  fa0, 0(a4)
    li   a7, 2
    ecall

    li   a0, 32
    li   a7, 11
    ecall

    addi t0, t0, 1
    j    print_loop

end_print:
    li   a0, 10
    li   a7, 11
    ecall

    # determinant
    jal  ra, det_approx

    la   a4, det_result
    flw  fa0, 0(a4)
    li   a7, 2
    ecall

    li   a0, 10
    li   a7, 11
    ecall

    # compare sum(Y) < det
    jal  ra, compare_sum_det

    la   a0, comparison_flag
    lw   a0, 0(a0)
    li   a7, 1
    ecall

    li   a0, 10
    li   a7, 11
    ecall

    li   a7, 10
    ecall


############################################################
# matvec_mul: Y = A * X
# Uses: t0,t1,t2 (x5-x7), s0,s1 (x8-x9), a0..a5 (x10-x15)
# Saves s0,s1 on stack; uses a5 as X-base temporary.
############################################################
.globl matvec_mul
matvec_mul:
    addi sp, sp, -8
    sw   s0, 4(sp)
    sw   s1, 0(sp)

    mv   s0, a0      # s0 = A base
    mv   s1, a2      # s1 = Y base
    mv   a5, a1      # a5 = X base (use a5 locally)
    mv   t2, a3      # t2 = N

    li   t0, 0       # row = 0

outer_mv:
    bge  t0, t2, mv_end

    fmv.s.x ft3, x0  # accumulator = 0.0

    li   t1, 0       # col = 0
inner_mv:
    bge  t1, t2, store_mv

    # load A[row][col]
    mul  a4, t0, t2
    add  a4, a4, t1
    slli a4, a4, 2
    add  a0, s0, a4
    flw  ft0, 0(a0)

    # load X[col]
    slli a1, t1, 2
    add  a1, a5, a1
    flw  ft1, 0(a1)

    fmul.s ft2, ft0, ft1
    fadd.s ft3, ft3, ft2

    addi t1, t1, 1
    j    inner_mv

store_mv:
    fmv.s fa0, ft3
    slli a4, t0, 2
    add  a4, s1, a4
    fsw  fa0, 0(a4)

    addi t0, t0, 1
    j    outer_mv

mv_end:
    lw   s1, 0(sp)
    lw   s0, 4(sp)
    addi sp, sp, 8
    jr   ra


############################################################
# det_approx: LU with partial pivoting
# Use only x5..x15 for integer work; s0 holds A base, s1 used as loop temp
# a0..a5 used as temporary integer/pointer registers
# Floating temporaries: ft0..ft9 as needed.
############################################################
.globl det_approx
det_approx:
    addi sp, sp, -8
    sw   s0, 4(sp)
    sw   s1, 0(sp)

    la   s0, A
    la   a3, N
    lw   t2, 0(a3)       # t2 = N

    li   a5, 1           # sign stored in a5 (x15)

    la   a0, one
    flw  ft9, 0(a0)      # ft9 = 1.0 (det accumulator)

    li   t0, 0           # i = 0

pivot_outer:
    bge  t0, t2, finish_det

    mv   a1, t0          # a1 = max_row

    # load pivot A[i][i]
    mul  a2, t0, t2
    add  a2, a2, t0
    slli a2, a2, 2
    add  a4, s0, a2
    flw  ft0, 0(a4)

    # ft1 = |pivot|
    fmv.x.w a2, ft0
    li      a3, 0x7fffffff
    and     a2, a2, a3
    fmv.w.x ft1, a2

    addi t1, t0, 1
pivot_search:
    bge  t1, t2, pivot_swap

    mul  a2, t1, t2
    add  a2, a2, t0
    slli a2, a2, 2
    add  a4, s0, a2
    flw  ft2, 0(a4)

    fmv.x.w a2, ft2
    and     a2, a2, a3
    fmv.w.x ft3, a2

    flt.s a2, ft1, ft3
    beq  a2, x0, no_pivot_upd
    mv   a1, t1
    fmv.s ft1, ft3
no_pivot_upd:
    addi t1, t1, 1
    j    pivot_search

pivot_swap:
    beq  a1, t0, no_swap

    li   s1, 0            # s1 = col index for swapping

swap_loop:
    bge  s1, t2, swap_done

    # addr = &A[i][s1]
    mul  a2, t0, t2
    add  a2, a2, s1
    slli a2, a2, 2
    add  a0, s0, a2
    flw  ft4, 0(a0)

    # addr2 = &A[a1][s1]
    mul  a3, a1, t2
    add  a3, a3, s1
    slli a3, a3, 2
    add  a4, s0, a3
    flw  ft5, 0(a4)

    fsw  ft4, 0(a4)
    fsw  ft5, 0(a0)

    addi s1, s1, 1
    j    swap_loop

swap_done:
    li   a2, -1
    mul  a5, a5, a2      # flip sign

    la   a4, neg_one
    flw  ft8, 0(a4)
    fmul.s ft9, ft9, ft8

no_swap:

    # reload pivot
    mul  a2, t0, t2
    add  a2, a2, t0
    slli a2, a2, 2
    add  a0, s0, a2
    flw  ft0, 0(a0)

    fmul.s ft9, ft9, ft0

    addi t1, t0, 1
elim_rows:
    bge  t1, t2, next_i

    # factor = A[k][i] / pivot
    mul  a2, t1, t2
    add  a2, a2, t0
    slli a2, a2, 2
    add  a0, s0, a2
    flw  ft4, 0(a0)

    fdiv.s ft5, ft4, ft0

    addi a1, t0, 1      # j = i+1
elim_cols:
    bge  a1, t2, elim_done

    # A[i][j]
    mul  a2, t0, t2
    add  a2, a2, a1
    slli a2, a2, 2
    add  a0, s0, a2
    flw  ft7, 0(a0)

    # A[k][j]
    mul  a2, t1, t2
    add  a2, a2, a1
    slli a2, a2, 2
    add  a4, s0, a2
    flw  ft8, 0(a4)

    fmul.s ft6, ft5, ft7
    fsub.s ft8, ft8, ft6
    fsw  ft8, 0(a4)

    addi a1, a1, 1
    j    elim_cols

elim_done:
    addi t1, t1, 1
    j    elim_rows

next_i:
    addi t0, t0, 1
    j    pivot_outer

finish_det:
    fmv.s fa0, ft9
    la   a0, det_result
    fsw  fa0, 0(a0)

    lw   s1, 0(sp)
    lw   s0, 4(sp)
    addi sp, sp, 8
    jr   ra


############################################################
# compare_sum_det: sum Y and compare to determinant
############################################################
.globl compare_sum_det
compare_sum_det:
    addi sp, sp, -8
    sw   s0, 4(sp)
    sw   s1, 0(sp)

    la   s1, Y
    la   a0, N
    lw   t2, 0(a0)

    fmv.s.x fa0, x0
    li   t0, 0

sum_loop:
    bge  t0, t2, sum_done

    slli a1, t0, 2
    add  a1, s1, a1
    flw  ft0, 0(a1)
    fadd.s fa0, fa0, ft0

    addi t0, t0, 1
    j    sum_loop

sum_done:
    la   a2, det_result
    flw  ft1, 0(a2)
    flt.s a0, fa0, ft1
    la   a1, comparison_flag
    sw   a0, 0(a1)

    lw   s1, 0(sp)
    lw   s0, 4(sp)
    addi sp, sp, 8
    jr   ra