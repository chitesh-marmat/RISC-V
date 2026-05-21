# matvec_det_optionB_fixed.S
# Option B — Academic-style fixed: all integer temporaries & addresses restricted to x5-x15.
# Careful allocation to avoid overwriting offsets mid-calculation and to prevent out-of-range addresses.

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

one:      .float  1.0
neg_one:  .float -1.0

.text
.globl main
main:
    # call matvec_mul(A, X, Y, N)
    la   a0, A          # a0 = A base (x10)
    la   a1, X          # a1 = X base (x11)
    la   a2, Y          # a2 = Y base (x12)
    la   t2, N          # t2 = addr N (x7 temporarily)
    lw   t2, 0(t2)      # t2 = N (x7)
    mv   a3, t2         # a3 = N (x13) - keep in a3 for call, but we'll pass in a3
    jal  ra, matvec_mul

    # print Y
    mv   s1, a2         # s1 = Y base (x9)
    mv   t2, a3         # t2 = N (x7)
    li   t0, 0          # t0 = index (x5)

print_loop:
    bge  t0, t2, end_y_print

    slli a4, t0, 2      # a4 = offset
    add  a4, s1, a4     # a4 = &Y[t0]
    flw  fa0, 0(a4)
    li   a7, 2          # print float syscall
    ecall

    # print space
    li   a0, 32
    li   a7, 11
    ecall

    addi t0, t0, 1
    j    print_loop

end_y_print:
    li   a0, 10
    li   a7, 11
    ecall

    # compute determinant
    jal  ra, det_approx

    # print determinant
    la   a4, det_result
    flw  fa0, 0(a4)
    li   a7, 2
    ecall

    li   a0, 10
    li   a7, 11
    ecall

    # compare sum(Y) < det ?
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

#############################################################
# matvec_mul — corrected version (X-base = a1)
############################################################
.globl matvec_mul
matvec_mul:
    addi sp, sp, -8
    sw   s0, 4(sp)
    sw   s1, 0(sp)

    mv   s0, a0        # s0 = A base (x8)
    mv   s1, a2        # s1 = Y base (x9)
    mv   t2, a3        # t2 = N (x7)
    # a1 already = X base from main

    li   t0, 0         # row = 0

outer_loop:
    bge  t0, t2, matvec_end

    mul  a4, t0, t2    # a4 = row*N
    fmv.s.x ft3, x0    # accumulator = 0.0

    li   t1, 0         # col = 0
inner_loop:
    bge  t1, t2, store_row

    add  a4, a4, t1
    slli a4, a4, 2
    add  a5, s0, a4    # A[row*N+col]
    flw  ft0, 0(a5)

    slli a4, t1, 2
    add  a4, a1, a4    # X[col]
    flw  ft1, 0(a4)

    fmul.s ft2, ft0, ft1
    fadd.s ft3, ft3, ft2

    addi t1, t1, 1
    mul  a4, t0, t2     # restore a4 = row*N
    j    inner_loop

store_row:
    fmv.s fa0, ft3
    slli a4, t0, 2
    add  a4, s1, a4
    fsw  fa0, 0(a4)

    addi t0, t0, 1
    j    outer_loop

matvec_end:
    lw   s1, 0(sp)
    lw   s0, 4(sp)
    addi sp, sp, 8
    jr   ra

############################################################
# det_approx — LU decomposition with partial pivoting
# Register plan (strict x5-x15 usage):
# t0 (x5): pivot i
# t1 (x6): inner counters (k/j) reused carefully
# t2 (x7): N
# s0 (x8): A base
# s1 (x9): scratch / Y base earlier not used here
# a0 (x10), a1 (x11), a2 (x12), a3 (x13), a4 (x14), a5 (x15): temps/offsets
#
# Note: a7 (x17) is NOT used here (only for syscalls).
############################################################
.globl det_approx
det_approx:
    addi sp, sp, -40
    sw   s0, 36(sp)
    sw   s1, 32(sp)
    sw   s2, 28(sp)
    sw   s3, 24(sp)
    sw   s4, 20(sp)

    la   s0, A          # s0 = A base (x8)
    la   a3, N
    lw   t2, 0(a3)      # t2 = N (x7)

    li   a5, 1          # a5 = sign (x15)

    la   a0, one
    flw  ft9, 0(a0)     # ft9 = 1.0 (det accumulator)

    li   t0, 0          # pivot i = 0

pivot_outer:
    bge  t0, t2, finish_det

    mv   a1, t0         # a1 = candidate max_row (x11)

    # load A[i][i]
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

    addi t1, t0, 1      # t1 = k = i+1
pivot_search:
    bge  t1, t2, pivot_swap

    # load A[k][i]
    mul  a2, t1, t2
    add  a2, a2, t0
    slli a2, a2, 2
    add  a4, s0, a2
    flw  ft2, 0(a4)

    # ft3 = |A[k][i]|
    fmv.x.w a2, ft2
    and     a2, a2, a3
    fmv.w.x ft3, a2

    # update pivot if ft3 > ft1
    flt.s a2, ft1, ft3
    beq  a2, x0, no_pivot_upd
    mv   a1, t1
    fmv.s ft1, ft3
no_pivot_upd:
    addi t1, t1, 1
    j    pivot_search

pivot_swap:
    beq  a1, t0, no_swap

    li   a2, 0          # a2 = col index for swapping
swap_loop:
    bge  a2, t2, swap_done

    # addr_i = (i*N + a2) * 4
    mul  a4, t0, t2
    add  a4, a4, a2
    slli a4, a4, 2
    add  a6, s0, a4
    flw  ft4, 0(a6)

    # addr_max = (max_row*N + a2) * 4
    mul  a4, a1, t2
    add  a4, a4, a2
    slli a4, a4, 2
    add  a7, s0, a4
    flw  ft5, 0(a7)

    fsw  ft4, 0(a7)
    fsw  ft5, 0(a6)

    addi a2, a2, 1
    j    swap_loop

swap_done:
    li   a2, -1
    mul  a5, a5, a2     # flip sign in a5

    la   a4, neg_one
    flw  ft8, 0(a4)
    fmul.s ft9, ft9, ft8

no_swap:

    # reload pivot = A[i][i] (after swap)
    mul  a4, t0, t2
    add  a4, a4, t0
    slli a4, a4, 2
    add  a6, s0, a4
    flw  ft0, 0(a6)

    # accumulate determinant *= pivot
    fmul.s ft9, ft9, ft0

    addi t1, t0, 1      # t1 = k = i+1
elim_rows:
    bge  t1, t2, next_i

    # factor = A[k][i] / pivot
    mul  a4, t1, t2
    add  a4, a4, t0
    slli a4, a4, 2
    add  a6, s0, a4
    flw  ft4, 0(a6)

    fdiv.s ft5, ft4, ft0    # factor in ft5

    addi a1, t0, 1          # a1 = j = i+1
elim_cols:
    bge  a1, t2, elim_done

    # load A[i][j] into ft7 (compute addr_i_j)
    mul  a4, t0, t2
    add  a4, a4, a1
    slli a4, a4, 2
    add  a6, s0, a4
    flw  ft7, 0(a6)

    # load A[k][j] into ft8 (compute addr_k_j)
    mul  a4, t1, t2
    add  a4, a4, a1
    slli a4, a4, 2
    add  a7, s0, a4
    flw  ft8, 0(a7)

    # A[k][j] -= factor * A[i][j]
    fmul.s ft6, ft5, ft7
    fsub.s ft8, ft8, ft6
    fsw   ft8, 0(a7)

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
    la   a4, det_result
    fsw  fa0, 0(a4)

    lw   s4, 20(sp)
    lw   s3, 24(sp)
    lw   s2, 28(sp)
    lw   s1, 32(sp)
    lw   s0, 36(sp)
    addi sp, sp, 40
    jr   ra

############################################################
# compare_sum_det
# computes sum(Y) in fa0, compares to det_result, stores 1/0 in comparison_flag
# uses t0..t2, a0..a2, s1.
############################################################
.globl compare_sum_det
compare_sum_det:
    addi sp, sp, -8
    sw   s0, 4(sp)
    sw   s1, 0(sp)

    la   s1, Y
    la   a0, N
    lw   t2, 0(a0)       # t2 = N

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
    