.data
N: .word 6

A: .float  1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
           2.0, 1.0, 2.0, 1.0, 2.0, 1.0,
           3.0, 0.0, 1.0, 2.0, 1.0, 0.0,
           1.0, 2.0, 3.0, 4.0, 3.0, 2.0,
           2.0, 1.0, 0.0, 1.0, 2.0, 1.0,
           1.0, 3.0, 2.0, 4.0, 1.0, 5.0

X: .float 1.0, 2.0, 3.0, 4.0, 5.0, 6.0
Y: .space 24

det_result: .float 0.0
comparison_flag: .word 0

one:      .float  1.0
neg_one:  .float -1.0
eps_const: .float 1.0e-7


############################################################
.text
.globl main
main:
    # call matvec_mul(A, X, Y, N)
    la   a0, A
    la   a1, X
    la   a2, Y
    la   t0, N
    lw   a3, 0(t0)
    jal  ra, matvec_mul

    # print Y
    la   s2, Y
    la   t0, N
    lw   t1, 0(t0)
    li   t2, 0

print_loop:
    bge  t2, t1, phase2
    slli t3, t2, 2
    add  t3, s2, t3
    flw  fa0, 0(t3)
    li   a7, 2
    ecall

    li   a7, 11
    li   a0, 10
    ecall

    addi t2, t2, 1
    j    print_loop


phase2:
    # compute determinant
    jal  ra, det_approx

    # print determinant
    la   t0, det_result
    flw  fa0, 0(t0)
    li   a7, 2
    ecall

    li   a7, 11
    li   a0, 10
    ecall

    # compare sum(Y) < det ?
    jal  ra, compare_sum_det

    la   t5, comparison_flag
    lw   a0, 0(t5)
    li   a7, 1
    ecall

    li   a7, 11
    li   a0, 10
    ecall

    li   a7, 10
    ecall


############################################################
# matvec_mul (your original code, untouched)
############################################################
.globl matvec_mul
matvec_mul:
    addi sp, sp, -8
    sw   s0, 4(sp)
    sw   s1, 0(sp)

    mv   a4, a0
    mv   a5, a1
    mv   s0, a2
    mv   s1, a3

    li   t0, 0

outer_loop:
    bge  t0, s1, matvec_end

    mul  t2, t0, s1
    fmv.s.x ft3, x0

    li   t1, 0
inner_loop:
    bge  t1, s1, store_row

    add  t3, t2, t1
    slli t3, t3, 2
    add  t4, a4, t3
    flw  ft0, 0(t4)

    slli t5, t1, 2
    add  t5, a5, t5
    flw  ft1, 0(t5)

    fmul.s ft2, ft0, ft1
    fadd.s ft3, ft3, ft2

    addi t1, t1, 1
    j    inner_loop

store_row:
    fmv.s fa0, ft3
    slli t6, t0, 2
    add  t6, s0, t6
    fsw  fa0, 0(t6)

    addi t0, t0, 1
    j    outer_loop

matvec_end:
    lw   s1, 0(sp)
    lw   s0, 4(sp)
    addi sp, sp, 8
    jr   ra


############################################################
# det_approx — CLEAN, WORKING, NO t7
############################################################
.globl det_approx
det_approx:
    addi sp, sp, -32
    sw   s0, 28(sp)
    sw   s1, 24(sp)
    sw   s2, 20(sp)
    sw   s3, 16(sp)

    la   s0, A
    la   t0, N
    lw   t1, 0(t0)

    li   t6, 1            # sign = +1
    li   t2, 0            # i = 0

pivot_outer:
    bge  t2, t1, compute_diag

    mv   t3, t2            # max_row = i

    # pivot A[i,i]
    mul  a2, t2, t1
    add  a2, a2, t2
    slli a2, a2, 2
    add  a3, s0, a2
    flw  ft0, 0(a3)

    # abs(pivot) -> ft1
    fmv.x.w a4, ft0
    li      a5, 0x7fffffff
    and     a4, a4, a5
    fmv.w.x ft1, a4

    # search pivot row
    addi s1, t2, 1
pivot_search:
    bge  s1, t1, pivot_swap

    mul  a2, s1, t1
    add  a2, a2, t2
    slli a2, a2, 2
    add  a3, s0, a2
    flw  ft2, 0(a3)

    fmv.x.w a4, ft2
    and     a4, a4, a5
    fmv.w.x ft3, a4

    flt.s a4, ft1, ft3
    beq   a4, x0, no_up
    mv    t3, s1
    fmv.s ft1, ft3
no_up:
    addi s1, s1, 1
    j    pivot_search

pivot_swap:
    beq  t3, t2, no_swap

    li   a6, 0
swap_loop:
    bge  a6, t1, swap_done

    mul  a2, t2, t1
    add  a2, a2, a6
    slli a2, a2, 2
    add  a3, s0, a2
    flw  ft4, 0(a3)

    mul  a2, t3, t1
    add  a2, a2, a6
    slli a2, a2, 2
    add  a5, s0, a2
    flw  ft5, 0(a5)

    fsw  ft4, 0(a5)
    fsw  ft5, 0(a3)

    addi a6, a6, 1
    j    swap_loop

swap_done:
    li   a6, -1
    mul  t6, t6, a6

no_swap:
    # reload new pivot
    mul  a2, t2, t1
    add  a2, a2, t2
    slli a2, a2, 2
    add  a3, s0, a2
    flw  ft0, 0(a3)

    # abs pivot in ft2
    fmv.x.w a4, ft0
    and     a4, a4, a5
    fmv.w.x ft2, a4

    la   a7, eps_const
    flw  ft6, 0(a7)
    #flt.s a4, ft2, ft6
    #bne   a4, x0, det_zero

    # elimination rows
    addi s1, t2, 1
elim_outer:
    bge  s1, t1, next_i

    mul  a2, s1, t1
    add  a2, a2, t2
    slli a2, a2, 2
    add  a3, s0, a2
    flw  ft4, 0(a3)

    fdiv.s ft5, ft4, ft0

    mv   a6, t2
elim_inner:
    bge  a6, t1, elim_next

    mul  a2, t2, t1
    add  a2, a2, a6
    slli a2, a2, 2
    add  a3, s0, a2
    flw  ft7, 0(a3)

    mul  a2, s1, t1
    add  a2, a2, a6
    slli a2, a2, 2
    add  a4, s0, a2
    flw  ft8, 0(a4)

    fmul.s ft9, ft5, ft7
    fsub.s ft8, ft8, ft9
    fsw    ft8, 0(a4)

    addi a6, a6, 1
    j    elim_inner

elim_next:
    addi s1, s1, 1
    j    elim_outer

next_i:
    addi t2, t2, 1
    j    pivot_outer

det_zero:
    fmv.s.x fa0, x0
    la   a3, det_result
    fsw  fa0, 0(a3)
    j    finish_det


############################################################
# multiply diagonal after elimination
############################################################
compute_diag:
    la   a3, one
    flw  ft0, 0(a3)

    li   t2, 0
diag_loop:
    bge  t2, t1, diag_done

    mul  a2, t2, t1
    add  a2, a2, t2
    slli a2, a2, 2
    add  a3, s0, a2
    flw  ft1, 0(a3)
    fmul.s ft0, ft0, ft1

    addi t2, t2, 1
    j    diag_loop

diag_done:
    bltz t6, diag_neg
    fmv.s fa0, ft0
    la   a3, det_result
    fsw  fa0, 0(a3)
    j    finish_det

diag_neg:
    la   a3, neg_one
    flw  ft2, 0(a3)
    fmul.s ft0, ft0, ft2
    fmv.s fa0, ft0
    la   a3, det_result
    fsw  fa0, 0(a3)

finish_det:
    lw   s3, 16(sp)
    lw   s2, 20(sp)
    lw   s1, 24(sp)
    lw   s0, 28(sp)
    addi sp, sp, 32
    jr   ra


############################################################
# compare_sum_det (unchanged)
############################################################
.globl compare_sum_det
compare_sum_det:
    addi sp, sp, -8
    sw   s0, 4(sp)
    sw   s1, 0(sp)

    la   s0, Y
    la   t0, N
    lw   t1, 0(t0)

    fmv.s.x fa0, x0
    li   s1, 0
sum_loop:
    bge  s1, t1, sum_done
    slli t2, s1, 2
    add  t2, s0, t2
    flw  ft0, 0(t2)
    fadd.s fa0, fa0, ft0
    addi s1, s1, 1
    j    sum_loop

sum_done:
    la   t2, det_result
    flw  ft1, 0(t2)
    flt.s a0, fa0, ft1
    la   t3, comparison_flag
    sw   a0, 0(t3)

    lw   s1, 0(sp)
    lw   s0, 4(sp)
    addi sp, sp, 8
    jr   ra