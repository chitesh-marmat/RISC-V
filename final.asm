.data
# storing matrix dimension N = 6
N: .word 6

# 6x6 input matrix A stored in row-major order
 
 
 
A: .float
  1.0,  2.0,  3.0,  0.5,  1.5,  2.5,
  0.5,  1.0,  1.5,  0.2,  0.8,  1.1,
  1.5,  3.0,  4.5,  0.7,  2.3,  3.6,
  2.0, -1.0,  0.5,  1.2,  0.4,  0.9,
 -0.5,  0.3,  1.1,  0.6,  1.7, -0.2,
  0.8,  1.4, -0.3,  0.9,  0.5,  1.2


X: .float 1.0, 2.0, 3.0, 4.0, 5.0, 6.0

# output vector Y (size 6 floats = 24 bytes)
 
Y: .space 24  

# storage for determinant value
det_result: .float 0.0

# boolean flag: 1 if sum(Y) < det(A), else 0
comparison_flag: .word 0

# floating-point constants  
one:     .float 1.0      
neg_one: .float -1.0       

.text
.globl main
main:
    # loading base address of A into a0
    la   a0, A
    # loading base address of X into a1
    la   a1, X
    # loading base address of Y into a2
    la   a2, Y
    # t2 gets address of N
    la   t2, N   
    # loading N                
    lw   t2, 0(t2) 
    # passing N to matvec_mul                
    mv   a3, t2                   

    # computing Y = A * X
    jal  ra, matvec_mul

    # printing Y vector (reload Y base)
    la   s1, Y
    mv   t2, a3                    # reloading N
    li   t0, 0                     # starting row index

print_loop:
    bge  t0, t2, end_print         # stopping if row == N
    
    slli a4, t0, 2                 # compute offset row*4
    add  a4, s1, a4                # computing Y[row]
    flw  fa0, 0(a4)                # loading Y[row] for printing
    li   a7, 2                     # print float
    ecall
    
    li   a0, 32                    # print space(for printing 6x1 output vector)
    li   a7, 11
    ecall

    addi t0, t0, 1                 # row++
    j    print_loop

end_print:
    li   a0, 10                    # newline after vector Y output
    li   a7, 11
    ecall

    # computing determinant det(A)
    jal  ra, det_approx

    # printing determinant
    la   a4, det_result
    flw  fa0, 0(a4)
    li   a7, 2
    ecall

    li   a0, 10
    li   a7, 11
    ecall

    # comparing sum(Y) < det
    jal  ra, compare_sum_det

    la   a0, comparison_flag       # load result flag
    lw   a0, 0(a0)
    li   a7, 1                     # print integer
    ecall

    li   a0, 10
    li   a7, 11
    ecall

    li   a7, 10                    # exit
    ecall

 
#Computing  Y[row] = Σ A[row][col] * X[col]
.globl matvec_mul
matvec_mul:
    addi sp, sp, -8
    sw   s0, 4(sp)                 # saving s0
    sw   s1, 0(sp)                 # saving s1

    # assigning s0 = A base
    mv   s0, a0  
    # assigning s1 = Y base                   
    mv   s1, a2     
    # assigning a5 = X base               
    mv   a5, a1    
    # t2 = N                
    mv   t2, a3                    

    li   t0, 0                     # row index = 0

outer_mv:
    bge  t0, t2, mv_end            # stopping if row == N

    fmv.s.x ft3, x0                # accumulator ft3 = 0.0

    li   t1, 0                     # column index = 0
inner_mv:
    bge  t1, t2, store_mv          # stopping if col == N

    # computing address of A[row][col]
    mul  a4, t0, t2
    add  a4, a4, t1
    slli a4, a4, 2
    add  a0, s0, a4                # a0 = A[row][col]
    flw  ft0, 0(a0)                # ft0 = A[row][col]

    # loading X[col]
    slli a1, t1, 2
    add  a1, a5, a1
    flw  ft1, 0(a1)                # ft1 = X[col]

    # multiplying A[row][col] with X[col]
    
    fmul.s ft2, ft0, ft1           # ft2 = A[row][col] * X[col]

    # accumulating into row sum
    fadd.s ft3, ft3, ft2           # ft3 += ft2

    addi t1, t1, 1                 # col++
    j    inner_mv

store_mv:
    # storing completed dot-product into Y[row]
    slli a4, t0, 2
    add  a4, s1, a4                # &Y[row]
    fmv.s fa0, ft3
    fsw  fa0, 0(a4)

    addi t0, t0, 1                 # row++
    j    outer_mv

mv_end:
    lw   s1, 0(sp)                 # restoring s1
    lw   s0, 4(sp)                 # restoring s0
    addi sp, sp, 8
    jr   ra


 
# LU decomposition with partial pivoting
# then computing determinant as product of diagonal pivots
 
.globl det_approx
det_approx:
    addi sp, sp, -8
    sw   s0, 4(sp)
    sw   s1, 0(sp)

    # assigning s0 = A base
    la   s0, A                     
    la   a3, N
    lw   t2, 0(a3)                 # loading N

    # sign correction factor
    li   a5, 1                     

    la   a0, one
     # determinant accumulator ft9 = 1.0
    flw  ft9, 0(a0)               

    # pivot row i = 0
    li   t0, 0                     

pivot_outer:
    bge  t0, t2, finish_det    # stopping if i == N

    mv   a1, t0             # assuming max_row = i (initially)

    # loading pivot A[i][i]
    mul  a2, t0, t2
    add  a2, a2, t0
    slli a2, a2, 2
    add  a4, s0, a2
    flw  ft0, 0(a4)                # ft0 = pivot

    # computing |pivot|
    fmv.x.w a2, ft0
    li      a3, 0x7fffffff
    and     a2, a2, a3
    fmv.w.x ft1, a2                # ft1 = |A[i][i]|

    addi t1, t0, 1                 # to start searching next row

pivot_search:
    bge  t1, t2, pivot_swap        # for stopping pivot search at end

    # loading A[t1][i]
    mul  a2, t1, t2
    add  a2, a2, t0
    slli a2, a2, 2
    add  a4, s0, a2
    flw  ft2, 0(a4)

    # computing |A[t1][i]|
    fmv.x.w a2, ft2
    and     a2, a3, a2
    fmv.w.x ft3, a2

    # checking if |A[t1][i]| > |pivot|
    
    flt.s a2, ft1, ft3             # if ft1 < ft3, update pivot
    beq  a2, x0, no_pivot_upd
    
    # updating max_row
    mv   a1, t1     
    # updating max pivot magnitude               
    fmv.s ft1, ft3                 

no_pivot_upd:
    addi t1, t1, 1
    j    pivot_search

pivot_swap:
    beq  a1, t0, no_swap           # no swap if best row = current row

    li   s1, 0                     # column index

swap_loop:
    bge  s1, t2, swap_done         # stopping swap if column == N

    # A[i][s1] address
    mul  a2, t0, t2
    add  a2, a2, s1
    slli a2, a2, 2
    add  a0, s0, a2
    flw  ft4, 0(a0)

    # A[a1][s1] address
    mul  a3, a1, t2
    add  a3, a3, s1
    slli a3, a3, 2
    add  a4, s0, a3
    flw  ft5, 0(a4)

    # swapping row i and row a1
    fsw  ft4, 0(a4)
    fsw  ft5, 0(a0)

    addi s1, s1, 1
    j    swap_loop

swap_done:
    # flip determinant sign from row swap
    li   a2, -1
    mul  a5, a5, a2

    la   a4, neg_one
    flw  ft8, 0(a4)
    fmul.s ft9, ft9, ft8           # multiply det by -1

no_swap:
    # reload pivot = A[i][i]
    mul  a2, t0, t2
    add  a2, a2, t0
    slli a2, a2, 2
    add  a0, s0, a2
    flw  ft0, 0(a0)

    # multiply pivot into determinant accumulator
    fmul.s ft9, ft9, ft0

    # eliminating rows below pivot row
    addi t1, t0, 1

elim_rows:
    bge  t1, t2, next_i            # stopping if k == N

    # factor = A[k][i] / pivot
    mul  a2, t1, t2
    add  a2, a2, t0
    slli a2, a2, 2
    add  a0, s0, a2
    flw  ft4, 0(a0)
    fdiv.s ft5, ft4, ft0           # elimination factor

    addi a1, t0, 1                 # start column j = i+1

elim_cols:
    bge  a1, t2, elim_done         # stop if j == N

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

    # A[k][j] -= factor * A[i][j]
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
    # storing final determinant
    fmv.s fa0, ft9
    la   a0, det_result
    fsw  fa0, 0(a0)

    lw   s1, 0(sp)
    lw   s0, 4(sp)
    addi sp, sp, 8
    jr   ra


 
# To compare sum v/s determinant
# Computes sum(Y) and checks if sum(Y) < det(A)
 
.globl compare_sum_det
compare_sum_det:
    addi sp, sp, -8
    sw   s0, 4(sp)                  
    sw   s1, 0(sp)                 

    la   s1, Y                     # assigning s1 = base address of Y vector
    la   a0, N
    lw   t2, 0(a0)                 # t2 = N (number of elements in Y)

    fmv.s.x fa0, x0                # fa0 = 0.0 (initializing running sum of Y)
    li   t0, 0                     # t0 = 0 (loop index for Y)

sum_loop:
    bge  t0, t2, sum_done          # stopping if t0 == N (finished summing all Y[i])

    slli a1, t0, 2                 # compute byte offset = index * 4
    add  a1, s1, a1                # compute address of Y[index]
    flw  ft0, 0(a1)                # loading Y[index] into ft0

    fadd.s fa0, fa0, ft0           # sum = sum + Y[index]

    addi t0, t0, 1                 # index++ to move to next element
    j    sum_loop                  # repeating until all Y elements are summed

sum_done:
    la   a2, det_result            # to load address of stored determinant
    flw  ft1, 0(a2)                # ft1 = determinant value

    # comparing sum(Y) < determinant
    # flt.s sets a0 = 1 if fa0 < ft1, else 0
    flt.s a0, fa0, ft1

    la   a1, comparison_flag
    sw   a0, 0(a1)                 # storing result flag (1 or 0)

    lw   s1, 0(sp)                 # restoring s1
    lw   s0, 4(sp)                 # restoring s0
    addi sp, sp, 8                 # restoring stack pointer
    jr   ra                        # return to caller