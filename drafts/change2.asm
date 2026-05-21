.data
# Matrix A (2x2)
A: .float 1.0, 2.0,
          3.0, 4.0

# Vector X (2x1)
X: .float 5.0, 6.0

# Output Vector Y
Y: .space 8           # 2 floats (8 bytes)

# Determinant
det_result: .float 0.0

.text
.globl main

########################################################
# main: Calls subroutines for matvec and determinant
########################################################
main:
    # Load addresses into argument registers
    la a0, A           # arg1: base of matrix A
    la a1, X           # arg2: base of vector X
    la a2, Y           # arg3: base of result vector Y

    # Call matvec_mul(A, X, Y)
    jal ra, matvec_mul

    # Call det_comp(A)
    la a0, A
    jal ra, det_comp

    # Store determinant into memory
    la t0, det_result
    fsw fa0, 0(t0)

    # Print results
    la t1, Y
    flw fa0, 0(t1)     # Y[0]
    li a7, 2
    ecall

    li a7, 11          # newline
    li a0, 10
    ecall

    flw fa0, 4(t1)     # Y[1]
    li a7, 2
    ecall

    li a7, 11          # newline
    li a0, 10
    ecall

    la t2, det_result
    flw fa0, 0(t2)     # determinant
    li a7, 2
    ecall

    # Exit
    li a7, 10
    ecall


########################################################
# matvec_mul(A, X, Y)
# Computes: Y = A * X for 2x2 matrix
# a0 = base of A, a1 = base of X, a2 = base of Y
########################################################
matvec_mul:
    # Y[0] = A[0]*X[0] + A[1]*X[1]
    flw f1, 0(a0)          # A[0][0]
    flw f2, 0(a1)          # X[0]
    fmul.s f3, f1, f2

    flw f4, 4(a0)          # A[0][1]
    flw f5, 4(a1)          # X[1]
    fmul.s f6, f4, f5

    fadd.s f7, f3, f6
    fsw f7, 0(a2)          # Y[0]

    # Y[1] = A[2]*X[0] + A[3]*X[1]
    flw f1, 8(a0)
    flw f2, 0(a1)
    fmul.s f3, f1, f2

    flw f4, 12(a0)
    flw f5, 4(a1)
    fmul.s f6, f4, f5

    fadd.s f7, f3, f6
    fsw f7, 4(a2)          # Y[1]

    ret


########################################################
# det_comp(A)
# Computes determinant of 2x2 matrix
# a0 = base of A
# returns det(A) in fa0
########################################################
det_comp:
    flw f1, 0(a0)          # A00
    flw f2, 12(a0)         # A11
    fmul.s f3, f1, f2      # A00 * A11

    flw f4, 4(a0)          # A01
    flw f5, 8(a0)          # A10
    fmul.s f6, f4, f5      # A01 * A10

    fsub.s fa0, f3, f6     # det = f3 - f6 (return value)
    ret