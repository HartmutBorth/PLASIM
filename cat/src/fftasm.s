# Fast Double Precision Fourier Transformation
# ============================================
# E. Kirk - 15-Apr-2015
# ============================================
# call ifft2s(a  ,c  ,trigs,n  ,lot)
#            (rdi,rsi,rdx  ,rcx,r8 )
#
# rdi: a(n,lot)  : real(8) source array
# rsi: c(n,lot)  : real(8) target array
# rdx: trigs(n)  : sine/cosine values
# rcx: n         : 1st. dimension
# r8 : lot       : 2nd. dimension
# --------------------------------------------

   .globl _ifft2s_
_ifft2s_:

   .globl ifft2s_
ifft2s_:

   pushq   %r13
   pushq   %r12
   movl    (%r8 ), %r8d              # r8    = lot
   movl    (%rcx), %ecx              # rcx   = n
   movq    $0x3fe0000000000000, %rax # rax   = 0.5
   movd    %rax, %xmm8               # xmm8  = 0.5
   movq    %rcx, %r9                 # r9    = n
   shrq    $2, %r9                   # r9    = n / 4
   dec     %r9                       # r9    = (n / 4) - 1
   movq    %r9 , %r10                # r10   = loop count
   movq    %rdx, %r11                # r11   = &trigs
   movq    %rdi, %r13                # r13   = &a

LOOP_LOT_2S:
   movapd  %xmm8 , %xmm0             # xmm0  = 0.5
   mulsd   (%rdi), %xmm0             # xmm0  = 0.5 * a(0)
   movsd   %xmm0 , (%rsi)            # c(0)  = xmm0
   movsd   %xmm0 ,8(%rsi)            # c(1)  = xmm0
   leaq    -16(%rdi,%rcx,8), %r12    # r12   = &a(n-2) : ib
   addq    $16, %rdi                 # rdi   = &a(2)   : ia
   addq    $16, %rsi                 # rsi   = &c(2)   : j

LOOP_IFFT2:
   addq    $16, %rdx                 # &trigs
   movsd   (%rdi), %xmm7             # xmm7  = a(ia)
   movsd   (%r12), %xmm3             # xmm3  = a(ib)
   movsd   8(%r12), %xmm6            # xmm6  = a(ib+1)
   movapd  %xmm7, %xmm1              # xmm1  = a(ia)
   movsd   8(%rdi), %xmm2            # xmm2  = a(ia+1)
   subsd   %xmm3, %xmm1              # xmm1  = a(ia) - a(ib)
   addsd   %xmm7, %xmm3              # xmm3  = a(ia) + a(ib)
   movsd   (%rdx), %xmm4             # xmm4  = trigs(ia)
   movapd  %xmm2, %xmm0              # xmm0  = a(ia+1)
   subsd   %xmm6, %xmm2              # xmm2  = a(ia+1) - a(ib+1)
   movsd   8(%rdx), %xmm5            # xmm5  = trigs(ia+1)
   addsd   %xmm6, %xmm0              # xmm0  = a(ia+1) + a(ib+1)
   movsd   %xmm3, (%rsi)             # c(j)  = xmm3
   movapd  %xmm5, %xmm3              # xmm3  = s1 = trigs(ia+1)
   movsd   %xmm2, 16(%rsi)           # c(j+2) = 
   movapd  %xmm4, %xmm2
   mulsd   %xmm1, %xmm2
   mulsd   %xmm0, %xmm3
   mulsd   %xmm5, %xmm1
   mulsd   %xmm4, %xmm0
   subsd   %xmm3, %xmm2
   addsd   %xmm1, %xmm0
   movsd   %xmm2,  8(%rsi)           # c(j+1)
   movsd   %xmm0, 24(%rsi)           # c(j+3)
   addq    $32, %rsi                 # j  = j  + 4
   addq    $16, %rdi                 # ia = ia + 2
   subq    $16, %r12                 # ib = ib - 2
   dec     %r9
   jnz     LOOP_IFFT2

   movsd   (%rdi), %xmm0             # xmm0   = a(ia)
   movsd   %xmm0, (%rsi)             # c(n-2) = a(ia)
   xorpd   %xmm0, %xmm0              # xmm0   = 0
   subsd   8(%rdi), %xmm0            # xmm0   = -a(ia+1)
   movsd   %xmm0, 8(%rsi)            # c(n-1) = -a(ia+1)

   leaq    (%r13,%rcx,8), %rdi       # advance &a + n
   movq    %rdi, %r13                # save new base
   addq    $16, %rsi                 # advance &c + 2
   movq    %r11, %rdx                # restore &trigs
   movq    %r10, %r9                 # count for inner loop
   dec     %r8                       # lot
   jnz     LOOP_LOT_2S

   popq    %r12
   popq    %r13
   ret


# Fast Double Precision Fourier Transformation
# ============================================
# E. Kirk - 08-Jun-2015
# ============================================
# call ifft4s(a  ,c  ,trigs,n  ,lot)
#            (rdi,rsi,rdx  ,rcx, r8)
#
# rdi: a(n,lot)  : real(8) source array
# rsi: c(n,lot)  : real(8) target array
# rdx: trigs(n)  : sine/cosine values
# rcx: n         : 1st. dimension
# r8 : lot       : 2nd. dimension
# --------------------------------------------

   .globl   _ifft4s_
_ifft4s_:

   .globl   ifft4s_
ifft4s_:

   pushq   %rbp
   movq    %rsp, %rbp
   call    mcount
   pushq   %r15
   pushq   %r14
   pushq   %r13
   pushq   %r12
   pushq   %rbx
   movl    (%r8 ), %r15d             # r15  = lot
   movl    (%rcx), %ecx              # ecx  = n
   movq    %rcx, %rax                # rax  = n
   sarq    %rax                      # rax  = n / 2
   movq    %rcx, %r14                # r14  = n
   sarq    $3, %r14                  # r14  = n / 8
   dec     %r14                      # r14  = n / 8 - 1
   shlq    $3, %rcx                  # rcx  = n * 8

LOOP_O_IFFT4S:

   leaq     16(%rdi)       , %r10    # &a(i0)  i0 = 2
   leaq     16(%rdi,%rax,8), %r11    # &a(i1)  i1 = 2 + n/2
   leaq    -16(%rdi,%rcx)  , %r12    # &a(i2)  i2 = n - 2
   leaq    -32(%r11)       , %r13    # &a(i3)  i3 = i2 - n/2

   movq   $0x3fe0000000000000, %r8   # r8    = 0.5
   movd        %r8    , %xmm6        # xmm6  = 0.5
   movsd      (%rdi)  , %xmm0        # xmm0  = a(0)
   movsd   -16(%r11)  , %xmm4        # xmm4  = a(i1)
   movsd    -8(%r11)  , %xmm5        # xmm5  = a(i1+1)

   mulsd   %xmm6 , %xmm0             # xmm0 = a(0) * 0.5
   movapd  %xmm0 , %xmm1             # xmm1 = a(0) * 0.5
   movapd  %xmm0 , %xmm2             # xmm2 = a(0) * 0.5
   movapd  %xmm0 , %xmm3             # xmm3 = a(0) * 0.5

   addsd   %xmm4 , %xmm0             # xmm0 = a(0) + a(i1)
   subsd   %xmm5 , %xmm1             # xmm1 = a(0) - a(i1+1)
   subsd   %xmm4 , %xmm2             # xmm2 = a(0) - a(i1)
   addsd   %xmm5 , %xmm3             # xmm3 = a(0) + a(i1+1)

   movsd   %xmm0,   (%rsi)           # c(0)   = xmm0
   movsd   %xmm1,  8(%rsi)           # c(1)   = xmm1
   movsd   %xmm2, 16(%rsi)           # c(2)   = xmm2
   movsd   %xmm3, 24(%rsi)           # c(3)   = xmm3

   movq    %rdx, %r8                 # r8   = trigs
   movq    %rdx, %r9                 # r9   = trigs

   pushq   %rdx                      # save trigs
   pushq   %r14                      # save inner loop count

LOOP_I_IFFT4S:                                 

   movsd   (%r10), %xmm2             # xmm2  = a(i0)
   movsd   (%r11), %xmm3             # xmm3  = a(i1)
   movsd   (%r12), %xmm4             # xmm4  = a(i2)
   movsd   (%r13), %xmm5             # xmm5  = a(i3)

   movsd  8(%r10), %xmm8             # xmm8  = a(i0+1)
   movsd  8(%r11), %xmm9             # xmm9  = a(i1+1)
   movsd  8(%r12), %xmm10            # xmm10 = a(i2+1)
   movsd  8(%r13), %xmm11            # xmm11 = a(i3+1)

   addq    $16, %rdx                 # rdx   = &trigs(2)
   addq    $32, %r8                  # r8    = &trigs(4)
   addq    $48, %r9                  # r9    = &trigs(6)
   addq    $64, %rsi                 # rsi   = &c

   movapd  %xmm2 , %xmm0             # xmm0  = a(i0)
   movapd  %xmm3 , %xmm1             # xmm1  = a(i1)
   addsd   %xmm4 , %xmm0             # xmm0  = a(i0) + a(i2)
   addsd   %xmm5 , %xmm1             # xmm1  = a(i1) + a(i3)
   subsd   %xmm4 , %xmm2             # xmm2  = a(i0) - a(i2)
   subsd   %xmm5 , %xmm3             # xmm3  = a(i1) - a(i3)

   movapd  %xmm8 , %xmm6             # xmm6  = a(i0+1)
   movapd  %xmm9 , %xmm7             # xmm7  = a(i1+1)
   addsd   %xmm10, %xmm6             # xmm6  = a(i0+1) + a(i2+1)
   addsd   %xmm11, %xmm7             # xmm7  = a(i1+1) + a(i3+1)
   subsd   %xmm10, %xmm8             # xmm8  = a(i0+1) - a(i2+1)
   subsd   %xmm11, %xmm9             # xmm9  = a(i1+1) - a(i3+1)

   movapd  %xmm0 , %xmm4             # xmm4  = a0p2
   addsd   %xmm1 , %xmm0             # xmm0  = a0p2 + a1p3
   subsd   %xmm1 , %xmm4             # xmm4  = a0p2 - a1p3

   movapd  %xmm8 , %xmm5             # xmm5  = a4m6
   addsd   %xmm9 , %xmm8             # xmm8  = a4m6 + a5m7
   subsd   %xmm9 , %xmm5             # xmm5  = a4m6 - a5m7

   movsd   %xmm0 , -32(%rsi)         # c(0)
   movsd   %xmm8 ,    (%rsi)         # c(4)

   movapd  %xmm2 , %xmm0             # xmm0  = a0m2
   addsd   %xmm7 , %xmm2             # xmm2  = a0m2 + a5p7
   subsd   %xmm7 , %xmm0             # xmm0  = a0m2 - a5p7

   movapd  %xmm6 , %xmm1             # xmm1  = a4p6
   addsd   %xmm3 , %xmm1             # xmm1  = a4p6 + a1m3
   subsd   %xmm3 , %xmm6             # xmm6  = a4p6 - a1m3

   movsd    (%r8), %xmm7             # xmm7  = c2
   movsd    (%r9), %xmm15            # xmm15 = c3
   movsd   8(%r8), %xmm9             # xmm9  = s2
   movsd   8(%r9), %xmm12            # xmm12 = s3

   movapd  %xmm4 , %xmm3             # a0p2m1p3
   mulsd   %xmm7 , %xmm4             # (c2,c3) *
   mulsd   %xmm9 , %xmm3             # (s2,s3) *

   movapd  %xmm2 , %xmm10            # a0m2p5p7
   mulsd   %xmm15, %xmm2             # (c2,c3) *
   mulsd   %xmm12, %xmm10            # (s2,s3) *

   movapd  %xmm5 , %xmm11            # a4m6m5m7
   mulsd   %xmm7 , %xmm5             # (c2,c3) *
   mulsd   %xmm9 , %xmm11            # (s2,s3) *

   movapd  %xmm6 , %xmm14            # a4p6m1m3
   mulsd   %xmm15, %xmm6             # (c2,c3) *
   mulsd   %xmm12, %xmm14            # (s2,s3) *

   subsd   %xmm11, %xmm4             # xmm4  = c(2)
   addsd   %xmm5 , %xmm3             # xmm3  = c(6)

   subsd   %xmm14, %xmm2             # xmm4  = c(3)
   addsd   %xmm6 , %xmm10            # xmm3  = c(7)

   movsd   %xmm4 ,-16(%rsi)          # c(2)
   movsd   %xmm3 , 16(%rsi)          # c(6)
   movsd   %xmm2 , -8(%rsi)          # c(3)
   movsd   %xmm10, 24(%rsi)          # c(7)

   movsd    (%rdx), %xmm7            # xmm7  = c1
   movsd   8(%rdx), %xmm9            # xmm9  = s1

   movapd  %xmm1 , %xmm2             # xmm2  = a4p6p1m3
   mulsd   %xmm7 , %xmm1             # (c1) *
   mulsd   %xmm9 , %xmm2             # (s1) *
   movapd  %xmm0 , %xmm3             # xmm3  = a0m2m5p7
   mulsd   %xmm7 , %xmm0             # (c1) *
   mulsd   %xmm9 , %xmm3             # (s1) *

   subsd   %xmm2 , %xmm0             # xmm0  = c(1)
   addsd   %xmm3 , %xmm1             # xmm1  = c(5)   

   movsd   %xmm0 , -24(%rsi)         # c(1)
   movsd   %xmm1 ,   8(%rsi)         # c(5)

   addq    $16, %r10                 # i0 = i0 + 2
   addq    $16, %r11                 # i1 = i1 + 2
   subq    $16, %r12                 # i2 = i2 - 2
   subq    $16, %r13                 # i3 = i3 - 2
   dec     %r14
   jnz     LOOP_I_IFFT4S

   movsd   (%r10), %xmm5             # xmm5   = a(i0)
   movsd   (%r11), %xmm8             # xmm8   = a(i1)
   movsd  8(%r10), %xmm3             # xmm3   = a(i0+1)
   movsd  8(%r11), %xmm11            # xmm11  = a(i1+1)

   movapd  %xmm5 , %xmm4             # xmm4  = a(i0)
   addsd   %xmm8 , %xmm4             # xmm4  = a(i0) + a(i1)
   subsd   %xmm8 , %xmm5             # xmm5  = a(i0) - a(i1)

   movapd  %xmm3 , %xmm9             # xmm6  = a(i0+1)
   addsd   %xmm11, %xmm9             # xmm6  = a(i0+1) + a(i1+1)
   subsd   %xmm3, %xmm11             # xmm11 = a(i1+1) - a(i0+1)

   movapd  %xmm5 , %xmm8             # xmm8  = a(i0) - a(i1)
   movsd   %xmm4 , 32(%rsi)          # c(n-3)
   subsd   %xmm9, %xmm8              # xmm8  = 

   movq  $0x3FE6A09E667F3BCD, %rdx   # rdx   = sin45 = qrt(0.5)
   movd    %rdx , %xmm4              # xmm4  = sin45
   addsd   %xmm9, %xmm5              # xmm9  = 
   mulsd   %xmm4, %xmm8              # xmm8  = c(n-2)
   mulsd   %xmm4, %xmm5              # xmm5  = c(n)
   movsd   %xmm8 , 40(%rsi)          # c(n-2)

   movq  $0x8000000000000000, %rdx   # rdx   = sign bit
   movd    %rdx , %xmm8              # xmm8  = sign bit

   movsd   %xmm11, 48(%rsi)          # c(n-1)
   xorpd   %xmm8 , %xmm5             # negate
   movsd   %xmm5 , 56(%rsi)          # c(n)

   addq    %rcx, %rdi                # advance &a
   addq    $64 , %rsi                # advance &c
   popq    %r14                      # restore loop count
   popq    %rdx                      # restore trigs
   movq    %rdx, %r8                 # r8   = trigs
   movq    %rdx, %r9                 # r9   = trigs

   dec     %r15
   jnz     LOOP_O_IFFT4S

   popq    %rbx
   popq    %r12
   popq    %r13
   popq    %r14
   popq    %r15
   movq    %rbp, %rsp
   popq    %rbp
   ret


# Fast Double Precision Fourier Transformation
# ============================================
# E. Kirk - 01-Jun-2015
# ============================================
# call ifft4m(a  ,c  ,trigs,n  ,la ,lot)
#            (rdi,rsi,rdx  ,rcx,r8 ,r9 )
#
# rdi: a(n,lot)  : real(8) source array
# rsi: c(n,lot)  : real(8) target array
# rdx: trigs(n)  : sine/cosine values
# rcx: n         : 1st. dimension
# r8 : la        : factor
# r9 : lot       : # of lines
# --------------------------------------------

   .globl _ifft4m_
_ifft4m_:

   .globl ifft4m_
ifft4m_:

   pushq   %r15
   pushq   %r14
   pushq   %r13
   pushq   %r12
   pushq   %rbp
   pushq   %rbx

   subq   $64, %rsp                  # reserve memory
   movl   (%r9), %r9d                # r9d = lot
   movq   %r9 , (%rsp)               # store lot
   movq   %rdx, %r13                 # &trigs       PERMANENT
   movl   (%rcx), %r11d              # r11 = n      PERMANENT
   movq   %r11, %rax                 # rax = n
   movq   %r11, %r15                 # r15 = n
   shrq   %r15                       # r15 = i5 = n / 2
   movl   (%r8), %r8d                # r8  = la     PERMANENT
   leaq   (,%r8,8), %r12             # r12 = la * 8 PERMANENT
   leaq   (%r8,%r8,2), %r14          # r14 = la * 3
   movq   %r8 , %rbx                 # rbx = la
   shrq   %rbx                       # rbx = la / 2
   movq   %r11, %r10                 # r10 = n
   subq   %r8 , %r10                 # r10 = i2 = n - la
   movq   %r15, %r9                  # r9  = n / 2
   subq   %r8 , %r9                  # r9  = i1 = n / 2 - la
   movq   %r9 ,  8(%rsp)             # i1
   movq   %r10, 16(%rsp)             # i2
   movq   %r15, 24(%rsp)             # i5
   movq   %r14, 32(%rsp)             # la * 3
   movq   %rdi, 40(%rsp)             # &a
   movq   %rsi, 48(%rsp)             # &c
   movq   %rbx, 56(%rsp)             # loop count la

LOOP_LOT_IFFT4M:
   movq    8(%rsp), %r9              # i1
   movq   16(%rsp), %r10             # i2
   movq   24(%rsp), %r15             # i5
   movq   32(%rsp), %r14             # la * 3
   movq   40(%rsp), %rdi             # &a
   movq   48(%rsp), %rsi             # &c
   movq   56(%rsp), %rbx             # loop count la

LOOP_1_I4FFTM:
   movapd  (%rdi)       , %xmm0      # xmm0  = a(i0)
   movapd  (%rdi,%r9 ,8), %xmm1      # xmm1  = a(i1)
   movapd  (%rdi,%r10,8), %xmm2      # xmm2  = a(i2)
   movapd  (%rdi,%r15,8), %xmm5      # xmm2  = a(i5)
   movapd  %xmm0, %xmm9              # xmm9  = a(i0)
   addpd   %xmm2, %xmm9              # xmm9  = a(i0) + a(i2)
   subpd   %xmm2, %xmm0              # xmm0  = a(i0) - a(i2)
   movapd  %xmm9, %xmm4              # xmm4  = a(i0) + a(i2)
   movapd  %xmm0, %xmm8              # xmm8  = a(i0) - a(i2)
   addpd   %xmm1, %xmm4              # xmm4  = a(i0) + a(i2) + a(i1)
   movapd  %xmm4, (%rsi)             # store   c(i)
   subpd   %xmm1, %xmm9              # xmm9  = a(i0) + a(i2) - a(i1)
   movapd  %xmm9, (%rsi,%r12,2)      # store   c(la*2+i)
   addpd   %xmm5, %xmm8              # xmm8  = a(i0) - a(i2) + a(i5)
   movapd  %xmm8, (%rsi,%r14,8)      # store   c(la*3+i)
   subpd   %xmm5, %xmm0              # xmm0  = a(i0) - a(i2) - a(i5)
   movapd  %xmm0, (%rsi,%r12)        # store   c(la+i)
   addq    $16,%rdi                  # advance &a
   addq    $16,%rsi                  # advance &c
   dec     %rbx
   jnz     LOOP_1_I4FFTM             # loop

   movq   %r11, %rbx                 # rbx = n
   leaq   (,%r8,4), %rax             # rax = la * 4
   subq   %rax, %rbx                 # rbx = i2 = n - la * 4
   movq   %rbx, %rdx                 # rdx = i2
   subq   %r15, %rdx                 # rdx = i3 = i2 - n / 2

   leaq   (%rdi,%r15,8), %r9         # r9    = &a(i1) : i1 = i0 + n / 2
   leaq   (%rdi,%rbx,8), %r10        # r10   = &a(i2) : i2 = n - 3 * la
   leaq   (%rdi,%rdx,8), %r15        # r15   = &a(i3) : i3 = n - 3 * la - n / 2
   leaq   (%r8 ,%r8 ,2), %rax        # rax   = la * 3 : rsi = &c(la)
   leaq   (%rsi,%rax,8), %rsi        # rsi   = &c(j0) : j0 = la * 4
   leaq   (%r8 ,%r8), %rbx           # rbx   = la * 2

   movq   %r11, %rax                 # rax = n
   shrq   $3, %rax                   # rax = n / 8
   movl   $0, %edx                   # rdx = 0
   divl   %r8d                       # rax = n / 8 / la
   dec    %rax                       # rax = n / 8 / la - 1
   movq   %rax, %rcx                 # store loop count
   jz     L0_I4FFTM                  # jump on zero loop

LOOP_O_I4FFTM:
   movddup  (%r13,%rbx,8), %xmm5     # xmm5  = c1 = trigs(kb)
   movddup 8(%r13,%rbx,8), %xmm6     # xmm6  = s1

   leaq     (%rbx,%rbx), %rax        # rax   = kc
   movddup  (%r13,%rax,8), %xmm7     # xmm7  = c2
   movddup 8(%r13,%rax,8), %xmm8     # xmm8  = s2

   leaq     (%rax,%rbx), %rax        # rax   = kd
   movddup  (%r13,%rax,8), %xmm9     # xmm9  = c3
   movddup 8(%r13,%rax,8), %xmm10    # xmm10 = s3

   movq   %r8 , %rax                 # rax   = loop count = la

LOOP_I_IFFT4M:
   movapd  (%rdi), %xmm1             # a(i0)
   movapd  (%r10), %xmm2             # a(i2)
   movapd  %xmm1 , %xmm15
   addpd   %xmm2 , %xmm15            # a(i0) + a(i2)
   subpd   %xmm2 , %xmm1             # a(i0) - a(i2)

   movapd  (%r9 ), %xmm11            # a(i1)
   movapd  (%r15), %xmm0             # a(i3)
   movapd  %xmm11, %xmm14
   addpd   %xmm0 , %xmm14            # a(i1) + a(i3)
   subpd   %xmm0 , %xmm11            # a(i1) - a(i3)

   movapd  (%rdi,%r12), %xmm13       # a(i4)
   movapd  (%r10,%r12), %xmm2        # a(i6)
   movapd  %xmm13, %xmm0
   addpd   %xmm2 , %xmm0             # a(i4) + a(i6)
   subpd   %xmm2 , %xmm13            # a(i4) - a(i6)

   movapd  (%r9 ,%r12), %xmm12       # a(i5)
   movapd  (%r15,%r12), %xmm3        # a(i7)
   movapd  %xmm12, %xmm2
   addpd   %xmm3 , %xmm2             # a(i5) + a(i7)
   subpd   %xmm3 , %xmm12            # a(i5) - a(i7)

   movapd  %xmm15, %xmm4             # xmm4  = a0p2
   subpd   %xmm14, %xmm4             # xmm4  = a0p2 - a1p3
   movapd  %xmm13, %xmm3             # xmm3  = a4m6
   subpd   %xmm12, %xmm3             # xmm3  = a4m6 - a5m7
   addpd   %xmm15, %xmm14            # xmm14 = a0p2 + a1p3

   movapd  %xmm14, (%rsi)            # c(j0) = a0p2 + a1p3

   addpd   %xmm13, %xmm12

   movapd  %xmm12, (%rsi,%r12,4)     # c(j4) = a4m6 - a5m7

   movapd  %xmm7, %xmm12
   mulpd   %xmm4, %xmm12
   movapd  %xmm8, %xmm13
   mulpd   %xmm3, %xmm13
   subpd   %xmm13, %xmm12

   movapd  %xmm12, (%rsi,%r12,2)     # c(j2) = c2 * a0p2m1p3 - s2 * a4m6m5m7

   mulpd   %xmm8, %xmm4
   mulpd   %xmm7, %xmm3
   addpd   %xmm4, %xmm3

   leaq    (%rsi,%r12,4), %rdx
   movapd  %xmm3, (%rdx,%r12,2)      # c(j6) = s2 * a0p2m1p3 + c2 * a4m6m5m7

   movapd  %xmm1, %xmm4
   subpd   %xmm2, %xmm4
   movapd  %xmm0, %xmm3
   addpd   %xmm11, %xmm3
   movapd  %xmm4, %xmm12
   mulpd   %xmm5, %xmm12
   movapd  %xmm3, %xmm13
   mulpd   %xmm6, %xmm13
   subpd   %xmm13, %xmm12
   movapd  %xmm12, (%rsi,%r12)       # c(j1)

   mulpd   %xmm6, %xmm4
   mulpd   %xmm5, %xmm3
   addpd   %xmm4, %xmm3

   movapd   %xmm3, (%rdx,%r12)        # c(j5)

   addpd   %xmm2, %xmm1
   subpd   %xmm11, %xmm0
   movapd  %xmm1, %xmm2
   mulpd   %xmm9, %xmm2
   movapd  %xmm0, %xmm3
   mulpd   %xmm10, %xmm3
   subpd   %xmm3, %xmm2

   leaq    (%rsi,%r12,2), %rdx
   movapd  %xmm2, (%rdx,%r12)        # c(j3)

   mulpd   %xmm10, %xmm1
   mulpd   %xmm9, %xmm0
   addpd   %xmm1, %xmm0

   leaq    (%rdx,%r12,4), %rdx
   movapd  %xmm0, (%rdx,%r12)        # c(j7)

   addq   $16, %rdi                  # &a(i0)++
   addq   $16, %r10                  # &a(i2)++
   addq   $16, %r9                   # &a(i1)++
   addq   $16, %r15                  # &a(i3)++
   addq   $16, %rsi                  # &c(j0)++
   subl   $2,%eax
   jnz   LOOP_I_IFFT4M

   leaq   (%rdi,%r12), %rdi          # adjust &a(i0)
   leaq   (%r9 ,%r12), %r9           # adjust &a(i1)
   leaq   (%r8 ,%r8 ,2), %rax        # rax = la * 3
   shlq   $3, %rax                   # rax * 8
   subq   %rax, %r10                 # adjust &a(i2)
   subq   %rax, %r15                 # adjust &a(i3)
   subq   %r12, %rsi
   leaq   (%rsi,%r12,8), %rsi        # adjust &c(j0)
   leaq   (%rbx,%r8,2), %rbx         # rbx += la * 2
   dec    %rcx                       # loop counter
   jnz   LOOP_O_I4FFTM

L0_I4FFTM:
   movq  $0x3FE6A09E667F3BCD, %rax   # rax   = sin45 = qrt(0.5)
   movd    %rax , %xmm5              # xmm5  = (sin45,   0 )
   movddup %xmm5, %xmm5              # xmm5  = (sin45,sin45)
   movq  $0x8000000000000000, %rax   # rax   = sign bit
   movd    %rax , %xmm6              # xmm6  = sign bit
   movddup %xmm6, %xmm6              # xmm6  = sign bits
   movq    %r8 , %rax                # loop count = la

LOOP_2_I4FFTM:
   movapd  (%rdi), %xmm1             # a(i0)
   movapd  (%r9 ), %xmm2             # a(i1)
   movapd  %xmm1, %xmm0
   addpd   %xmm2, %xmm0
   movapd  %xmm0, (%rsi)             # c(j0)
   subpd   %xmm2, %xmm1
   movapd  (%rdi,%r12), %xmm4        # a(i0+la)
   movapd  (%r9 ,%r12), %xmm2        # a(i1+la)
   movapd  %xmm4, %xmm0
   addpd   %xmm2, %xmm0
   movapd  %xmm1, %xmm3
   subpd   %xmm0, %xmm3
   mulpd   %xmm5, %xmm3
   movapd  %xmm3, (%rsi,%r12)        # c(j1)
   subpd   %xmm4, %xmm2
   movapd  %xmm2, (%rsi,%r12,2)      # c(j2)
   addpd   %xmm1, %xmm0
   mulpd   %xmm5, %xmm0
   xorpd   %xmm6, %xmm0
   leaq    (%rsi,%r12,2), %r14
   movapd  %xmm0, (%r14,%r12)        # c(j3)
   addq   $16, %rdi                  # &a(i0)++
   addq   $16, %r9                   # &a(i1)++
   addq   $16, %rsi                  # &c(j0)++
   subl   $2, %eax
   jnz    LOOP_2_I4FFTM

   movq   40(%rsp), %rdi             # &a
   movq   48(%rsp), %rsi             # &c
   leaq   (%rdi,%r11,8), %rdi        # &a += n * 8
   leaq   (%rsi,%r11,8), %rsi        # &c += n * 8
   movq   %rdi, 40(%rsp)             # &a
   movq   %rsi, 48(%rsp)             # &c

   decl    (%rsp)                     # --lot
   jnz    LOOP_LOT_IFFT4M

   addq   $64, %rsp                   # restore %rsp
   popq   %rbx
   popq   %rbp
   popq   %r12
   popq   %r13
   popq   %r14
   popq   %r15
   ret

# Fast Double Precision Fourier Transformation
# ============================================
# E. Kirk - 19-Mar-2015
# ============================================
# call ifft8(a  ,b  ,n  ,lot)
#           (rdi,rsi,rdx,rcx)
#
# rdi: a(n,lot)  : real(8) array
# rsi: n         : 1st. dimension
# rdx: lot       : 2nd. dimension
# ---------------------------------------


.globl _ifft8e_                      # MAC OSX names
_ifft8e_:

.globl ifft8e_                       # Linux names
ifft8e_:

    push    %rbp                     # save frame pointer
    movq    %rsp, %rbp
    call    mcount                   # enable profiling
    push    %r15
    push    %r14
    push    %r13
    push    %r12

    movl    (%rcx), %ecx             # lot
    movl    (%rdx), %edx             # n
    movl    %edx, %eax               # n
    sarl    $3, %eax                 # n / 8
    cltq
    movq    %rax, %r14               # i1 = la
    movq    %r14, %r13
    addq    %r14, %r13               # i2
    movq    %r13, %r12
    addq    %r14, %r12               # i3
    movq    %r12, %r11
    addq    %r14, %r11               # i4
    movq    %r11, %r10
    addq    %r14, %r10               # i5
    movq    %r10, %r9 
    addq    %r14, %r9                # i6
    movq    %r9 , %r8 
    addq    %r14, %r8                # i7
    movabsq $4609047870845172684, %rax
    movd    %rax, %xmm12             # SQRT(2.0)
    movlhps %xmm12, %xmm12           # expand

LOOP_O_IFFT8:
    movq    %r14, %r15               # loop count = la

LOOP_I_IFFT8:

# a0p7 = a(i0) + a(i7)
# a1m5 = a(i1) - a(i5)
#---------------------
    movupd  (%rdi)       , %xmm0     # a(i0)
    movapd  %xmm0, %xmm8             # a(i0)
    movupd  (%rdi,%r8 ,8), %xmm7     # a(i7)
    addpd    %xmm7, %xmm0            # a(i0) + a(i7)
    subpd    %xmm7, %xmm8            # a(i0) - a(i7)

# a1p5 = a(i1) + a(i5)
# a1m5 = a(i1) - a(i5)
#---------------------
    movupd  (%rdi,%r14,8), %xmm1     # a(i1)
    movapd  %xmm1, %xmm9             # a(i1)
    movupd  (%rdi,%r10,8), %xmm5     # a(i5)
    addpd   %xmm5, %xmm1             # a(i1) + a(i5)
    subpd   %xmm5, %xmm9             # a(i1) - a(i5)

# a2p6 = a(i2) + a(i6)
# a2m6 = a(i2) - a(i6)
#---------------------
    movupd  (%rdi,%r13,8), %xmm2     # a(i2)
    movapd  %xmm2, %xmm10            # a(i2)
    movupd  (%rdi,%r9 ,8), %xmm6     # a(i6)
    addpd   %xmm6, %xmm2             # a(i2) + a(i6)
    subpd   %xmm6, %xmm10            # a(i2) - a(i6)

# a0p7p3 = a0p7 + a(i3)
# a0p7m3 = a0p7 - a(i3)
#----------------------
    movupd  (%rdi,%r12,8), %xmm3     # a(i3)
    movapd  %xmm0, %xmm6             # a0p7
    movupd  (%rdi,%r11,8), %xmm4     # a(i4)
    addpd   %xmm3, %xmm0             # a0p7 + a(i3)
    subpd   %xmm3, %xmm6             # a0p7 - a(i3)

# a0m7p4   = 2.0 * (a0m7 + a(i4))
# a0m7m4   = 2.0 * (a0m7 - a(i4))
#--------------------------------
    movapd  %xmm8 , %xmm11           # a0m7
    addpd   %xmm4 , %xmm8            # a0m7 + a(i4)
    subpd   %xmm4 , %xmm11           # a0m7 - a(i4)
    addpd   %xmm8 , %xmm8            # * 2
    addpd   %xmm11, %xmm11           # * 2
        
# a1m5p2p6 = SQRT2 * (a1m5 + a2p6)
# a1m5m2p6 = SQRT2 * (a1m5 - a2p6)
#---------------------------------
    movapd  %xmm9 , %xmm14            # a1m5
    addpd   %xmm2 , %xmm9             # a1m5 + a2p6
    subpd   %xmm2 , %xmm14            # a1m5 - a2p6
    mulpd   %xmm12, %xmm9             # * SQRT(2.0)
    mulpd   %xmm12, %xmm14            # * SQRT(2.0)

# a(i0)  = 2.0 * (a0p7p3 + a1p5)
# a(i4)  = 2.0 * (a0p7p3 - a1p5)
#-------------------------------
    movapd  %xmm0, %xmm4             # a0p7p3
    addpd   %xmm1, %xmm0             # a0p7p3 + a1p5
    subpd   %xmm1, %xmm4             # a0p7p3 - a1p5
    addpd   %xmm0, %xmm0             # * 2
    addpd   %xmm4, %xmm4             # * 2
    movupd  %xmm0, (%rsi)            # store a(i0)
    movupd  %xmm4, (%rsi,%r11,8)     # store a(i4)

# a(i6)  = 2.0 * (a0p7m3 + a2m6)
# a(i2)  = 2.0 * (a0p7m3 - a2m6)
#-------------------------------
    movapd  %xmm6 , %xmm2            # a0p7m3
    addpd   %xmm10, %xmm6            # a0p7m3 + a2m6
    subpd   %xmm10, %xmm2            # a0p7m3 - a2m6
    addpd   %xmm6 , %xmm6            # * 2
    addpd   %xmm2 , %xmm2            # * 2
    movupd  %xmm6 , (%rsi,%r9 ,8)    # store a(i6)
    movupd  %xmm2 , (%rsi,%r13,8)    # store a(i2)

# a(i1)  = a0m7m4 + a1m5m2p6
# a(i5)  = a0m7m4 - a1m5m2p6
#---------------------------
    movapd  %xmm11, %xmm5            # a0m7m4
    addpd   %xmm14, %xmm11           # a0m7m4 + a1m5m2p6
    subpd   %xmm14, %xmm5            # a0m7m4 - a1m5m2p6
    movupd  %xmm11, (%rsi,%r14,8)    # store a(i1)
    movupd  %xmm5 , (%rsi,%r10,8)    # store a(i5)

# a(i3)  = a0m7p4 - a1m5p2p6
# a(i7)  = a0m7p4 + a1m5p2p6
#---------------------------
    movapd  %xmm8, %xmm3             # a0m7p4
    addpd   %xmm9, %xmm8             # a0m7p4 - a1m5p2p6
    subpd   %xmm9, %xmm3             # a0m7p4 - a1m5p2p6
    movupd  %xmm8, (%rsi,%r8 ,8)     # store a(i7)
    movupd  %xmm3, (%rsi,%r12,8)     # store a(i3)

    addq    $16, %rdi                # shift base of a
    addq    $16, %rsi                # shift base of b
    dec     %r15
    dec     %r15
    jnz     LOOP_I_IFFT8

    movq    %rdx, %rax               # n
    subq    %r14, %rax               # n - la
    leaq    (%rdi,%rax,8), %rdi      # next line
    leaq    (%rsi,%rax,8), %rsi      # next line
    dec     %ecx                     # --lot
    jnz     LOOP_O_IFFT8

    pop     %r12
    pop     %r13
    pop     %r14
    pop     %r15
    leave
    ret                              # finito


# Fast Double Precision Matrix Transposition
# ==========================================
# E. Kirk - 26-Jun-2015
# ==========================================
# call fast_mtp(a  ,n  )
#              (rdi,rsi)
#
# rdi: a(n,n)    : real(8) matrix
# rsi: n         : dimension
# ------------------------------------------

   .globl _fast_mtp
_fast_mtp:
   .globl fast_mtp_
fast_mtp_:
   .globl _fast_mtp_
_fast_mtp_:

   movl    (%rsi), %esi              # rsi = n

   movq    $8, %rax                  # (1,0) : 8
   leaq    (,%rsi,8), %r8            # (0,1) : 8 * n
   leaq    -1(%rsi), %rdx            # outer loop : n-1

LOOP_O_FAST_MTP:
   movq    %rax, %r10                # (i,j)
   movq    %r8 , %r11                # (j,i)
   movq    %rdx, %rcx                # inner loop count

LOOP_I_FAST_MTP:
   movsd   (%rdi,%r10), %xmm0        # a(i,j)
   movsd   (%rdi,%r11), %xmm1        # a(j,i)
   movsd   %xmm0, (%rdi,%r11)
   movsd   %xmm1, (%rdi,%r10)

   addq    $8, %r10                  # += 1
   leaq    (%r11,%rsi,8), %r11       # += n

   dec     %rcx                      # inner loop counter
   jnz     LOOP_I_FAST_MTP

   leaq   8(%rax,%rsi,8), %rax       # rax += n + 1
   leaq   8(%r8 ,%rsi,8), %r8        # r8  += n + 1
   
   dec     %rdx                      # outer loop counter
   jnz     LOOP_O_FAST_MTP

   ret


# Fast Double Precision Fourier Transposition
# ==========================================
# E. Kirk - 25-Jun-2015
# ==========================================
# call fast_ftp(a  ,b  ,n  )
#              (rdi,rsi,rdx)
#
# rdi: a(*,*)    : complex fourier coefficients
# rsi: b(n,n)    : square matrix
# rdx: n         : dimension
# ------------------------------------------
# rax: multi purpose
# rbx: multi purpose, loop length
# rcx: loop counter
# rdx: 3rd. parameter, division rest
# rsp: only chnaged by push / pop
# rbp: unused
# rdi: array pointer a
# rsi: array pointer b
# r8 : byte count for one row of matrix a
# r9 : index register for a(i,k)
# r10: column address for a(k,i)
# r11: n/2 - k (# of fill positions at end of row)
# r12: n/2
# r13: index for accessing a(k,i)
# r14: byte count for one row of matrix b (n * 8)
# r15: loop index rows
# ------------------------------------------
# xmm0 : 0.0
# xmm5 : 0.5 , 0.5

# source array example for T2
#             0           1           2
# -------------------------------------
#  0   R00    0    R10  I10    R20  I20
#  1   R01  I01    R11  I11    R21  I21
#  2   R02  I02    R12  I12    R22  I22
# -2   P02  Q02    P12  Q12    P22  Q22
# -1   P01  Q01    P11  Q11    P21  Q21
#
# target array example for N8
#        0    1      2    3      4    5      6    7
# -------------------------------------------------
#  0   X00    0    X10  Y10    X20  Y20      0    0
#  1     0    0      0    0      0    0      0    0
#  2   X02    0    X12  Y12    X22  Y22      0    0
#  3   X03    0    X13  Y13    X23  Y23      0    0
#  4   X04    0    X14  Y14    X24  Y24      0    0
#  5   X05    0    X15  Y15    X25  Y25      0    0
#  6     0    0      0    0      0    0      0    0
#  7     0    0      0    0      0    0      0    0

# (X00,  0) = (R00,  0)
# (X10,Y10) = (R01,I01)   (X20,Y20) = (R02,I02)
# (X02,  0) = (R10,  0)   (X03,  0) = (I10,  0)   ...
# (X12,Y12) = 0.5 * (R11,I11) +- (P11,Q11))
# (X13,Y13) = 0.5 * (P11,Q11) -+ (R11,I11))
# (X22,Y22) = 0.5 * (R12,I12) +- (P12,Q12))
# (X23,Y23) = 0.5 * (P12,Q12) -+ (R12,I12))


   .globl _fast_ftp
_fast_ftp:
   .globl fast_ftp_
fast_ftp_:
   .globl _fast_ftp_
_fast_ftp_:


   push    %rbx
   push    %r12
   push    %r13
   push    %r14
   push    %r15

   xorpd   %xmm0, %xmm0              # xmm0 = 0
   movq    $0x3fe0000000000000, %rax # rax  = 0.5
   movd    %rax, %xmm5               # xmm5 = 0.5
   movddup %xmm5, %xmm5              # xmm5 = 0.5 , 0.5
   movl    (%rdx), %eax              # rax  = n
   movq    %rax, %r11                # r11  = n
   movq    %rax, %r14                # r14  = n
   shl     $3, %r14                  # r14  = n * 8 (byte count)
   movl    $3, %ebx                  # rbx  = 3
   xorq    %rdx, %rdx                # rdx  = 0
   divl    %ebx                      # rax  = n/3
   inc     %rax                      # rax  = k = (n/3) + 1
   sarq    %r11                      # r11  = n/2
   subq    %rax, %r11                # r11  = n/2 - k
   movq    %rax, %r8                 # r8   = k
   movq    %rax, %rbx                # rbx  = k
   dec     %rbx                      # wavenumber excluding 0
   shlq    $4, %r8                   # r8  = row increment [bytes]

# the imaginary part of mode (0,0) is zero by definition

   movsd   (%rdi), %xmm1             # xmm1 = a(0,0)
   movapd  %xmm1, (%rsi)             # b(0,0) = xmm1
   movapd  %xmm0, (%rsi,%r14)        # b(0,1) = 0.0
   mov     %r8, %r12                 # r12 = row
   addq    $16, %rsi                 # &b += 16 bytes (one complex)
   movq    %rax, %rcx                # rcx = k
   dec     %rcx

# transpose row 0 and set row 1 to zero

LOOP_1_FAST_FTP:
   movapd  (%rdi,%r12), %xmm1        # xmm1 = a(0,i)
   movapd  %xmm1, (%rsi)             # b(i,0) = xmm1
   movapd  %xmm0, (%rsi,%r14)        # b(i,1) = 0.0
   addq    %r8, %r12                 # r12+= row
   addq    $16, %rsi                 # &b += 16 bytes (one complex)
   dec     %rcx
   jnz     LOOP_1_FAST_FTP

# clear rest of rows

   movq    %r11, %rcx                # rcx = n/2 - k

LOOP_2_FAST_FTP:
   movapd  %xmm0, (%rsi)             # b(i,0) = 0.0
   movapd  %xmm0, (%rsi,%r14)        # b(i,1) = 0.0
   addq    $16, %rsi                 # &b += 16 bytes (one complex)
   dec     %rcx
   jnz     LOOP_2_FAST_FTP
   addq    %r14, %rsi                # advance one row

# compute column 0

   movq    $16, %r10                 # r10 = index j
   movq    %r8, %r13                 # r13 = rowsize
   sal     $1, %rax                  # rax = k * 2
   subl    $2, %eax                  # rax = k * 2 - 3
   mul     %r13d                     # rax = rowsize * (k-1)
   leaq    (%r10, %rax), %r9         # r9  = index k

   movq    %rbx, %r15                # loop count rows

LOOP_O_FAST_FTP:
   movq    %r10, %r12                # j
   movq    %r9 , %r13                # k
   movsd   (%rdi,%r12), %xmm1        # xmm1 = a(i  ,0)
   movsd   8(%rdi,%r12), %xmm2       # xmm1 = a(i+1,0)
   movapd  %xmm1, (%rsi)             # b(0,i  ) = xmm1
   movapd  %xmm2, (%rsi,%r14)        # b(0,i+1) = xmm2
   addq    $16, %rsi                 # &b += 16 bytes (one complex)

# compute fourier elemnt for wavenumber i,j (4 values)

   movq    %rbx, %rcx                # loop count = k
   addq    %r8, %r12                 # r12+= row
   

LOOP_I_FAST_FTP:
   movapd  (%rdi,%r12), %xmm1        # xmm1 = a(j,i)
   movapd  %xmm1, %xmm7              # xmm7 = a(j,i)
   movapd  (%rdi,%r13), %xmm6        # xmm6 = a(j,k)
   xorpd   %xmm2, %xmm2              # xmm2 = 0
   subpd   %xmm6, %xmm2              # xmm2 = -a(j,k)
   addsubpd %xmm2, %xmm1             # xmm1 = a(j,i) -+ a(j,k)
   mulpd   %xmm5, %xmm1              # xmm1 * 0.5
   movapd  %xmm1, (%rsi)             # b(i,j)
   addsubpd %xmm7, %xmm6             # xmm6 = a(j,k) -+ a(j,i)
   mulpd   %xmm5, %xmm6              # xmm6 * 0.5
   shufpd  $1, %xmm6, %xmm6          # swap cat words
   movapd  %xmm6, (%rsi,%r14)        # b(i,j+n)

   addq    %r8, %r12                 # r12 += row
   subq    %r8, %r13                 # r13 -= row
   addq    $16, %rsi                 # &b += 16 bytes (one complex)
   dec     %rcx
   jnz     LOOP_I_FAST_FTP

# fill rest of row

   movq    %r11, %rcx                # rcx = n/2 - k

LOOP_3_FAST_FTP:
   movapd  %xmm0, (%rsi)             # b(i,j  ) = 0
   movapd  %xmm0, (%rsi,%r14)        # b(i,j+1) = 0
   addq    $16, %rsi                 # &b += 16 bytes (one complex)
   dec     %rcx
   jnz     LOOP_3_FAST_FTP

   addq    %r14,  %rsi               # advance row
   addq    $16, %r10                 # advance column
   addq    $16, %r9

   dec     %r15
   jnz     LOOP_O_FAST_FTP

# fill rest of target array

   movq    %r11, %rax                # rax = n/2 - k
   mul     %r14                      # rax = (n*8) * (n/2 - k)
   sar     $3, %rax                  # rax = n * (n/2 - k)
   
LOOP_4_FAST_FTP:
   movapd  %xmm0, (%rsi)             # b(i,j  ) = 0
   addq    $16, %rsi                 # &b += 16 bytes (one complex)
   dec     %rax
   jnz     LOOP_4_FAST_FTP

   pop     %r15
   pop     %r14
   pop     %r13
   pop     %r12
   pop     %rbx

   ret

# Fast Double Precision Fourier Transformation
# ============================================
# E. Kirk - 08-Jul-2015
# ============================================
# call dfft2(a  ,c  ,trg,n  ,lot)
#           (rdi,rsi,rdx,rcx,r8 )
#
# rdi: a(n,lot)  : real(8) source
# rsi: c(n,lot)  : real(8) traget
# rdx: trg(n)    : real(8) trigs
# rcx: n         : 1st. dimension
# r8 : lot       : 2nd. dimension
# --------------------------------------------

   .globl _dfft2e_
_dfft2e_:

   .globl dfft2e_
dfft2e_:

   pushq   %rbx
   pushq   %r12
   pushq   %r13
   pushq   %r14

   movq    %rdx, %r13                # r13   = &trigs
   movl    (%r8 ), %r14d             # r14   = lot
   movl    (%rcx), %r10d             # r10   = n
   movq    %r10, %rbx                # rbx   = n
   shrq    $2, %rbx                  # rbx   = n / 4
   dec     %rbx                      # rbx   = (n / 4) - 1
   movq    %rbx, %r12                # r12   = loop count
   movq    $0x8000000000000000, %rax # rax   = sign bit
   movd    %rax, %xmm7               # xmm7  = sign bit
   shufpd  $1, %xmm7, %xmm7          # swap cat words

LOOP_O_DFFT2:
   movq    $0, 8(%rsi)               # c(1)  = 0
   movsd   (%rdi), %xmm0             # xmm0  = a(0)
   addsd   8(%rdi), %xmm0            # xmm0  = a(0) + a(1)
   movsd   %xmm0, (%rsi)             # c(0)
   leaq     16(%rsi), %r9            # r9    = &c(2)
   leaq    -16(%rsi,%r10,8), %r8     # r8    = &c(n-2)
   addq    $16, %rdi                 # rdi   = &a(2)

LOOP_I_DFFT2:
   addq    $16, %rdx                 # advance &trigs
   movsd    (%rdx), %xmm5            # xmm5  = cos = trigs(ja)
   movsd   8(%rdx), %xmm4            # xmm4  = sin = trigs(ja+1)
   movsd   8(%rdi), %xmm3            # xmm3  = a(i+1)
   movapd  %xmm5, %xmm2              # xmm2  = cos
   movsd   24(%rdi), %xmm0           # xmm0  = a(i+3)
   movapd  %xmm4, %xmm1              # xmm1  = sin
   mulsd   %xmm3, %xmm2              # xmm2  = cos * a(i+1) 
   mulsd   %xmm0, %xmm1              # xmm1  = sin * a(i+3)
   mulsd   %xmm4, %xmm3              # xmm3  = sin * a(i+1)
   mulsd   %xmm5, %xmm0              # xmm0  = cos * a(i+3)
   addsd   %xmm1, %xmm2              # xmm2  = a1p3
   movsd   (%rdi), %xmm1             # a(i)
   subsd   %xmm3, %xmm0              # xmm0  = a3m1
   movapd  %xmm2, %xmm3              # xmm3  = a1p3
   addsd   %xmm1, %xmm3              # xmm3  = a(i) + a1p3
   subsd   %xmm2, %xmm1              # xmm1  = a(i) - a1p3
   movapd  %xmm0, %xmm2              # xmm2  = a3m1
   movsd   %xmm3, (%r9 )             # c(ja)
   movsd   %xmm1, (%r8)              # c(jb)
   movsd   16(%rdi), %xmm1           # a(i+2)
   addsd   %xmm1, %xmm2              # xmm1  = a3m1 + a(i+2)
   subsd   %xmm1, %xmm0              # xmm0  = a3m1 - a(i+2)
   movsd   %xmm2, 8(%r9 )            # c(ja+1)
   movsd   %xmm0, 8(%r8)             # c(jb+1)
   addq    $32, %rdi                 # &a(i ) += 4
   addq    $16, %r9                  # &c(ja) += 2
   subq    $16, %r8                  # &c(jb) -= 2
   dec    %r12
   jnz    LOOP_I_DFFT2

   movapd  (%rdi), %xmm0             # a(n-2), a(n-1)
   xorpd   %xmm7, %xmm0              # negate  a(n-1)
   movapd  %xmm0, (%r9 )             # c(ja), c(ja+1)

   addq    $16, %rdi                 # advance &a += 2
   leaq    (%rsi,%r10,8), %rsi       # advance &c += n
   movq    %r13, %rdx                # restore trigs
   movq    %rbx, %r12                # inner loop count
   dec     %r14                      # lot--
   jnz      LOOP_O_DFFT2

   popq    %r14
   popq    %r13
   popq    %r12
   popq    %rbx
   ret


# Fast Double Precision Fourier Transformation
# ============================================
# E. Kirk - 27-Aug-2015
# ============================================
# call dfft8(a  ,c  ,n  ,lot)
#           (rdi,rsi,rdx,rcx)
#
# rdi: a(n,lot)  : real(8) array
# rsi: c(n,lot)  : real(8) array
# rdx: n         : 1st. dimension
# rcx: lot       : 2nd. dimension
# --------------------------------------------

   .globl _dfft8s_
_dfft8s_:

   .globl dfft8s_
dfft8s_:

   pushq   %r15
   pushq   %r14
   pushq   %r13
   pushq   %r12
   pushq   %rbp
   pushq   %rbx

   movq    $0x3ff0000000000000, %rax # rax   = 1.0
   movd    %rax  , %xmm14            # xmm14 = 1.0
   movddup %xmm14, %xmm14            # xmm14 = 1.0 , 1.0
   movq    $0x3FE6A09E667F3BCD, %rax # rax   = sqrt(0.5)
   movd    %rax  , %xmm15            # xmm15 = sqrt(0.5)
   movddup %xmm15, %xmm15            # xmm15 = sqrt(0.5) , sqrt(0.5)

   movl   (%rdx), %edx               # rdx   = n
   movl   (%rcx), %ecx               # rcx   = lot
   cvtsi2sd  %edx, %xmm13            # xmm13 = real(n)
   movddup %xmm13, %xmm13            # xmm13 = real(n) , real(n)
   movq    %rdx, %rbp                # rbp   = n
   sarq    $3  , %rbp                # rbp   = la = n/8
   divpd   %xmm13, %xmm14            # xmm14 = 1.0 / n
   mulpd   %xmm14, %xmm15            # xmm15 = sqrt(0.5) / n
   leaq    (    ,%rbp,8), %r11       # r11   = i1
   leaq    (    ,%r11,2), %r12       # r12   = i2
   leaq    (%r11,%r11,2), %r13       # r13   = i3
   leaq    (    ,%r11,4), %r14       # r14   = i4
   leaq    (%r11,%r11,4), %r15       # r15   = i5
   leaq    (    ,%r13,2), %r8        # r8    = i6
   leaq    (%r8 ,%r11  ), %r9        # r9    = i7
   subq    %rbp, %rdx                # rdx   = n - la
   sarq    %rbp                      # rbp   = la / 2

LOOP_O_DFFT8:
   movq    %rbp, %r10                # r10   = loop count = la / 2

LOOP_I_DFFT8:

   movapd  (%rdi     ), %xmm0        # xmm0 = a(i0)
   movapd  (%rdi,%r11), %xmm1        # xmm1 = a(i1)
   movapd  (%rdi,%r12), %xmm2        # xmm2 = a(i2)
   movapd  (%rdi,%r13), %xmm3        # xmm3 = a(i3)
   movapd  (%rdi,%r14), %xmm4        # xmm4 = a(i4)
   movapd  (%rdi,%r15), %xmm5        # xmm5 = a(i5)
   movapd  (%rdi,%r8 ), %xmm6        # xmm6 = a(i6)
   movapd  (%rdi,%r9 ), %xmm7        # xmm7 = a(i7)

   movapd   %xmm0 , %xmm8            #  a(0)
   addpd    %xmm4 , %xmm8            #  a(0) + a(4)
   subpd    %xmm4 , %xmm0            #  a(0) - a(4)      xmm4  free
   mulpd    %xmm14, %xmm0            # (a(0) - a(4)) / n

   movapd   %xmm6 , %xmm9            #  a(6)
   addpd    %xmm2 , %xmm9            #  a(6) + a(2)
   subpd    %xmm2 , %xmm6            #  a(6) - a(2)      xmm2  free
   mulpd    %xmm14, %xmm6            # (a(6) - a(2)) / n

   movapd   %xmm5 , %xmm10           #  a(5)
   addpd    %xmm1 , %xmm10           #  a(5) + a(1)
   subpd    %xmm1 , %xmm5            #  a(5) - a(1)      xmm1  free

   movapd   %xmm7 , %xmm4            #  a(7)             xmm4  reuse
   addpd    %xmm3 , %xmm4            #  a(7) + a(3)
   subpd    %xmm3 , %xmm7            #  a(7) - a(3)      xmm3  free

   movapd   %xmm8 , %xmm1            #  a0+a4            xmm1  reuse
   addpd    %xmm9 , %xmm1            #  a0+a4 + a6+a2
   subpd    %xmm9 , %xmm8            #  a0+a4 - a6+a2    xmm9  free
   mulpd    %xmm14, %xmm8            # (a0+a4 - a6+a2) / n

   movapd   %xmm4 , %xmm2            #  a7+a3            xmm2  reuse
   addpd    %xmm10, %xmm2            #  a7+a3 + a5+a1
   subpd    %xmm10, %xmm4            #  a7+a3 - a5+a1    xmm10 free
   mulpd    %xmm14, %xmm4            # (a7+a3 - a5+a1) / n

   movapd   %xmm7 , %xmm3            #  a7-a3            xmm3  reuse
   addpd    %xmm5 , %xmm3            #  a7-a3 + a5-a1
   subpd    %xmm5 , %xmm7            #  a7-a3 - a5-a1    xmm5  free
   mulpd    %xmm15, %xmm3            # (a7-a3 + a5-a1) * sin(45) / n
   mulpd    %xmm15, %xmm7            # (a7-a3 - a5-a1) * sin(45) / n

   movapd   %xmm0 , %xmm5            # a0-a4                xmm5  reuse
   addpd    %xmm7 , %xmm5            # a0-a4 + a7-a3-a5-a1
   subpd    %xmm7 , %xmm0            # a0-a4 - a7-a3-a5-a1  xmm7  free

   movapd   %xmm1 , %xmm7            # a0+a4+a6+a2               xmm7  reuse
   addpd    %xmm2 , %xmm7            # a0+a4+a6+a2 + a7+a3+a5+a
   subpd    %xmm2 , %xmm1            # a0+a4+a6+a2 - a7+a3+a5+a  xmm2  free
   mulpd    %xmm14, %xmm7            # xmm7 / n
   mulpd    %xmm14, %xmm1            # xmm1 / n

   movapd   %xmm3 , %xmm2            # a7-a3+a5-a1          xmm2  reuse
   addpd    %xmm6 , %xmm2            # a7-a3+a5-a1 + a6m2
   subpd    %xmm6 , %xmm3            # a7-a3+a5-a1 - a6m2   xmm6  free

   movapd   %xmm7 , (%rsi)           # c(i0)
   movapd   %xmm5 , (%rsi,%r11)      # c(i1)
   movapd   %xmm2 , (%rsi,%r12)      # c(i2)
   movapd   %xmm8 , (%rsi,%r13)      # c(i3)
   movapd   %xmm4 , (%rsi,%r14)      # c(i4)
   movapd   %xmm0 , (%rsi,%r15)      # c(i5)
   movapd   %xmm3 , (%rsi,%r8 )      # c(i6)
   movapd   %xmm1 , (%rsi,%r9)       # c(i7)

   addq    $16, %rdi                 # shift base of a
   addq    $16, %rsi                 # shift base of c
   dec    %r10
   jnz   LOOP_I_DFFT8

   leaq   (%rdi,%rdx,8), %rdi        # &a += n - la
   leaq   (%rsi,%rdx,8), %rsi        # &c += n - la
   dec    %rcx                       # lot--
   jne   LOOP_O_DFFT8

   popq   %rbx
   popq   %rbp
   popq   %r12
   popq   %r13
   popq   %r14
   popq   %r15
   ret


# Fast Double Precision Fourier Transformation
# ============================================
# E. Kirk - 14-Jul-2015
# ============================================
# call dfft4(a  ,c  ,trg,n  ,lot)
#           (rdi,rsi,rdx,rcx,r8 )
#
# rdi: a(n,lot)  : real(8) source
# rsi: c(n,lot)  : real(8) traget
# rdx: trg(n)    : real(8) trigs
# rcx: n         : 1st. dimension
# r8 : lot       : 2nd. dimension
# --------------------------------------------

   .globl _dfft4e_
_dfft4e_:

   .globl dfft4e_
dfft4e_:

   pushq   %rbp
   push    %rbx

   movq    $0x3FE6A09E667F3BCD, %rbx # rbx   = sqrt(0.5)
   movl    (%r8) , %r8d              # r8    = lot
   movl    (%rcx), %ecx              # rcx   = n
   leaq    -8(%rcx), %rbp            # rbp   = n-8
   sarq    $3, %rbp                  # rbp   = n/8 - 1

LOOP_O_DFFT4E:
   movsd     (%rdi), %xmm2           # xmm2  = a(0)
   movsd    8(%rdi), %xmm3           # xmm3  = a(1)
   movsd   16(%rdi), %xmm4           # xmm4  = a(2)
   movsd   24(%rdi), %xmm1           # xmm1  = a(3)

   movapd  %xmm2, %xmm0              # xmm0  = a(0)
   addsd   %xmm3, %xmm0              # xmm0  = a(0) + a(1)
   subsd   %xmm4, %xmm2              # xmm2  = a(0) - a(2)
   addsd   %xmm4, %xmm0              # xmm0  = a(0) + a(1) + a(2)
   addsd   %xmm1, %xmm0              # xmm0  = a(0) + a(1) + a(2) + a(3)
   subsd   %xmm3, %xmm1              # xmm1  = a(3) - a(1)

   movsd   %xmm0,  (%rsi)            # c(0)
   movq    $0   , 8(%rsi)            # c(1)
   movsd   %xmm2,  (%rsi,%rcx,4)     # c(n/2)
   movsd   %xmm1, 8(%rsi,%rcx,4)     # c(n/2+1)

   leaq   -16(%rsi,%rcx,4), %r10     # r10   = &c(n/2-2)  j3
   addq   $32, %rdi                  # &rdi += 4
   addq   $16, %rsi                  # &rsi += 2
   movq   $4, %r9                    # r9    = 4
   movq   %rbp, %r11                 # r11   = n/8 - 1 loop count

LOOP_I_DFFT4E:
   leaq       (%rdx,%r9,8), %rax     # rax   = &trg(kc)
   movsd      (%rax,%r9,4), %xmm1    # xmm1  = c3 = trg(kd  )
   movsd     8(%rax,%r9,4), %xmm9    # xmm9  = s3 = trg(kd+1)

   movlpd     (%rdx,%r9,4), %xmm15   # xmm15 = c1 = trg(kb  )
   movhpd     (%rdx,%r9,8), %xmm15   # xmm15 = c2 = trg(kc  )
   movlpd    8(%rdx,%r9,4), %xmm14   # xmm14 = s1 = trg(kb+1)
   movhpd    8(%rdx,%r9,8), %xmm14   # xmm12 = s2 = trg(kc+1)

   movapd   %xmm15, %xmm10           # xmm10 = c1 , c2
   movapd   %xmm14, %xmm6            # xmm6  = s1 , s2

   movupd    8(%rdi), %xmm0          # xmm0  = a(1) , a(2)
   mulpd    %xmm0 , %xmm10           # xmm10 = c1 * a(1) , c2 * a(2)
   mulpd    %xmm14, %xmm0            # xmm0  = s1 * a(1) , s2 * a(2)

   movupd   40(%rdi), %xmm3          # xmm3  = a(5),a(6)
   mulpd    %xmm3 , %xmm6            # xmm6  = s1 * a(5) , s2 * a(6)
   mulpd    %xmm15, %xmm3            # xmm3  = c1 * a(5) , c2 * a(6)

   addpd    %xmm6 , %xmm10           # xmm10 = a1p5 , a2p6
   subpd    %xmm0 , %xmm3            # xmm3  = a5m1 , a6m2

   movsd    24(%rdi), %xmm2          # xmm2  = a(3)
   movsd    56(%rdi), %xmm7          # xmm7  = a(7)
   movapd   %xmm1 , %xmm8            # xmm8  = c3
   mulsd    %xmm7 , %xmm1            # xmm1  = c3 * a(7)
   mulsd    %xmm2 , %xmm8            # xmm8  = c3 * a(3)
   mulsd    %xmm9 , %xmm2            # xmm2  = s3 * a(3)
   movapd   %xmm9 , %xmm6            # xmm6  = s3
   mulsd    %xmm7 , %xmm6            # xmm6  = s3 * a(7)

   movapd   %xmm10, %xmm7            # xmm7  = a1p5 , a2p6
   subsd    %xmm2 , %xmm1            # xmm1  = a7m3
   movapd   %xmm3 , %xmm4            # xmm4  = a5m1 , a6m2
   movhpd   32(%rdi), %xmm1          # xmm1  = a7m3 , a(4)
   addsd    %xmm6 , %xmm8            # xmm8  = a3p7
   movhpd   (%rdi), %xmm8            # xmm8  = a3p7 , a(0)
   subpd    %xmm1 , %xmm3            # xmm3  = b3 , b2

   addpd    %xmm8 , %xmm7            # xmm7  = a1 , a0
   shufpd   $1, %xmm7, %xmm7         # xmm7  = a0 , a1
   movapd   %xmm7 , %xmm9            # xmm9  = a0 , a1
   addpd    %xmm1 , %xmm4            # xmm4  = b1 , b0
   haddpd   %xmm4 , %xmm7            # xmm7  = a0 + a1 , b1 + b0
   hsubpd   %xmm4 , %xmm9            # xmm9  = a0 - a1 , b1 - b0
   movapd   %xmm7 , (%rsi)           # c(j0)
   movapd   %xmm9 , (%r10,%rcx,4)    # c(j2)

   subpd    %xmm10, %xmm8            # xmm8  = a3 , a2
   shufpd   $1, %xmm8, %xmm8         # xmm8  = a2 , a3
   movapd   %xmm8 , %xmm4            # xmm4  = a2 , a3
   addpd    %xmm3 , %xmm8            # xmm8  = a2 + b3 , a3 + b2
   subpd    %xmm3 , %xmm4            # xmm4  = a2 - b3 , a3 - b2

   movsd    %xmm8 , (%rsi,%rcx,4)    # c(j1)
   movsd    %xmm4 , (%r10)           # c(j3)
   movhpd   %xmm4 , 8(%rsi,%rcx,4)   # c(j1+1)
   movhpd   %xmm8 , 8(%r10)          # c(j3+1)

   addq    $64, %rdi                 # rdi  += 8  &a(i)
   addq    $16, %rsi                 # rsi  += 2  &c(j0)
   subq    $16, %r10                 # r10  -= 2  &c(j3)
   addq    $4 , %r9
   dec    %r11
   jnz   LOOP_I_DFFT4E

   movsd     (%rdi), %xmm4           # xmm4  = a(0)
   movsd    8(%rdi), %xmm0           # xmm0  = a(1)
   movsd   16(%rdi), %xmm5           # xmm5  = a(2)
   movsd   24(%rdi), %xmm3           # xmm3  = a(3)
   movd    %rbx  , %xmm1             # xmm1  = sqrt(0.5)
   movapd   %xmm0, %xmm2             # xmm2  = a(1)
   subsd    %xmm3, %xmm0             # xmm0  = a(1) - a(3)
   addsd    %xmm3, %xmm2             # xmm2  = a(2) + a(3)
   mulsd    %xmm1, %xmm0             # xmm0  = sin45 * (a(1)-a(3))
   mulsd    %xmm1, %xmm2             # xmm2  = sin45 * (a(1)+a(3))

   movapd   %xmm0, %xmm3             # xmm3  = a1m3
   addsd    %xmm4, %xmm3             # xmm3  = a(0) + a1m3
   subsd    %xmm0, %xmm4             # xmm4  = a(0) - a1m3
   xorpd    %xmm0, %xmm0             # xmm0  = 0.0
   subsd    %xmm5, %xmm0             # xmm0  = -a(2) 
   subsd    %xmm2, %xmm5             # xmm5  =  a(2) - a1p3
   subsd    %xmm2, %xmm0             # xmm0  = -a(2) - a1p3

   movsd    %xmm3, (%rsi)            # c(j0)
   movsd    %xmm4, (%rsi,%rcx,4)     # c(j1)
   movsd    %xmm0, 8(%rsi)           # c(j0+1)
   movsd    %xmm5, 8(%rsi,%rcx,4)    # c(j1+1)

   addq     $32, %rdi                # &rdi += 4
   leaq     (%rsi,%rcx,4), %rsi      # &rsi += n/2
   leaq     (%rsi,%rcx,2), %rsi      # &rsi += n/4

   dec      %r8
   jnz      LOOP_O_DFFT4E

   popq     %rbx
   popq     %rbp
   ret


# Fast Double Precision Fourier Transformation
# ============================================
# E. Kirk - 27-Aug-2015
# ============================================
# call dfft4m(a  ,c  ,trg,n  ,lot,la)
#            (rdi,rsi,rdx,rcx,r8 ,r9)
#
# rdi: a(n,lot)  : real(8) source
# rsi: c(n,lot)  : real(8) traget
# rdx: trg(n)    : real(8) trigs
# rcx: n         : 1st. dimension
# r8 : lot       : # of transforms
# r9 : la        : current factor
# --------------------------------------------

   .globl _dfft4m_
_dfft4m_:

   .globl dfft4m_
dfft4m_:

    pushq    %r15
    pushq    %r14
    pushq    %r13
    pushq    %r12
    pushq    %rbp
    pushq    %rbx

    movl    (%r8 ), %r15d            # r15   = lot
    movl    (%rcx), %r11d            # r11   = n

    movl    (%r9), %ebx              # rbx   = la
    leaq    (,%rbx,8), %rbp          # rbp   = la*8
    leaq    (%rbp,%rbp,2), %r13      # r13   = la*8*3
    leaq    (%r13,%rbp,4), %r12      # r12   = la*8*7

LOOP_DFFT4M_LOT:

    leaq    (,%r11,4), %r9           # r9    = n/2*8
    subq    %rbp, %r9                # r9    = n/2*8-la*8 : j1
    leaq    (,%r11,8), %r10          # r10   = n*8
    subq    %rbp, %r10               # r10   = n*8-la*8   : j2

    movq    %rbx, %r14               # r14   = la : loop index

LOOP_DFFT4M_A:
    movapd  (%rdi       ), %xmm0     # xmm0  = a(i     )
    movapd  (%rdi,%rbp  ), %xmm1     # xmm1  = a(i+  la)
    movapd  (%rdi,%rbp,2), %xmm2     # xmm2  = a(i+2*la)
    movapd  (%rdi,%r13  ), %xmm3     # xmm3  = a(i+3*la)

    movapd  %xmm0, %xmm4             # xmm4  = a(i)
    addpd   %xmm2, %xmm4             # xmm4  = a0p2
    movapd  %xmm1, %xmm5             # xmm5  = a(i+la)
    subpd   %xmm2, %xmm0             # xmm0  = a(i) - a(i+2*la)
    addpd   %xmm3, %xmm5             # xmm5  = a1p3
    movapd  %xmm4, %xmm6             # xmm6  = a0p2
    subpd   %xmm1, %xmm3             # xmm3  = a(i+3*la) - a(i+la)
    addpd   %xmm5, %xmm6             # xmm6  = a0p2 + a1p3
    subpd   %xmm5, %xmm4             # xmm4  = a0p2 - a1p3

    movapd  %xmm6, (%rsi       )     # c(i   )
    movapd  %xmm0, (%rsi,%r9   )     # c(i+j1)   j1 = n/2-la
    movapd  %xmm4, (%rsi,%r10  )     # c(i+j2)   j2 = n-la
    movapd  %xmm3, (%rsi,%r11,4)     # c(i+j5)   j5 = n/2

    addq    $16, %rdi                # &a += 2
    addq    $16, %rsi                # &c += 2
    subq    $2 , %r14                # loopcount -= 2
    jnz     LOOP_DFFT4M_A


    addq    %r13, %rdi               # rdi   = &a(la*4)
    leaq    (    ,%rbp,4), %rax      # rax   = la*8*4
    leaq    (%rsi,%r11,4), %r8       # r8    = &c(j0+n/2)
    subq     %rax,%r8                # r8    = &c(j3)
    leaq    (%r8 ,%r11,4), %r10      # r10   = &c(j2)
    leaq    (%rsi,%r11,4), %r9       # r9    = &c(j1) 

    leaq    (%rbp,%rbp,4), %rax      # rax   = la*8*5
    movq    %rbp, %r14               # r14   = la*8 (outer index)


LOOP_DFFT4M_O:

    leaq      (%r14,%r14,2), %rcx    # rcx   = r14*3
    movddup   (%rdx,%r14,2), %xmm10  # xmm10 = trg(kb)    : c1
    movddup  8(%rdx,%r14,2), %xmm11  # xmm11 = trg(kb+1)  : s1
    movddup   (%rdx,%r14,4), %xmm12  # xmm12 = trg(kc)    : c2
    movddup  8(%rdx,%r14,4), %xmm13  # xmm13 = trg(kc+1)  : s2
    movddup   (%rdx,%rcx,2), %xmm14  # xmm14 = trg(kd)    : c3
    movddup  8(%rdx,%rcx,2), %xmm15  # xmm15 = trg(kd+1)  : s3

    movq    %rbx, %rcx               # rcx   = loop count = la

# registers used in inner loop
# ==================================
# %rbx : constant   : la
# %rbp : constant   : la*8
# %r11 : constant   : n
# %rdx : constant   : &trg
#-----------------------------------
# %r15 : outer loop : lot
# %rcx : inner loop : la/2
# %r14 : middl loop : 
#-----------------------------------
# %rdi : base a(i0) : += 16
# %rsi : base c(j0) : += 16
#-----------------------------------
# %rbp : la*8   : index for i1,i2,i4
# %r13 : la*8*3 : index for i3,i6
# %r12 : la*8*7 : index for i7
#-----------------------------------
# %r8  :        : index for j3
# %r9  : n/2    : index for j1
# %r10 : n-la*4 : index for j2


LOOP_DFFT4M_I:

    movapd   (%rdi,%rbp  ), %xmm1    # xmm1  = a(i1)
    movapd   (%rdi,%rbp,2), %xmm2    # xmm2  = a(i2)
    movapd   (%rdi,%r13  ), %xmm3    # xmm3  = a(i3)
    movapd   (%rdi,%rax  ), %xmm5    # xmm5  = a(i5)
    movapd   (%rdi,%r13,2), %xmm6    # xmm6  = a(i6)
    movapd   (%rdi,%r12  ), %xmm7    # xmm7  = a(i7)

    movapd    %xmm10, %xmm9          # xmm9  = c1
    movapd    %xmm11, %xmm0          # xmm0  = s1
    movapd    %xmm12, %xmm8          # xmm8  = c2
    movapd    %xmm14, %xmm4          # xmm4  = c3

    mulpd     %xmm1 , %xmm9          # xmm9  = c1 * a(i1)
    mulpd     %xmm5 , %xmm0          # xmm0  = s1 * a(i5)
    mulpd     %xmm2 , %xmm8          # xmm8  = c2 * a(i2)
    mulpd     %xmm13, %xmm2          # xmm2  = s2 * a(i6)
    mulpd     %xmm3 , %xmm4          # xmm4  = c3 * a(i3)
    mulpd     %xmm15, %xmm3          # xmm3  = s3 * a(i3)
    addpd     %xmm0 , %xmm9          # xmm9  = a1p5
    mulpd     %xmm11, %xmm1          # xmm1  = s1 * a(i1)
    movapd    %xmm13, %xmm0          # xmm0  = s2
    mulpd     %xmm6 , %xmm0          # xmm0  = s2 * a(i6)
    mulpd     %xmm12, %xmm6          # xmm6  = c2 * a(i6)
    mulpd     %xmm10, %xmm5          # xmm5  = c1 * a(i5)
    addpd     %xmm0 , %xmm8          # xmm8  = a2p6
    movapd    %xmm15, %xmm0          # xmm0  = s3
    mulpd     %xmm7 , %xmm0          # xmm0  = s3 * a(i7)
    subpd     %xmm2 , %xmm6          # xmm6  = a6m2
    movapd    (%rdi), %xmm2          # xmm2  = a(i0)
    mulpd     %xmm14, %xmm7          # xmm7  = c3 * a(i7)
    subpd     %xmm1 , %xmm5          # xmm5  = a5m1
    movapd    %xmm6 , %xmm1          # xmm1  = a6m2
    addpd     %xmm4 , %xmm0          # xmm0  = a3p7
    subpd     %xmm3 , %xmm7          # xmm7  = a7m3
    movapd    %xmm8 , %xmm3          # xmm3  = a2p6
    addpd     %xmm2 , %xmm3          # xmm3  = a(i0) + a2p6
    subpd     %xmm8 , %xmm2          # xmm2  = a(i0) - a2p6
    movapd    %xmm2 , %xmm8          # xmm8  = a2
    movapd    %xmm9 , %xmm2          # xmm2  = a1p5
    addpd     %xmm0 , %xmm2          # xmm2  = a1
    subpd     %xmm9 , %xmm0          # xmm0  = a3
    movapd    (%rdi,%rbp,4), %xmm9   # xmm9  = a(i4)
    addpd     %xmm9 , %xmm1          # xmm1  = b0
    subpd     %xmm6 , %xmm9          # xmm9  = b2
    movapd    %xmm9 , %xmm6          # xmm6  = b2
    movapd    %xmm5 , %xmm9          # xmm9  = a5m1
    addpd     %xmm7 , %xmm9          # xmm9  = b1
    subpd     %xmm7 , %xmm5          # xmm5  = b3
    movapd    %xmm3 , %xmm7          # xmm7  = a0
    addpd     %xmm2 , %xmm7          # xmm7  = a0+a1
    subpd     %xmm2 , %xmm3          # xmm3  = a0-a1
    movapd    %xmm7 , (%rsi)         # c(j0)
    movapd    %xmm1 , %xmm7          # xmm7  = b0
    addpd     %xmm9 , %xmm7          # xmm7  = b0+b1
    movapd    %xmm3 , (%r10)         # c(j2)
    subpd     %xmm1 , %xmm9          # xmm9  = b1-b0
    movapd    %xmm7 , (%rsi,%rbp)    # c(j4) = c(j0+la)
    movapd    %xmm9 , (%r10,%rbp)    # c(j6)
    movapd    %xmm8 , %xmm9          # xmm9  = a2
    subpd     %xmm5 , %xmm8          # xmm8  = a2-b3
    addpd     %xmm5 , %xmm9          # xmm9  = a2+b3
    movapd    %xmm9 , (%r9)          # c(j1)
    movapd    %xmm0 , %xmm9          # xmm9  = a3
    subpd     %xmm6 , %xmm0          # xmm0  = a3-b2
    addpd     %xmm6 , %xmm9          # xmm9  = a3+b2
    movapd    %xmm9 , (%r9 ,%rbp)    # c(j5)
    movapd    %xmm8 , (%r8)          # c(j3)
    movapd    %xmm0 , (%r8 ,%rbp)    # c(j7)
    addq      $16, %rdi
    addq      $16, %rsi
    addq      $16, %r8
    addq      $16, %r9
    addq      $16, %r10
    subq      $2 , %rcx
    jnz       LOOP_DFFT4M_I

    addq    %r12, %rdi               # rdi    += la*8*7
    addq    %rbp, %rsi               # &c(j0) += la
    addq    %rbp, %r9                # &c(j1) += la
    subq    %r13, %r10               # r10    -= la*8*3
    subq    %r13, %r8                # r8     -= la*8*3

    addq    %rbp, %r14
    cmpq    %r14, %r11

    jg      LOOP_DFFT4M_O

    movq   $0x3FE6A09E667F3BCD, %rax # rax   = sin45 = qrt(0.5)
    movd    %rax , %xmm3             # xmm3  = sin45
    movddup %xmm3, %xmm3
    movq   $0x8000000000000000, %rax # rax   = sign bit
    movd    %rax , %xmm5             # xmm5  = sign bit
    movddup %xmm5, %xmm5             # xmm5  = sign bits 
    movq    %rbx, %rcx               # rcx   = la

LOOP_DFFT4M_B:

    movapd   (%rdi,%r13), %xmm2      # xmm2  = a(i3)
    movapd   (%rdi,%rbp), %xmm0      # xmm0  = a(i1)
    movapd   %xmm0, %xmm1            # xmm1  = a(i1)
    subpd    %xmm2, %xmm0            # xmm0  = a(i1) - a(i3)
    addpd    %xmm2, %xmm1            # xmm1  = a(i1) + a(i3)
    movapd   (%rdi), %xmm2           # xmm2  = a(i0)
    mulpd    %xmm3, %xmm0            # xmm0  = a1m3
    mulpd    %xmm3, %xmm1            # xmm1  = a1p3
    movapd   %xmm0, %xmm4            # xmm4  = a1m3
    addpd    %xmm2, %xmm4            # xmm4  = a(i0) + a1m3
    subpd    %xmm0, %xmm2            # xmm2  = a(i0) - a1m3
    movapd   (%rdi,%rbp,2), %xmm0    # xmm0  = a(i2)
    movapd   %xmm4, (%rsi)           # c(j0)
    movapd   %xmm2, (%r9 )           # c(j1)
    movapd   %xmm0, %xmm2            # xmm2  = a(i2)
    subpd    %xmm1, %xmm0            # xmm0  = a(i2) - a1p3
    xorpd    %xmm5, %xmm2            # xmm2  =-a(i2)
    subpd    %xmm1, %xmm2            # xmm2  =-a(i2) - a1p3 
    movapd   %xmm2, (%rsi,%rbp)      # c(j4)
    movapd   %xmm0, (%r9 ,%rbp)      # c(j5)
    addq     $16, %rdi
    addq     $16, %rsi
    addq     $16, %r9 
    subq     $2, %rcx
    jne      LOOP_DFFT4M_B

    addq     %r13, %rdi              # rdi   = next line
    leaq     (%r11,%r11,2), %rax     # rax   = n*3
    leaq     (%rsi,%rax,2), %rsi     # rsi   = nextline

    dec      %r15
    jnz      LOOP_DFFT4M_LOT

    popq     %rbx
    popq     %rbp
    popq     %r12
    popq     %r13
    popq     %r14
    popq     %r15
    ret


# Fast Double Precision Fourier Transposition
# ==========================================
# E. Kirk - 31-Aug-2015
# ==========================================
# call fast_gtp(a  ,b  ,n  )
#              (rdi,rsi,rdx)
#
# rdi: a(*,*)    : complex fourier coefficients
# rsi: b(n,n)    : square matrix
# rdx: n         : dimension
# ------------------------------------------

   .globl _fast_gtp
_fast_gtp:
   .globl fast_gtp_
fast_gtp_:
   .globl _fast_gtp_
_fast_gtp_:

   push    %rbx
   push    %r12
   push    %r13
   push    %r14
   push    %r15

   movq    $0x8000000000000000, %rax # rax   = sign bit
   movd    %rax, %xmm7               # xmm7  = sign bit
   movapd  %xmm7, %xmm6              # xmm6  = sign bit
   shufpd  $1, %xmm7, %xmm7          # swap cat words
   movddup %xmm6, %xmm6              # xmm6  = both sign bits set

   movl    (%rdx), %eax              # rax  = n
   movq    %rax, %r11                # r11  = n
   leaq    (,%rax,8), %r8            # r8   = n*8
   movl    $3, %ebx                  # rbx  = 3
   xorq    %rdx, %rdx                # rdx  = 0
   divl    %ebx                      # rax  = n/3
   movq    %rax, %rbx                # rbx  = k
   inc     %rax                      # rax  = k + 1
   leaq    (%rax,%rax), %r10         # r10  = k*2
   leaq    (,%r10,8), %r9            # r9   = rowsize fc

   leaq    (%rbx,%rbx), %rax         # rax  = k*2
   mulq    %r9                       # rax  * rowsize
   leaq    (%rsi,%rax), %r15         # r15  = &fc(lastrow)

   movq    %rdi, %r12
   leaq    2(%rbx,%rbx), %rcx        # rcx  = loop count k*2+2

# transpose col 0 to row 0

LOOP_1_FAST_GTP:
   movsd   (%r12), %xmm1             # xmm1 = a(0,i)
   movsd   %xmm1, (%rsi)             # b(i,0) = xmm1
   addq    %r8, %r12                 # r12+= row
   addq    $8, %rsi                  # &b += 8 bytes
   dec     %rcx
   jnz     LOOP_1_FAST_GTP


# transpose row 0 to col 0

   addq    $16, %rdi                 # start with gp(3-4,1)
   movq    %rsi, %r13                # start with fc(0-1,1)
   movq    %r15,  %r14               # start with fc(0-1,2k-1)
   movq    %rbx, %rcx                # loop count k

LOOP_2_FAST_GTP:
   movapd  (%rdi), %xmm1             # xmm1 = gp(j,1)
   movapd  %xmm1, (%r13)             # fc(0,j) = xmm1
   xorpd   %xmm7, %xmm1              # conjugate
   movapd  %xmm1, (%r14)             # fc(0,j) = xmm1

   addq    $16, %rdi                 # next number
   addq    %r9, %r13                 # row++
   subq    %r9, %r14                 # row--

   dec     %rcx
   jnz     LOOP_2_FAST_GTP

# transpose remaining matrix
# --------------------------
# rdi,r12 : read  index 
# r13,r14 : write index
# r8      : rowsize gp
# r9      : rowsize fc
# rcx     : inner loop
# rdx     : outer loop
# rbx     : k-1
# r10     : k*2
# r11     : n
# r15     : start of fc(lastrow)
# rax     : (n - k*2) * 8 + rowsize

   movq    %r11, %rax                # n
   subq    %r10, %rax                # n - k*2
   leaq    16(%r8 ,%rax,8), %rax     # row skip value
   addq    %rax, %rdi                # gp(3-4,j)
   leaq    (%rdi,%r8 ), %r12         # gp(3-4,j+1)
   movq    %rbx, %rdx                # loop count k
    
LOOP_O_FAST_GTP:
   leaq    16(%rsi), %r13            # start with fc(j,2-3)
   leaq    16(%r15), %r14            # start with fc(j,2-3)
   movq    %rbx, %rcx                # loop count k

LOOP_I_FAST_GTP:
   movapd   (%rdi), %xmm1            # xmm1 = gp(j,i)
   movapd   (%r12), %xmm2            # xmm2 = gp(j,i+1)
   movapd   %xmm1 , %xmm3            # xmm3 = gp(j,i)
   shufpd   $1, %xmm2, %xmm2         # swap cat words
   addsubpd %xmm2, %xmm1             # xmm1 = +k
   xorpd    %xmm6, %xmm3             # xmm3 = -gp(j,i)
   addsubpd %xmm3, %xmm2             # xmm2 = -k
   movapd   %xmm1, (%r13)
   movapd   %xmm2, (%r14)

   addq    $16, %rdi                 # next number
   addq    $16, %r12                 # next number
   addq    %r9, %r13                 # row++
   subq    %r9, %r14                 # row--

   dec     %rcx
   jnz     LOOP_I_FAST_GTP

   addq    %rax, %rdi
   addq    %rax, %r12
   addq    $16 , %rsi                # col++
   addq    $16 , %r15                # col++

   dec     %rdx
   jnz     LOOP_O_FAST_GTP

   pop     %r15
   pop     %r14
   pop     %r13
   pop     %r12
   pop     %rbx

   ret


