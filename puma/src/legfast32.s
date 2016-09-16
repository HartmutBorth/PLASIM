# Fast Direct Legendre Transformation
# ===================================
# E. Kirk - 13-Dec-2011
# ===================================
# call fc2sp(fc,sp)  single precision
#           (di,si)     REAL (kind=4)
#
# rdi: complex fc(nlon/2,nhpp)
# rsi: complex sp(ncsp)
# rdx: real    qc(ncsp,nhpp)
# ecx: m loop index
# eax: 2 * (nlon-ntp1) * sizeof(REAL)
# ebx: n loop index
# ebp: ntp1
#
# The xmm register are used to store four
# single precision floating point numbers
# in the following order:
# ---------------------------------------
# | 0 |     symmetric mode real      part
# | 1 |     symmetric mode imaginary part
# | 2 | antisymmetric mode real      part
# | 3 | antisymmetric mode imaginary part
# ---------------------------------------

# Define two entry points to satisfy linker hich generate
# a leading and trailing underscore or a trailing one only.

.globl _fc2sp_
_fc2sp_:
.globl fc2sp_
fc2sp_:
        push      %rbx
        push      %rbp
        mov       LEG_QC(%rip),%rdx
        mov       LEG_NHPP(%rip),%r8d
        mov       LEG_NTP1(%rip),%ebp
        mov       LEG_NLON(%rip),%r9d
        mov       LEG_SPM1(%rip),%ecx
        xorps     %xmm1,%xmm1      # xmm1 = 0
        xorps     %xmm2,%xmm2      # xmm2 = 0
L_INIT_FC2SP:
        movlps    %xmm1,(%rsi,%rcx,8)
        dec       %rcx
        jns       L_INIT_FC2SP
        movl      %r9d,%eax       # eax  = nlon
        shlq      $2,%r9          # r9   = nlon * 4
        sub       %ebp,%eax       # eax  = nlon - ntp1
        shl       $3,%eax         # eax  = 8 * (nlon-ntp1)
LOOP_LAT_FC2SP:
        push     %rsi             # save for next lat
        movl     %ebp,%ecx        # reset m = ntp1
LOOP_M_FC2SP:
        movddup  (%rdi),%xmm3     # xmm3 = North / North
        movlps   (%rdi,%r9),%xmm1 # xmm1 = South /   0
        movlhps  %xmm1,%xmm2      # xmm2 =  0    / South
        addps    %xmm1,%xmm3      # nc + sc (lo)
        subps    %xmm2,%xmm3      # nc - sc (hi)
        add      $8,%rdi          # fc += 2
        movl     %ecx, %ebx       # index n = m
LOOP_N_FC2SP:
        dec      %ebx             # --n
        jz       L_SINGLE_FC2SP
        movups  (%rsi),%xmm4      # sp(w)
        movlps  (%rdx),%xmm5      # qc(w) & qc(w+1)
        unpcklps %xmm5,%xmm5      # duplicate & combine qc's
        mulps    %xmm3,%xmm5      # qc(w) * (fc(m) + fc(m+nlat))
        addps    %xmm5,%xmm4      # sp(w) + qc(w) * fc
        movups   %xmm4,(%rsi)     # two modes
        add      $16,%rsi         # sp += 4
        add      $8,%rdx          # qc++
        dec      %ebx             # --n
        jnz      LOOP_N_FC2SP
        loop     LOOP_M_FC2SP
        jmp      LOOP_M_FC2SP_EXIT
L_SINGLE_FC2SP:
        movss   (%rsi),%xmm4      # sp(w) real
        movss   (%rdx),%xmm5      # qc(w)
        mulss    %xmm3,%xmm5      # qc(w) * (fc(m) + fc(m+nlat))
        addss    %xmm5,%xmm4      # sp(w) + qc(w) * fc
        movss    %xmm4,(%rsi)     # single mode
        addq     $4,%rsi          # sp += 2
        movss   (%rsi),%xmm4      # sp(w) imag
        movss   (%rdx),%xmm5      # qc(w)
        mulss    %xmm3,%xmm5      # qc(w) * (fc(m) + fc(m+nlat))
        addss    %xmm5,%xmm4      # sp(w) + qc(w) * fc
        movss    %xmm4,(%rsi)     # single mode
        addq     $4,%rsi          # sp += 2
        add      $4,%rdx          # qc++
        loop     LOOP_M_FC2SP     # while (m != 0)
LOOP_M_FC2SP_EXIT:
        add      %rax,%rdi        # fc += 2 * (nlon-ntp1)
        pop      %rsi             # reset sp
        dec      %r8d             # --l
        jnz      LOOP_LAT_FC2SP   # while (l != 0)
        pop      %rbp
        pop      %rbx
        ret

# Fast Indirect Legendre Transformation
# =====================================
# E. Kirk - 17-Nov-2011
# =====================================
# call sp2fc(sp,fc)
#           (di,si)
#
# rsi: complex sp(ncsp)
# rdi: complex fc(nlon/2,nlpp)
# rdx: real    qi(ncsp,nlpp)
# ebx: n loop index
# ecx: m loop index
# ebp: ntp1
# eax: nlon / 2 - ntp1
#
#do l = 1 , nhpp
#  w = 1   
#  do m = 1 ,ntp1
#    e = w + ntp1 - m 
#    fs = dot_product(qi(w  :e:2,l),sp(w  :e:2))
#    fa = dot_product(qi(w+1:e:2,l),sp(w+1:e:2))
#    fc(m     ,l) = fs + fa
#    fc(m+nlat,l) = fs - fa
#    w = e + 1 
#  enddo ! m 
#enddo ! l 

.globl _sp2fc_
_sp2fc_:
.globl sp2fc_
sp2fc_:
        push     %rbx             # callee saved
        push     %rbp             # callee saved
        xchg     %rdi,%rsi        # si = sp, di = fc
        mov      LEG_QI(%rip),%rdx
        mov      LEG_NTP1(%rip),%ebp
        mov      LEG_NLON(%rip),%r9d
        mov      LEG_NHPP(%rip),%r8d
        movl     %r9d,%eax        # eax = nlon
        shl      $2,%r9           # nlon * 4
        shr      $1,%eax          # nlon / 2
        sub      %ebp,%eax        # eax = nlon / 2 - ntp1
        xorps    %xmm0,%xmm0      # remains zero
LOOP_LAT_SP2FC:
        push     %rsi             # save sp(:)
        movl     %ebp,%ecx        # reset m = ntp1
LOOP_M_SP2FC:
        xorps    %Xmm3,%Xmm3      # prepare sum
        movl     %ecx, %ebx       # index n = m
LOOP_N_SP2FC:
        dec      %ebx             # --n
        jz      L_SINGLE_SP2FC
        movups  (%rsi),%xmm2      # sp(w) & sp(w+1)
        movlps  (%rdx),%xmm1      # qi(w) & qi(w+1)
        unpcklps %xmm1,%xmm1      # duplicate & combine qi's
        mulps    %xmm1,%xmm2      # qi(w) * sp(w)
        addps    %xmm2,%Xmm3      # sum += qi(w) * sp(w)
        add      $16,%rsi         # sp += 4
        add      $8,%rdx          # qi++
        dec      %ebx             # --n
        jnz      LOOP_N_SP2FC
        jmp      LOOP_N_SP2FC_EXIT
L_SINGLE_SP2FC:                   # last single mode in row
        movss   (%rdx),%xmm1      # qi(w)
        unpcklps %xmm1,%xmm1      # duplicate qi(w)
        movlps  (%rsi),%xmm2      # sp(w)
        mulps    %xmm1,%xmm2      # qi(w) * sp(w)
        shufps  $228,%xmm0,%xmm2  # fa = 0 (upper quad)
        addps   %xmm2,%Xmm3       # sum += qi(w) * sp(w)
        add     $8,%rsi           # sp += 2
        add     $4,%rdx           # qc++
LOOP_N_SP2FC_EXIT:
        movhlps %Xmm3,%xmm1       # fa
        movq    %Xmm3,%xmm2       # fs
        addps   %xmm1,%xmm2       # fs + fa
        subps   %xmm1,%Xmm3       # fs - fa
        movlps  %xmm2,(%rdi)      # fc(m     ,l) = fs + fa
        movlps  %Xmm3,(%rdi,%r9)  # fc(m+nlon,l) = fs - fa
        add     $8,%rdi 
        loop    LOOP_M_SP2FC      # while (m != 0)
        mov     %eax,%ecx         # set remaining fc's to zero
LOOP_CLEAR_SP2FC:
        movlps  %xmm0,(%rdi)      # fc north = 0.0
        movlps  %xmm0,(%rdi,%r9)  # fc south = 0.0
        add     $8,%rdi           # sizeof(complex)
        loop    LOOP_CLEAR_SP2FC
        add     %r9,%rdi          # skip southern lat
        pop     %rsi              # reset sp(:)
        dec     %r8d              # --l
        jnz     LOOP_LAT_SP2FC    # while (l != 0)
        pop     %rbp
        pop     %rbx
        ret
# Fast Indirect Legendre Transformation
# =====================================
# E. Kirk - 17-Nov-2011
# =====================================
# call sp2fcdmu(sp,fc)
#              (di,si)
#
# rsi: complex sp(ncsp)
# rdi: complex fc(nlon/2,nlpp)
# rdx: real    qi(ncsp,nlpp)
# ebx: n loop index
# ecx: m loop index
# ebp: ntp1
# eax: nlon / 2 - ntp1
#
#do l = 1 , nhpp
#  w = 1   
#  do m = 1 ,ntp1
#    e = w + ntp1 - m 
#    fs = dot_product(qi(w  :e:2,l),sp(w  :e:2))
#    fa = dot_product(qi(w+1:e:2,l),sp(w+1:e:2))
#    fc(m     ,l) = fs + fa
#    fc(m+nlat,l) = fs - fa
#    w = e + 1 
#  enddo ! m 
#enddo ! l 

.globl _sp2fcdmu_
_sp2fcdmu_:
.globl sp2fcdmu_
sp2fcdmu_:
        push    %rbx              # callee saved
        push    %rbp              # callee saved
        xchg    %rdi,%rsi         # si = sp, di = fc
        mov     LEG_QJ(%rip),%rdx
        mov     LEG_NTP1(%rip),%ebp
        mov     LEG_NLON(%rip),%r9d
        mov     LEG_NHPP(%rip),%r8d
        mov     %r9d,%eax         # eax = nlon
        shl     $2,%r9            # nlon * 4
        shr     $1,%eax           # nlon / 2
        sub     %ebp,%eax         # eax = nlon / 2 - ntp1
        xorps   %xmm0,%xmm0       # remains zero
LOOP_LAT_DMU:
        push     %rsi             # save for next lat
        movl     %ebp,%ecx        # reset m = ntp1
LOOP_M_DMU:
        xorps    %Xmm3,%Xmm3      # prepare sum
        movl     %ecx, %ebx       # index n = m
LOOP_N_DMU:
        dec      %ebx             # --n
        jz      L_SINGLE_DMU
        movlps  (%rdx),%xmm1      # qj(w) & qj(w+1)
        unpcklps %xmm1,%xmm1      # duplicate & combine qj's
        movups  (%rsi),%xmm2      # sp(w) & sp(w+1)
        mulps    %xmm1,%xmm2      # qj(w) * sp(w)
        addps    %xmm2,%Xmm3      # sum += qj(w) * sp(w)
        add      $16,%rsi         # sp += 4
        add      $8,%rdx          # qj++
        dec      %ebx             # --n
        jnz      LOOP_N_DMU
        jmp      LOOP_N_DMU_EXIT
L_SINGLE_DMU:                     # last single mode in row
        movss   (%rdx),%xmm1      # qj(w)
        unpcklps %xmm1,%xmm1      # duplicate qj(w)
        movlps  (%rsi),%xmm2      # sp(w)
        mulps    %xmm1,%xmm2      # qj(w) * sp(w)
        shufps  $228,%xmm0,%xmm2  # fa = 0 (upper quad)
        addps   %xmm2,%Xmm3       # sum += qj(w) * sp(w)
        add     $8,%rsi           # sp += 2
        add     $4,%rdx           # qj++
LOOP_N_DMU_EXIT:
        movhlps %Xmm3,%xmm1       # fa
        movq    %Xmm3,%xmm2       # fs
        addps   %xmm1,%xmm2       # fs + fa
        subps   %xmm3,%Xmm1       # fa - fs
        movlps  %xmm2,(%rdi)      # fc(m     ,l) = fs + fa
        movlps  %Xmm1,(%rdi,%r9)  # fc(m+nlon,l) = fa - fs
        add     $8,%rdi 
        loop    LOOP_M_DMU        # while (m != 0)
        mov     %eax,%ecx         # set remaining fc's to zero
LOOP_CLEAR_DMU:
        movlps  %xmm0,(%rdi)      # fc north = 0.0
        movlps  %xmm0,(%rdi,%r9)  # fc south = 0.0
        add     $8,%rdi           # sizeof(complex)
        loop    LOOP_CLEAR_DMU
        add     %r9,%rdi          # skip southern lat
        pop     %rsi              # reset sp(:)
        dec     %r8d              # --l
        jnz     LOOP_LAT_DMU      # while (l != 0)
        pop     %rbp
        pop     %rbx
        ret

# Fast Legendre Transformation
# ============================
# E. Kirk - 17-Nov-2011
# ============================
# call dv2uv(pd,pz,pu,pv,db)
# register: (di,si,dx,cx,r8)
#
# rdi: complex pd(ncsp)
# rsi: complex pz(ncsp)
# r12: complex pu(nlon,nhpp)
# r13: complex pv(nlon,nhpp)
# r14: real    qu(ncsp,nhpp)
# r15: real    qv(ncsp,nhpp)
# ecx: m loop index
# ebp: ntp1
# eax: nlon / 2 - ntp1
#
# xmm0: always 0
# xmm1: save pz(3) = absolute vorticity

.globl _dv2uv_
_dv2uv_:
.globl dv2uv_
dv2uv_:
        push    %r15
        push    %r14
        push    %r13
        push    %r12
        push    %r11
        push    %r10
        push    %rbx              # callee saved
        push    %rbp              # callee saved
        mov     %rdx,%r12         # r12 = pu
        mov     %rcx,%r13         # r13 = pv

        mov     LEG_QU(%rip),%r14
        mov     LEG_QV(%rip),%r15
        mov     LEG_NTP1(%rip),%ebp
        mov     LEG_NHPP(%rip),%r8d
        mov     LEG_NLON(%rip),%r9d
        movss   LEG_PVOR(%rip),%xmm2

        movss   8(%rsi),%xmm1     # xmm1 = pz(3)
        movss   %xmm1,%xmm3       # xmm3 = pz(3)
        subss   %xmm2,%xmm3       # xmm3 = pz(3) - plavor
        movss   %xmm3,8(%rsi)     # pz(3)= pz(3) - plavor

        mov     %r9d,%eax         # eax = nlon
        shl     $2,%r9            # nlon * 4
        shr     $1,%eax           # nlon / 2
        sub     %ebp,%eax         # eax = nlon / 2 - ntp1
        xorps   %xmm0,%xmm0       # remains zero
LOOP_LAT_DV2UV:                   # do r8d = nhpp , 1 , -1
        xor      %r10,%r10        # index = 0
        mov      %ebp,%ecx
LOOP_M_DV2UV:                     # do ecx = ntp1 , 1 , -1
        xorps    %Xmm2,%Xmm2      # sum(ud) = 0
        xorps    %Xmm3,%Xmm3      # sum(vd) = 0
        xorps    %Xmm4,%Xmm4      # sum(uz) = 0
        xorps    %Xmm5,%Xmm5      # sum(vz) = 0
        mov      %ecx, %ebx       # index n = m

LOOP_N_DV2UV:                     # do ebx = ecx  , 1 , -1
        cmp      $1,%ebx          # last single mode ?
        jne      L_DUAL_DV2UV
        xorps    %xmm6,%xmm6      # xmm6 = 0
        xorps    %xmm7,%xmm7      # xmm7 = 0
        xorps    %xmm8,%xmm8      # xmm8 = 0
        xorps    %xmm9,%xmm9      # xmm9 = 0
        movss    (%r14),%xmm6     # xmm6 = qu (lo)
        movss    (%r15),%xmm7     # xmm7 = qv (lo)
        movlps   (%rdi,%r10,8),%xmm8 # d(:)
        movlps   (%rsi,%r10,8),%xmm9 # z(:)
        add      $4,%r14          # qu++
        add      $4,%r15          # qv++
        inc      %r10
        jmp      L_STEP_N_DV2UV
L_DUAL_DV2UV:
        movlps   (%r14),%xmm6     # xmm6 = qu
        movlps   (%r15),%xmm7     # xmm7 = qv
        movups   (%rdi,%r10,8),%xmm8 # d(:)
        movups   (%rsi,%r10,8),%xmm9 # z(:)
        add      $8,%r14          # qu += 2
        add      $8,%r15          # qv += 2
        add      $2,%r10
        dec      %ebx
L_STEP_N_DV2UV:
        unpcklps %xmm6,%xmm6      # duplicate & combine qu's
        unpcklps %xmm7,%xmm7      # duplicate & combine qv's
        movaps   %xmm6,%xmm10     # u
        mulps    %xmm8,%xmm10     # ud
        addsubps %xmm10,%xmm2     # sum(ud) +ai -ar +si -sr
        mulps    %xmm7,%xmm8      # vd
        subps    %xmm8,%xmm3      # sum(-vd)
        mulps    %xmm9,%xmm6      # uz
        addsubps %xmm6,%xmm4      # sum(uz) +ai -ar +si -sr
        mulps    %xmm9,%xmm7      # vz
        addps    %xmm7,%xmm5      # sum(vz)
#       ------------------------- #
        dec      %ebx             # --n
        jnz      LOOP_N_DV2UV
#       ------------------------- #
        shufps   $114,%xmm5,%xmm5 # vz: +si +ai +sr +ar 1302 = 72
        shufps   $114,%xmm3,%xmm3 # vd: -si -ai -sr -ar 1302 = 72
        movaps   %xmm5,%xmm6      # vz
        haddps   %xmm3,%xmm5      # no: v(-ud) & u(vz)
        hsubps   %xmm3,%xmm6      # so: v(-ud) & u(vz)
#       ------------------------- #
        shufps   $141,%xmm2,%xmm2 # ud: -ar -sr +ai +si 2031
        shufps   $141,%xmm4,%xmm4 # uz: -ar -sr +ai +si 2031
        movaps   %xmm2,%xmm3
        haddps   %xmm4,%xmm3      # north ud & uz
        hsubps   %xmm4,%xmm2      # south ud & uz
        addps    %xmm3,%xmm5      # north  u & v
        addps    %xmm2,%xmm6      # south  u & v
#       ------------------------- #
        movlps   %xmm5,(%r12)     # pu(m     ,l)
        movhps   %xmm5,(%r13)     # pv(m     ,l)
        movlps   %xmm6,(%r12,%r9) # pu(m+nlon,l)
        movhps   %xmm6,(%r13,%r9) # pv(m+nlon,l)
        add      $8,%r12          # pu += 2
        add      $8,%r13          # pv += 2
        dec      %ecx
        jnz      LOOP_M_DV2UV     # while (m != 0)
        mov      %eax,%ecx        # set remaining fc's to zero
LOOP_CLEAR_DV2UV:
        movlps   %xmm0,(%r12)     
        movlps   %xmm0,(%r13)     
        movlps   %xmm0,(%r12,%r9)
        movlps   %xmm0,(%r13,%r9)
        add      $8,%r12          # pu += 2
        add      $8,%r13          # pv += 2
        loop     LOOP_CLEAR_DV2UV
        add      %r9,%r12         # skip southern lat
        add      %r9,%r13         # skip southern lat
        dec      %r8d             # --l
        jnz      LOOP_LAT_DV2UV         # while (l != 0)
        movss    %xmm1,8(%rsi)    # restore pz(3)
        pop      %rbp
        pop      %rbx
        pop      %r10
        pop      %r11
        pop      %r12
        pop      %r13
        pop      %r14
        pop      %r15
        ret

# Fast Direct Legendre Transformation
# =======================================
# E. Kirk - 17-Nov_2011
# =======================================
# call mktend(d ,t ,z ,tn,fu,fv,ke,ut,vt)
#            (di,si,dx,cx,r8,r9,16,24,32)
#
# arg name  par (bp)  xreg
# -------------------------
#   1   d   rdi  -48  xmm7
#   2   t   rsi  -40  xmm8
#   3   z   rdx  -32  xmm9
#   4   tn  rcx  -24  xmm10
#   5   fu  r8   -16  xmm11
#   6   fv  r9    -8  xmm12
#   7   ke        16  xmm13
#   8   ut        24  xmm14
#   9   vt        32  xmm15
#
# ecx: m loop index
# eax: 2 * (nlon-ntp1) * sizeof(float)
# ebx: n loop index
# r15: ntp1
#
.globl _mktend_
_mktend_:
.globl mktend_
mktend_:
        push   %rbp               # frame pointer
        mov    %rsp,%rbp
        sub    $48,%rsp           # 6 local pointer
        push   %r15
        push   %r14
        push   %r13
        push   %r12
        push   %r11
        push   %r10
        push   %rbx

        mov    %rdi,-48(%rbp)     #  d(:)
        mov    %rsi,-40(%rbp)     #  t(:)
        mov    %rdx,-32(%rbp)     #  z(:)
        mov    %rcx,-24(%rbp)     # tn(:)
        mov    %r8 ,-16(%rbp)     # fu(:)
        mov    %r9 , -8(%rbp)     # fv(:)

        mov     LEG_QQ(%rip),%r14
        mov     LEG_QM(%rip),%r11
        mov     LEG_QE(%rip),%r12
        mov     LEG_QC(%rip),%r13
        mov     LEG_NHPP(%rip),%r8d
        mov     LEG_NLON(%rip),%r9d
        mov     LEG_NTP1(%rip),%r15d
        mov     LEG_SPM1(%rip),%ecx
#       ------------------------- #
        xorps     %xmm0,%xmm0
L_INIT_MKTEND:
        movlps    %xmm0,(%rdi,%rcx,8) # d(:) = 0.0
        movlps    %xmm0,(%rsi,%rcx,8) # t(:) = 0.0
        movlps    %xmm0,(%rdx,%rcx,8) # z(:) = 0.0
        dec       %rcx
        jns       L_INIT_MKTEND
#       ------------------------- #
        movl      %r9d,%eax       # eax = nlon
        shlq      $2,%r9          # r9  =nlon * 4
        sub       %r15d,%eax      # eax = nlon - ntp1
        shl       $3,%eax         # eax = 8 * (nlon-ntp1)
        xorps     %xmm1,%xmm1     # xmm1 = 0
        xorps     %xmm2,%xmm2     # xmm2 = 0
LOOP_LAT_MKTEND:
        mov     -48(%rbp),%rdi    # reset d(:)
        mov     -40(%rbp),%rsi    # reset t(:)
        mov     -32(%rbp),%rdx    # reset z(:)
        movl     %r15d,%ecx       # reset m = ntp1
LOOP_M_MKTEND:
#       ------- tn -------------- #
        mov     -24(%rbp),%r10    # tn(:)
        movddup (%r10),%xmm10     # xmm10 = tn: North / North
        movlps  (%r10,%r9),%xmm1  # xmm1  = tn: South /   0
        movlhps %xmm1,%xmm2       # xmm2  = tn:  0    / South
        addps   %xmm1,%xmm10      # nc + sc (lo)
        subps   %xmm2,%xmm10      # nc - sc (hi)
        addq    $8,-24(%rbp)      # tn += 2
#       ------- fu -------------- #
        mov     -16(%rbp),%r10    # fu(:)
        movddup (%r10),%xmm11     # xmm11 = fu: North / North
        movlps  (%r10,%r9),%xmm1  # xmm1 = fu: South /   0
        movlhps %xmm1,%xmm2       # xmm2 = fu:  0    / South
        addps   %xmm1,%xmm11      # nc + sc (lo)
        subps   %xmm2,%xmm11      # nc - sc (hi)
        addq    $8,-16(%rbp)      # fu += 2
#       ------- fv -------------- #
        mov     -8(%rbp),%r10     # fv(:)
        movddup (%r10),%xmm12     # xmm12 = fv: North / North
        movlps  (%r10,%r9),%xmm1  # xmm1 = fv: South /   0
        movlhps %xmm1,%xmm2       # xmm2 = fv:  0    / South
        subps   %xmm1,%xmm12      # nc - sc (lo)
        addps   %xmm2,%xmm12      # nc + sc (hi)
        addq    $8,-8(%rbp)       # fv += 2
#       ------- ke -------------- #
        mov     16(%rbp),%r10     # ke(:)
        movddup (%r10),%xmm13     # xmm13 = ke: North / North
        movlps  (%r10,%r9),%xmm1  # xmm1 = ke: South /   0
        movlhps %xmm1,%xmm2       # xmm2 = ke:  0    / South
        addps   %xmm1,%xmm13      # nc + sc (lo)
        subps   %xmm2,%xmm13      # nc - sc (hi)
        addq    $8,16(%rbp)       # ke += 2
#       ------- ut -------------- #
        mov     24(%rbp),%r10     # ut(:)
        movddup (%r10),%xmm5      # xmm5 = ut: North / North
        movlps  (%r10,%r9),%xmm1  # xmm1 = ut: South /   0
        movlhps %xmm1,%xmm2       # xmm2 = ut:  0    / South
        addps   %xmm1,%xmm5       # nc + sc (lo)
        subps   %xmm2,%xmm5       # nc - sc (hi)
        shufps  $177,%xmm5,%xmm5  # real <-> imag  2:3:0:1 = b1
        xorps   %xmm14,%xmm14     # xmm14 = ut
        subps   %xmm5,%xmm14
        addq    $8,24(%rbp)       # ut += 2
#       ------- vt -------------- #
        mov     32(%rbp),%r10     # vt(:)
        movddup (%r10),%xmm15     # xmm15 = vt: North / North
        movlps  (%r10,%r9),%xmm1  # xmm1  = vt: South /   0
        movlhps %xmm1,%xmm2       # xmm2  = vt:  0    / South
        subps   %xmm1,%xmm15      # nc - sc (lo)
        addps   %xmm2,%xmm15      # nc + sc (hi)
        addq    $8,32(%rbp)       # vt += 2
#       ------------------------- #
        movl    %ecx, %ebx        # index n = m
        cmp     $1, %ecx          # last mode ?
        jne     LOOP_N_MKTEND     # proceed
#       ------------------------- #
        movsd   (%rdi),%xmm7      # d(w)
        movsd   (%rsi),%xmm8      # t(w)
        movsd   (%rdx),%xmm9      # z(w)
        jmp     LOAD_QE
        
LOOP_N_MKTEND:
        movups  (%rdi),%xmm7      # d(w) / d(w+1)
        movups  (%rsi),%xmm8      # t(w) / t(w+1)
        movups  (%rdx),%xmm9      # z(w) / z(w+1)
LOAD_QE:
#       ------- load qe---------- #
        movlps  (%r12),%xmm3      # qe(w) & qe(w+1)
        unpcklps %xmm3,%xmm3      # duplicate & combine qe's
#       ------- load qm --------- #
        movlps  (%r11),%xmm4      # qm(w) & qm(w+1)
        unpcklps %xmm4,%xmm4      # duplicate & combine qm's
#       ------- qq * ke --------- #
        movlps  (%r14),%xmm5      # qq(w) & qq(w+1)
        unpcklps %xmm5,%xmm5      # duplicate & combine qq's
        mulps    %xmm13,%xmm5     # qq(w) * (ke(m) + ke(m+nlat))
        addps    %xmm5,%xmm7      # d(w) + qq(w) * ke
#       ------- qe * fv---------- #
        movaps   %xmm3,%xmm5      # qe(w) & qe(w+1)
        mulps    %xmm12,%xmm5     # qe(w) * (fv(m) + fv(m+nlat))
        subps    %xmm5,%xmm7      # d(w) - qe(w) * fv(a)
#       ------- qm * fu --------- #
        movaps   %xmm4,%xmm5      # qm(w) & qm(w+1)
        mulps    %xmm11,%xmm5     # qm(w) * (ke(m) + ke(m+nlat))
        shufps  $177,%xmm5,%xmm5  # real <-> imag  2:3:0:1 = b1
        addsubps %xmm5,%xmm7      # d(w) -+-+ qm(w) * fu
#       ------- qe * vt---------- #
        movaps   %xmm3,%xmm5      # qe(w) & qe(w+1)
        mulps    %xmm15,%xmm5     # qe(w) * (vt(m) + vt(m+nlat))
        addps    %xmm5,%xmm8      # t(w) + qe(w) * vt
#       ------- qc * tn --------- #
        movlps  (%r13),%xmm5      # qc(w) & qc(w+1)
        unpcklps %xmm5,%xmm5      # duplicate & combine qc's
        mulps    %xmm10,%xmm5     # qc(w) * (tn(m) + tn(m+nlat))
        addps    %xmm5,%xmm8      # t(w) + qc(w) * tn
#       ------- qm * ut --------- #
        movaps   %xmm4,%xmm5      # qm(w) & qm(w+1)
        mulps    %xmm14,%xmm5     # qm(w) * (ut(m) + ut(m+nlat))
        addsubps %xmm5,%xmm8      # t(w) -+-+ qm(w) * ut
#       ------- qe * fu---------- #
        movaps   %xmm11,%xmm5     # fu
        shufps  $78,%xmm5,%xmm5   # symm <-> asym 1:0:3:2 = 4e 
        mulps    %xmm3,%xmm5      # qe(w) * (fu(m) + fu(m+nlat))
        addps    %xmm5,%xmm9      # z(w) + qe(w) * fu(a)
#       ------- qm * fv --------- #
        movaps   %xmm12,%xmm5     # fv
        shufps  $27,%xmm5,%xmm5   # 0:1:2:3 = 1b
        mulps    %xmm4,%xmm5      # qm(w) * (fv(m) + fv(m+nlat))
        addsubps %xmm5,%xmm9      # z(w) -+-+ zm(w) * fv
#       ------------------------- #
        dec      %ebx             # --n
        jz      L_SINGLE_MKTEND
        movups  %xmm7,(%rdi)      # two modes of d(:)
        movups  %xmm8,(%rsi)      # two modes of t(:)
        movups  %xmm9,(%rdx)      # two modes of z(:)
        add     $16,%rdi          # d += 4
        add     $16,%rsi          # t += 4
        add     $16,%rdx          # z += 4
        add     $8,%r14           # qq++
        add     $8,%r13           # qc++
        add     $8,%r12           # qe++
        add     $8,%r11           # qm++
        dec     %ebx              # --n
        jnz     LOOP_N_MKTEND
        dec     %ecx
        jnz     LOOP_M_MKTEND
        jmp     LOOP_M_MKTEND_EXIT
L_SINGLE_MKTEND:
        movlps  %xmm7,(%rdi)      # single mode of d(:)
        movlps  %xmm8,(%rsi)      # single mode of t(:)
        movlps  %xmm9,(%rdx)      # single mode of z(:)
        add     $8,%rdi           # d += 2
        add     $8,%rsi           # t += 2
        add     $8,%rdx           # z += 2
        add     $4,%r14           # qq++
        add     $4,%r13           # qc++
        add     $4,%r12           # qe++
        add     $4,%r11           # qm++
        dec     %ecx
        jnz     LOOP_M_MKTEND            # while (m != 0)
LOOP_M_MKTEND_EXIT:
        add     %rax,-24(%rbp)    # tn += 2 * (nlon-ntp1)
        add     %rax,-16(%rbp)    # fu += 2 * (nlon-ntp1)
        add     %rax, -8(%rbp)    # fv += 2 * (nlon-ntp1)
        add     %rax, 16(%rbp)    # ke += 2 * (nlon-ntp1)
        add     %rax, 24(%rbp)    # ut += 2 * (nlon-ntp1)
        add     %rax, 32(%rbp)    # vt += 2 * (nlon-ntp1)
        dec     %r8d              # --l
        jnz     LOOP_LAT_MKTEND   # while (l != 0)
        pop     %rbx              # restore register
        pop     %r10
        pop     %r11
        pop     %r12
        pop     %r13
        pop     %r14
        pop     %r15
        leave
        ret

# Fast Direct Legendre Transformation
# =======================================
# E. Kirk - 17-Nov_2011
# =======================================
# call qtend(q ,qn,uq,vq)
#           (di,si,dx,cx)
#
# arg name  par (bp)  xreg
# -------------------------
#   1   q   rdi
#   2   qn  rcx
#   3   uq  rdx
#   4   vq  rcx
#
# ecx: m loop index
# eax: 2 * (nlon-ntp1) * sizeof(float)
# ebx: n loop index
# r15: ntp1
#
.globl _qtend_
_qtend_:
.globl qtend_
qtend_:
        push   %r15
        push   %r14
        push   %r13
        push   %r12
        push   %r11
        push   %r10
        push   %rbx

        mov    %rcx,%r14          # vq(:)

        mov     LEG_QE(%rip),%r12
        mov     LEG_QM(%rip),%r11
        mov     LEG_QC(%rip),%r13
        mov     LEG_NHPP(%rip),%r8d
        mov     LEG_NLON(%rip),%r9d
        mov     LEG_SPM1(%rip),%ecx
        mov     LEG_NTP1(%rip),%r15d

        xorps     %xmm0,%xmm0
L_INIT_QTEND:
        movlps    %xmm0,(%rdi,%rcx,8) # q(:) = 0.0
        dec       %rcx
        jns       L_INIT_QTEND
        movl    %r9d,%eax         # eax = nlon
        shlq    $2,%r9            # r9  =nlon * 4
        sub     %r15d,%eax        # eax = nlon - ntp1
        shl     $3,%eax           # eax = 8 * (nlon-ntp1)
        xorps   %xmm1,%xmm1       # xmm1 = 0
        xorps   %xmm2,%xmm2       # xmm2 = 0
        xor     %r10,%r10         # r10  = 0
LOOP_LAT_QTEND:
        push    %rdi              # save q(:)
        movl     %r15d,%ecx       # reset m = ntp1
LOOP_M_QTEND:
#       ------- qn -------------- #
        movddup (%rsi),%xmm10     # xmm10 = qn: North / North
        movlps  (%rsi,%r9),%xmm1  # xmm1  = qn: South /   0
        movlhps %xmm1,%xmm2       # xmm2  = qn:  0    / South
        addps   %xmm1,%xmm10      # nc + sc (lo)
        subps   %xmm2,%xmm10      # nc - sc (hi)
        add     $8,%rsi           # qn += 2
#       ------- uq -------------- #
        movddup (%rdx),%xmm5      # xmm5 = uq: North / North
        movlps  (%rdx,%r9),%xmm1  # xmm1 = uq: Souqh /   0
        movlhps %xmm1,%xmm2       # xmm2 = uq:  0    / Souqh
        addps   %xmm1,%xmm5       # nc + sc (lo)
        subps   %xmm2,%xmm5       # nc - sc (hi)
        shufps  $177,%xmm5,%xmm5  # real <-> imag  2:3:0:1 = b1
        xorps   %xmm14,%xmm14     # xmm14 = uq
        subps   %xmm5,%xmm14
        add     $8,%rdx            # uq += 2
#       ------- vq -------------- #
        movddup (%r14),%xmm15     # xmm15 = vq: North / North
        movlps  (%r14,%r9),%xmm1  # xmm1  = vq: South /   0
        movlhps %xmm1,%xmm2       # xmm2  = vq:  0    / South
        subps   %xmm1,%xmm15      # nc - sc (lo)
        addps   %xmm2,%xmm15      # nc + sc (hi)
        add     $8,%r14           # vq += 2
#       ------------------------- #
        movl    %ecx, %ebx        # index n = m
LOOP_N_QTEND:
        movups  (%rdi),%xmm8      # q(w)
#       ------- qe * vq---------- #
        movlps  (%r12,%r10,4),%xmm3 # qe(w) & qe(w+1)
        unpcklps %xmm3,%xmm3      # duplicate & combine qe's
        mulps    %xmm15,%xmm3     # qe(w) * (vt(m) + vt(m+nlat))
        addps    %xmm3,%xmm8      # q(w) + qe(w) * vt
#       ------- qc * qn --------- #
        movlps  (%r13,%r10,4),%xmm4 # qc(w) & qc(w+1)
        unpcklps %xmm4,%xmm4      # duplicate & combine qc's
        mulps    %xmm10,%xmm4     # qc(w) * (tn(m) + tn(m+nlat))
        addps    %xmm4,%xmm8      # q(w) + qc(w) * tn
#       ------- qm * uq --------- #
        movlps  (%r11,%r10,4),%xmm5 # qm(w) & qm(w+1)
        unpcklps %xmm5,%xmm5      # duplicate & combine qm's
        mulps    %xmm14,%xmm5     # qm(w) * (ut(m) + ut(m+nlat))
        addsubps %xmm5,%xmm8      # q(w) -+-+ qm(w) * ut
#       ------------------------- #
        dec      %ebx             # --n
        jz      L_SINGLE_QTEND
        movups  %xmm8,(%rdi)      # two modes of q(:)
        add     $16,%rdi          # q += 4
        add     $2,%r10
        dec     %ebx              # --n
        jnz     LOOP_N_QTEND
        dec     %ecx
        jnz     LOOP_M_QTEND
        jmp     LOOP_M_QTEND_EXIT
L_SINGLE_QTEND:
        movlps  %xmm8,(%rdi)      # single mode of q(:)
        add     $8,%rdi           # q += 2
        inc     %r10
        dec     %ecx
        jnz     LOOP_M_QTEND            # while (m != 0)
LOOP_M_QTEND_EXIT:
        add     %rax,%rsi         # qn += 2 * (nlon-ntp1)
        add     %rax,%rdx         # uq += 2 * (nlon-ntp1)
        add     %rax,%r14         # vq += 2 * (nlon-ntp1)
        pop     %rdi              # reset q(:)
        dec     %r8d              # --l
        jnz     LOOP_LAT_QTEND   # while (l != 0)
        pop     %rbx              # restore register
        pop     %r10
        pop     %r11
        pop     %r12
        pop     %r13
        pop     %r14
        pop     %r15
        ret

# Fast Direct Legendre Transformation
# =======================================
# E. Kirk - 12-Dec_2011
# =======================================
# call uv2dv(u ,v ,d , z)
#           (di,si,dx,cx)
#
# arg name  par (bp)  xreg
# -------------------------
#   1   u   rdi
#   2   v   rsi
#   3   d   rdx
#   4   z   rcx
#
# eax: 2 * (nlon-ntp1) * sizeof(float)
# ebx: n loop index
# r14: m loop index
# r15: ntp1
#
.globl _uv2dv_
_uv2dv_:
.globl uv2dv_
uv2dv_:
        push   %r15
        push   %r14
        push   %r13
        push   %r12
        push   %r11
        push   %r10
        push   %rbx

        mov     LEG_QE(%rip),%r12
        mov     LEG_QM(%rip),%r11
        mov     LEG_NHPP(%rip),%r8d
        mov     LEG_NLON(%rip),%r9d
        mov     LEG_NTP1(%rip),%r15d
        mov     LEG_SPM1(%rip),%ebx
#       ------------------------- #
        xorps     %xmm0,%xmm0
L_INIT_UV2DV:
        movlps    %xmm0,(%rdx,%rbx,8) # d(:) = 0.0
        movlps    %xmm0,(%rcx,%rbx,8) # z(:) = 0.0
        dec       %rbx
        jns       L_INIT_UV2DV
#       ------------------------- #
        mov       %r9d,%eax       # eax = nlon
        shlq      $2,%r9          # r9  =nlon * 4
        sub       %r15d,%eax      # eax = nlon - ntp1
        shl       $3,%eax         # eax = 8 * (nlon-ntp1)
        xorps     %xmm1,%xmm1     # xmm1 = 0
        xorps     %xmm2,%xmm2     # xmm2 = 0
LOOP_LAT_UV2DV:
        push     %rdx             # save d(:)
        push     %rcx             # save z(:)
        mov      %r15d,%r14d      # reset m = ntp1
LOOP_M_UV2DV:
#       ------- u --------------- #
        movddup (%rdi),%xmm11     # xmm11 = fu: North / North
        movlps  (%rdi,%r9),%xmm1  # xmm1 = fu: South /   0
        movlhps %xmm1,%xmm2       # xmm2 = fu:  0    / South
        addps   %xmm1,%xmm11      # nc + sc (lo)
        subps   %xmm2,%xmm11      # nc - sc (hi)
        add     $8,%rdi           # fu += 2
#       ------- v --------------- #
        movddup (%rsi),%xmm12     # xmm12 = fv: North / North
        movlps  (%rsi,%r9),%xmm1  # xmm1 = fv: South /   0
        movlhps %xmm1,%xmm2       # xmm2 = fv:  0    / South
        addps   %xmm1,%xmm12      # nc + sc (lo)
        subps   %xmm2,%xmm12      # nc - sc (hi)
        add     $8,%rsi           # fv += 2
#       ------------------------- #
        movl    %r14d, %ebx       # index n = m
LOOP_N_UV2DV:
        movups  (%rdx),%xmm7      # d(w)
        movups  (%rcx),%xmm9      # z(w)
#       ------- load qe---------- #
        movlps  (%r12),%xmm3      # qe(w) & qe(w+1)
        unpcklps %xmm3,%xmm3      # duplicate & combine qe's
#       ------- load qm --------- #
        movlps  (%r11),%xmm4      # qm(w) & qm(w+1)
        unpcklps %xmm4,%xmm4      # duplicate & combine qm's
#       ------- qm * vs --------- #
        movaps   %xmm4,%xmm5      # qm(w) & qm(w+1)
        mulps    %xmm12,%xmm5     # qm(w) * v
        shufps  $177,%xmm5,%xmm5  # real <-> imag  2:3:0:1 = b1
        addsubps %xmm5,%xmm9      # z(w) -+-+ qm(w) * v(s)
#       ------- qm * us --------- #
        mulps    %xmm11,%xmm4     # qm(w) * u
        shufps  $177,%xmm4,%xmm4  # real <-> imag  2:3:0:1 = b1
        addsubps %xmm4,%xmm7      # d(w) -+-+ qm(w) * u(s)
#       ------- qe * ua---------- #
        movaps   %xmm3,%xmm5      # qe(w) & qe(w+1)
        movaps   %xmm11,%xmm6     #
        shufps  $78,%xmm6,%xmm6   # symm <-> asym 1:0:3:2 = 4e 
        mulps    %xmm6,%xmm5      # qe(w) * u(a)
        addps    %xmm5,%xmm9      # z(w) + qe(w) * u(a)
#       ------- qe * va---------- #
        movaps   %xmm12,%xmm6     #
        shufps  $78,%xmm6,%xmm6   # symm <-> asym 1:0:3:2 = 4e 
        mulps    %xmm6,%xmm3      # qe(w) * v
        subps    %xmm3,%xmm7      # d(w) - qe(w) * v(a)
#       ------------------------- #
        dec      %ebx             # --n
        jz      L_SINGLE_UV2DV
        movups  %xmm7,(%rdx)      # two modes of d(:)
        movups  %xmm9,(%rcx)      # two modes of z(:)
        add     $16,%rdx          # d += 4
        add     $16,%rcx          # z += 4
        add     $8,%r12           # qe++
        add     $8,%r11           # qm++
        dec     %ebx              # --n
        jnz     LOOP_N_UV2DV
        dec     %r14d
        jnz     LOOP_M_UV2DV
        jmp     LOOP_M_UV2DV_EXIT
L_SINGLE_UV2DV:
        movlps  %xmm7,(%rdx)      # single mode of d(:)
        movlps  %xmm9,(%rcx)      # single mode of z(:)
        add     $8,%rdx           # d += 2
        add     $8,%rcx           # z += 2
        add     $4,%r12           # qe++
        add     $4,%r11           # qm++
        dec     %r14d
        jnz     LOOP_M_UV2DV      # while (m != 0)
LOOP_M_UV2DV_EXIT:
        add     %rax,%rdi         # u += 2 * (nlon-ntp1)
        add     %rax,%rsi         # v += 2 * (nlon-ntp1)
        pop     %rcx              # restore z(:)
        pop     %rdx              # restore d(:)
        dec     %r8d              # --l
        jnz     LOOP_LAT_UV2DV    # while (l != 0)
        pop     %rbx              # restore register
        pop     %r10
        pop     %r11
        pop     %r12
        pop     %r13
        pop     %r14
        pop     %r15
        ret

# FORTRAN - Assember data interface
# ===================================
# E. Kirk - 13-Dec-2011
# ===================================
# subroutine legpar stores all needed
# constants.
#
# call legpar(ntp1,nlon,nhpp,spm1,pvor)
#            (rdi ,rsi ,rdx ,rcx ,r8  )

.globl _legpar_
_legpar_:
.globl legpar_
legpar_:

        mov     (%rdi),%eax
        mov     %eax,LEG_NTP1(%rip)  # truncation + 1
        mov     (%rsi),%eax
        mov     %eax,LEG_NLON(%rip)  # longitudes
        mov     (%rdx),%eax
        mov     %eax,LEG_NHPP(%rip)  # longitudes
        mov     (%rcx),%eax
        mov     %eax,LEG_SPM1(%rip)  # spectral modes - 1
        mov     (%r8),%eax
        mov     %eax,LEG_PVOR(%rip)  # planetary vorticity
        ret

# FORTRAN - Assember data interface
# ===================================
# E. Kirk - 13-Dec-2011
# ===================================
# subroutine legpolx stores addresses
# of polynomial arrays
#
# call legpola(qi ,qj ,qc ,spm1)
#             (rdi,rsi,rdx,rcx )

.globl _legpola_
_legpola_:
.globl legpola_
legpola_:

        mov     %rdi,LEG_QI(%rip)
        mov     %rsi,LEG_QJ(%rip)
        mov     %rdx,LEG_QC(%rip)
        mov     %rcx,LEG_QE(%rip)
        ret

.globl _legpolb_
_legpolb_:
.globl legpolb_
legpolb_:

        mov     %rdi,LEG_QM(%rip)
        mov     %rsi,LEG_QQ(%rip)
        mov     %rdx,LEG_QU(%rip)
        mov     %rcx,LEG_QV(%rip)
        ret

        .data
        .balign 8
LEG_QI: .quad 0
LEG_QJ: .quad 0
LEG_QC: .quad 0
LEG_QE: .quad 0
LEG_QM: .quad 0
LEG_QQ: .quad 0
LEG_QU: .quad 0
LEG_QV: .quad 0

LEG_PVOR:
        .long 0
LEG_NTP1:
        .long 0
LEG_NLON:
        .long 0
LEG_NHPP:
        .long 0
LEG_SPM1:
        .long 0
