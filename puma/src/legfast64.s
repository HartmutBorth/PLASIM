# Fast Double Precision Legendre Transformation
# =============================================
# E. Kirk - 14-Dec-2011
# =============================================
# call fc2sp(fc,sp,)
#            (di,si)
#
# rdi: fc(nlon,nlpp)
# rsi: sp(nesp)
# rdx: qc(ncsp,nlpp)
# ecx: m loop index
# eax: 2 * (nlon-ntp1) * sizeof(float)
# ebx: n loop index
# ebp: ntp1
#
# usage of the MME register xmm0 - xmm15:
# ---------------------------------------
# | 0 | spectral mode real      part
# | 1 | spectral mode imaginary part
# ---------------------------------------

.globl _fc2sp_
_fc2sp_:
.globl fc2sp_
fc2sp_:
        push     %rbx                
        push     %rbp                
        push     %r10
        mov      LEG_QC(%rip),%rdx   
        mov      LEG_NHPP(%rip),%r8d 
        mov      LEG_NLON(%rip),%r9d 
        mov      LEG_NTP1(%rip),%ebp 
        mov      LEG_SPM1(%rip),%ecx 
        shl      $1,%ecx             
        xorps    %xmm1,%xmm1         # xmm1 = 0
L_INIT_FC2SP:
        movups   %xmm1,(%rsi,%rcx,8) 
        sub      $2,%rcx             
        jns      L_INIT_FC2SP        
        movl     %r9d,%eax           # eax = nlon
        shlq     $3,%r9              # r9  =nlon * 8
        sub      %ebp,%eax           # eax = nlon - ntp1
        shl      $4,%eax             # eax = 16 * (nlon-ntp1)
LOOP_LAT_FC2SP:
        xor      %r10,%r10           # mode index
        mov      %ebp,%ecx           # reset m = ntp1
LOOP_M_FC2SP:
        movupd   (%rdi),%xmm1        # xmm1 = North fc
        movupd   (%rdi,%r9),%xmm3    # xmm3 = South fc
        movapd   %xmm1,%xmm2         # xmm2 = North fc
        addpd    %xmm3,%xmm1         # xmm1 = North + South
        subpd    %xmm3,%xmm2         # xmm2 = North - South
        add      $16,%rdi            # fc += 2
        mov      %ecx,%ebx           # index n = m
LOOP_N_FC2SP:
#       ---------------------------- # symmetric mode
        movddup  (%rdx,%r10),%xmm3   # xmm3 = qc(w)
        movupd   (%rsi,%r10,2),%xmm4 # xmm4 = sp(w)
        mulpd    %xmm1,%xmm3         # xmm3 = qc(w) * (fc(m) + fc(m+nlat))
        addpd    %xmm3,%xmm4         # xmm4 = sp(w) + qc(w) * fc(symm)
        movupd   %xmm4,(%rsi,%r10,2) # store sp(w)
        add      $8,%r10             # qc++
        dec      %ebx                # --n
        jz       L_SINGLE_FC2SP      
#       ---------------------------- # antisymmetric mode
        movddup  (%rdx,%r10),%xmm3   # xmm3 = qc(w)
        movupd   (%rsi,%r10,2),%xmm4 # xmm4 = sp(w)
        mulpd    %xmm2,%xmm3         # xmm3 = qc(w) * (fc(m) - fc(m+nlat))
        addpd    %xmm3,%xmm4         # xmm4 = sp(w) + qc(w) * fc(anti)
        movupd   %xmm4,(%rsi,%r10,2) # store sp(w)
        add      $8,%r10             # qc++
        dec      %ebx                # --n
        jnz      LOOP_N_FC2SP        
#       ---------------------------- # 
L_SINGLE_FC2SP:
        loop     LOOP_M_FC2SP        # while (m != 0)
LOOP_M_FC2SP_EXIT:
        add      %rax,%rdi           # fc += 2 * (nlon-ntp1)
        add      %r10,%rdx
        dec      %r8d                # --l
        jnz      LOOP_LAT_FC2SP      # while (l != 0)
        pop      %r10
        pop      %rbp                
        pop      %rbx                
        ret                          

# Fast Indirect Legendre Transformation
# =====================================
# E. Kirk - 25-Nov-2011
# =====================================
# call sp2fc(sp,fc) ,qi,ntp1,nhpp,nlon)
#           (di,si) ,dx,  cx,  r8,  r9)
#
# rsi: sp(nesp)
# rdi: fc(nlon,nlpp)
# rdx: qi(ncsp,nlpp)
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
        push     %rbx                
        push     %rbp                
        xchg     %rdi,%rsi           # si = sp, di = fc
        mov      LEG_QI(%rip),%rdx   
        mov      LEG_NTP1(%rip),%ebp 
        mov      LEG_NLON(%rip),%r9d 
        mov      LEG_NHPP(%rip),%r8d 
        movl     %r9d,%eax           # eax = nlon
        shl      $3,%r9              # nlon * 8
        shr      $1,%eax             # nlon / 2
        sub      %ebp,%eax           # eax = nlon / 2 - ntp1
        xorpd    %xmm0,%xmm0         # remains zero
LOOP_LAT_SP2FC:
        push     %rsi                # save sp(:)
        movl     %ebp,%ecx           # reset m = ntp1
LOOP_M_SP2FC:
        xorpd    %Xmm3,%Xmm3         # prepare sum symm
        xorpd    %Xmm4,%Xmm4         # prepare sum anti
        movl     %ecx,%ebx           # index n = m
LOOP_N_SP2FC:
#       ---------------------------- # symmetric mode
        movddup  (%rdx),%xmm1        # xmm1 = qi(w)
        movupd   (%rsi),%xmm2        # xmm2 = sp(w)
        mulpd    %xmm1,%xmm2         # xmm1 = qi(w) * sp(w)
        addpd    %xmm2,%xmm3         # sum += qi(w) * sp(w)
        add      $8,%rdx             # qi++
        add      $16,%rsi            # sp += 2
        dec      %ebx                # --n
        jz       LOOP_N_SP2FC_EXIT   
#       ---------------------------- # antisymmetric mode
        movddup  (%rdx),%xmm1        # xmm1 = qi(w)
        movupd   (%rsi),%xmm2        # xmm2 = sp(w)
        mulpd    %xmm1,%xmm2         # xmm1 = qi(w) * sp(w)
        addpd    %xmm2,%xmm4         # sum += qi(w) * sp(w)
        add      $8,%rdx             # qi++
        add      $16,%rsi            # sp += 2
        dec      %ebx                # --n
        jnz      LOOP_N_SP2FC        
LOOP_N_SP2FC_EXIT:
        movapd   %xmm3,%xmm1         # fs
        addpd    %xmm4,%xmm1         # fs + fa
        subpd    %xmm4,%xmm3         # fs - fa
        movupd   %xmm1,(%rdi)        # fc(m     ,l) = fs + fa
        movupd   %Xmm3,(%rdi,%r9)    # fc(m+nlon,l) = fs - fa
        add      $16,%rdi            # fc += 2
        loop     LOOP_M_SP2FC        # while (m != 0)
        mov      %eax,%ecx           # set remaining fc's to zero
LOOP_CLEAR_SP2FC:
        movupd   %xmm0,(%rdi)        # fc north = 0.0
        movupd   %xmm0,(%rdi,%r9)    # fc south = 0.0
        add      $16,%rdi            # sizeof(complex)
        loop     LOOP_CLEAR_SP2FC    
        add      %r9,%rdi            # skip southern lat
        pop      %rsi                # reset sp(:)
        dec      %r8d                # --l
        jnz      LOOP_LAT_SP2FC      # while (l != 0)
        pop      %rbp                
        pop      %rbx                
        ret                          

# Fast Indirect Legendre Transformation
# =====================================
# E. Kirk - 17-Nov-2011
# =====================================
# call sp2fcdmu(sp,fc)
#              (di,si)

.globl _sp2fcdmu_
_sp2fcdmu_:
.globl sp2fcdmu_
sp2fcdmu_:
        push     %rbx                
        push     %rbp                
        xchg     %rdi,%rsi           # si = sp, di = fc
        mov      LEG_QJ(%rip),%rdx   
        mov      LEG_NTP1(%rip),%ebp 
        mov      LEG_NLON(%rip),%r9d 
        mov      LEG_NHPP(%rip),%r8d 
        movl     %r9d,%eax           # eax = nlon
        shl      $3,%r9              # nlon * 8
        shr      $1,%eax             # nlon / 2
        sub      %ebp,%eax           # eax = nlon / 2 - ntp1
        xorpd    %xmm0,%xmm0         # remains zero
LOOP_LAT_SP2FCDMU:
        push     %rsi                # save sp(:)
        movl     %ebp,%ecx           # reset m = ntp1
LOOP_M_SP2FCDMU:
        xorpd    %Xmm3,%Xmm3         # prepare sum symm
        xorpd    %Xmm4,%Xmm4         # prepare sum anti
        movl     %ecx,%ebx           # index n = m
LOOP_N_SP2FCDMU:
#       ---------------------------- # symmetric mode
        movddup  (%rdx),%xmm1        # xmm1 = qi(w)
        movupd   (%rsi),%xmm2        # xmm2 = sp(w)
        mulpd    %xmm1,%xmm2         # xmm1 = qi(w) * sp(w)
        addpd    %xmm2,%xmm3         # sum += qi(w) * sp(w)
        add      $8,%rdx             # qi++
        add      $16,%rsi            # sp += 2
        dec      %ebx                # --n
        jz       LOOP_N_SP2FCDMU_EXIT 
#       ---------------------------- # antisymmetric mode
        movddup  (%rdx),%xmm1        # xmm1 = qi(w)
        movupd   (%rsi),%xmm2        # xmm2 = sp(w)
        mulpd    %xmm1,%xmm2         # xmm1 = qi(w) * sp(w)
        addpd    %xmm2,%xmm4         # sum += qi(w) * sp(w)
        add      $8,%rdx             # qi++
        add      $16,%rsi            # sp += 2
        dec      %ebx                # --n
        jnz      LOOP_N_SP2FCDMU     
LOOP_N_SP2FCDMU_EXIT:
        movapd   %xmm4,%xmm1         # fa
        addpd    %xmm3,%xmm1         # fs + fa
        subpd    %xmm3,%xmm4         # fa - fs
        movupd   %xmm1,(%rdi)        # fc(m     ,l) = fa + fs
        movupd   %Xmm4,(%rdi,%r9)    # fc(m+nlon,l) = fa - fs
        add      $16,%rdi            # fc += 2
        loop     LOOP_M_SP2FCDMU     # while (m != 0)
        mov      %eax,%ecx           # set remaining fc's to zero
LOOP_CLEAR_SP2FCDMU:
        movupd   %xmm0,(%rdi)        # fc north = 0.0
        movupd   %xmm0,(%rdi,%r9)    # fc south = 0.0
        add      $16,%rdi            # sizeof(complex)
        loop     LOOP_CLEAR_SP2FCDMU 
        add      %r9,%rdi            # skip southern lat
        pop      %rsi                # reset sp(:)
        dec      %r8d                # --l
        jnz      LOOP_LAT_SP2FCDMU   
        pop      %rbp                
        pop      %rbx                
        ret                          

# Fast Legendre Transformation
# ============================
# E. Kirk - 17-Nov-2011
# ============================
# call dv2uv(pd,pz,pu,pv,db)
# register: (di,si,dx,cx,r8)
#
# rdi: pd(nesp)
# rsi: pz(nesp)
# r12: pu(nlon,nlpp)
# r13: pv(nlon,nlpp)
# r14: qu(ncsp,nhpp)
# r15: qv(ncsp,nhpp)
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
        push     %r15                
        push     %r14                
        push     %r13                
        push     %r12                
        push     %r11                
        push     %r10                
        push     %rbx                # callee saved
        push     %rbp                # callee saved
        mov      %rdx,%r12           # r12 = pu
        mov      %rcx,%r13           # r13 = pv

        mov      LEG_QU(%rip),%r14   
        mov      LEG_QV(%rip),%r15   
        mov      LEG_NTP1(%rip),%ebp 
        mov      LEG_NHPP(%rip),%r8d 
        mov      LEG_NLON(%rip),%r9d 
        movsd    LEG_PVOR(%rip),%xmm2 # xmm2 = plavor

        movsd    16(%rsi),%xmm1      # xmm1 = pz(3)
        movsd    %xmm1,%xmm3         # xmm3 = pz(3)
        subsd    %xmm2,%xmm3         # xmm3 = pz(3) - plavor
        movsd    %xmm3,16(%rsi)      # pz(3)= pz(3) - plavor

        mov      %r9d,%eax           # eax = nlon
        shl      $3,%r9              # nlon * 8
        shr      $1,%eax             # nlon / 2
        sub      %ebp,%eax           # eax = nlon / 2 - ntp1
        xorpd    %xmm0,%xmm0         # remains zero
LOOP_LAT_DV2UV:                      # do r8d = nhpp , 1 , -1
        xor      %r10,%r10           # index = 0
        mov      %ebp,%ecx           
LOOP_M_DV2UV:                        # do ecx = ntp1 , 1 , -1
        xorpd    %Xmm2,%Xmm2         # symm(ud) = 0
        xorpd    %Xmm3,%Xmm3         # symm(vd) = 0
        xorpd    %Xmm4,%Xmm4         # symm(uz) = 0
        xorpd    %Xmm5,%Xmm5         # symm(vz) = 0
        xorpd    %Xmm6,%Xmm6         # anti(ud) = 0
        xorpd    %Xmm7,%Xmm7         # anti(vd) = 0
        xorpd    %Xmm8,%Xmm8         # anti(uz) = 0
        xorpd    %Xmm9,%Xmm9         # anti(vz) = 0
        mov      %ecx,%ebx           # index n = m

LOOP_N_DV2UV:                        # do ebx = ecx  , 1 , -1
#       -------- symm -------------- #
        movupd   (%rdi,%r10),%xmm12  # d(:)
        movupd   (%rsi,%r10),%xmm13  # z(:)
        movddup  (%r14),%xmm14       # xmm14 = qu
        movddup  (%r15),%xmm15       # xmm15 = qv
        add      $8,%r14             # qu += 1
        add      $8,%r15             # qv += 1
        add      $16,%r10            # size of (complex)
#       ---------------------------- #
        movapd   %xmm14,%xmm10       # u
        mulpd    %xmm12,%xmm10       # ud
        addsubpd %xmm10,%xmm2        # sum(ud) +si -sr
        mulpd    %xmm15,%xmm12       # vd
        subpd    %xmm12,%xmm3        # sum(-vd)
        mulpd    %xmm13,%xmm14       # uz
        addsubpd %xmm14,%xmm4        # sum(uz) +si -sr
        mulpd    %xmm13,%xmm15       # vz
        addpd    %xmm15,%xmm5        # sum(vz)
#       ---------------------------- #
        dec      %ebx                # --n
        jz       L_N_EXIT_DV2UV      
#       -------- anti -------------- #
        movupd   (%rdi,%r10),%xmm12  # d(:)
        movupd   (%rsi,%r10),%xmm13  # z(:)
        movddup  (%r14),%xmm14       # xmm14 = qu
        movddup  (%r15),%xmm15       # xmm15 = qv
        add      $8,%r14             # qu += 1
        add      $8,%r15             # qv += 1
        add      $16,%r10            # size of (complex)
#       ---------------------------- #
        movapd   %xmm14,%xmm10       # u
        mulpd    %xmm12,%xmm10       # ud
        addsubpd %xmm10,%xmm6        # sum(ud) +si -sr
        mulpd    %xmm15,%xmm12       # vd
        subpd    %xmm12,%xmm7        # sum(-vd)
        mulpd    %xmm13,%xmm14       # uz
        addsubpd %xmm14,%xmm8        # sum(uz) +si -sr
        mulpd    %xmm13,%xmm15       # vz
        addpd    %xmm15,%xmm9        # sum(vz)
#       ---------------------------- #
        dec      %ebx                # --n
        jnz      LOOP_N_DV2UV        
L_N_EXIT_DV2UV:
#       ---------------------------- #
        shufpd   $1,%xmm2,%xmm2      # swap ud(s)
        shufpd   $1,%xmm4,%xmm4      # swap uz(s)
        shufpd   $1,%xmm6,%xmm6      # swap ud(a)
        shufpd   $1,%xmm8,%xmm8      # swap uz(a)
#       ---------------------------- #
        movapd   %xmm5,%xmm10        # un = vzs
        addpd    %xmm2,%xmm10        # un = vzs + uds
        addpd    %xmm9,%xmm10        # un = vzs + uds + vza
        addpd    %xmm6,%xmm10        # un = vzs + uds + vza + uda
#       ---------------------------- #
        movapd   %xmm4,%xmm11        # vn = uzs
        addpd    %xmm3,%xmm11        # vn = uzs - vds
        addpd    %xmm8,%xmm11        # vn = uzs - vds + uza
        addpd    %xmm7,%xmm11        # vn = uzs - vds + uza - vda
#       ---------------------------- #
        subpd    %xmm5,%xmm2         # us = uds - vzs
        addpd    %xmm9,%xmm2         # us = uds - vzs + vza
        subpd    %xmm6,%xmm2         # us = uds - vzs + vza - uda
#       ---------------------------- #
        subpd    %xmm3,%xmm4         # vs = uzs + vds
        subpd    %xmm8,%xmm4         # vs = uzs + vds - uza
        addpd    %xmm7,%xmm4         # vs = uzs + vds - uza - vda
#       ---------------------------- #
        movupd   %xmm10,(%r12)       # pu(m     ,l)
        movupd   %xmm11,(%r13)       # pv(m     ,l)
        movupd   %xmm2,(%r12,%r9)    # pu(m+nlon,l)
        movupd   %xmm4,(%r13,%r9)    # pv(m+nlon,l)
        add      $16,%r12            # pu += 2
        add      $16,%r13            # pv += 2
        dec      %ecx                
        jnz      LOOP_M_DV2UV        # while (m != 0)
        mov      %eax,%ecx           # set remaining fc's to zero
LOOP_CLEAR_DV2UV:
        movupd   %xmm0,(%r12)        
        movupd   %xmm0,(%r13)        
        movupd   %xmm0,(%r12,%r9)    
        movupd   %xmm0,(%r13,%r9)    
        add      $16,%r12            # pu += 2
        add      $16,%r13            # pv += 2
        loop     LOOP_CLEAR_DV2UV    
        add      %r9,%r12            # skip southern lat
        add      %r9,%r13            # skip southern lat
        dec      %r8d                # --l
        jnz      LOOP_LAT_DV2UV      # while (l != 0)
        movsd    %xmm1,16(%rsi)      # restore pz(3)
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
        push     %rbp                # frame pointer
        mov      %rsp,%rbp           
        sub      $48,%rsp            # 6 local pointer
        push     %r15                
        push     %r14                
        push     %r13                
        push     %r12                
        push     %r11                
        push     %r10                
        push     %rbx                

        mov      %rdi,-48(%rbp)      #  d(:)
        mov      %rsi,-40(%rbp)      #  t(:)
        mov      %rdx,-32(%rbp)      #  z(:)
        mov      %rcx,-24(%rbp)      # tn(:)
        mov      %r8,-16(%rbp)       # fu(:)
        mov      %r9,-8(%rbp)        # fv(:)

        mov      LEG_QQ(%rip),%r14   
        mov      LEG_QC(%rip),%r13   
        mov      LEG_QE(%rip),%r12   
        mov      LEG_QM(%rip),%r11   
        mov      LEG_SPM1(%rip),%ecx 
        mov      LEG_NHPP(%rip),%r8d 
        mov      LEG_NLON(%rip),%r9d 
        mov      LEG_NTP1(%rip),%r15d 
        shl      $1,%ecx             
        xorps    %xmm0,%xmm0         # xmm0 = 0
L_INIT_MKTEND:
        movups   %xmm0,(%rdi,%rcx,8) # d(:) = 0.0
        movups   %xmm0,(%rsi,%rcx,8) # t(:) = 0.0
        movups   %xmm0,(%rdx,%rcx,8) # z(:) = 0.0
        sub      $2,%rcx             
        jns      L_INIT_MKTEND       
        movl     %r9d,%eax           # eax = nlon
        shlq     $3,%r9              # r9  =nlon * 8
        sub      %r15d,%eax          # eax = nlon - ntp1
        shl      $4,%eax             # eax = 16 * (nlon-ntp1)
LOOP_LAT_MKTEND:
        mov      -48(%rbp),%rdi      # reset d(:)
        mov      -40(%rbp),%rsi      # reset t(:)
        mov      -32(%rbp),%rdx      # reset z(:)
        movl     %r15d,%ecx          # reset m = ntp1
LOOP_M_MKTEND:
#       ---------------------------- #
        mov      16(%rbp),%r10       # ke(:)
        movupd   (%r10),%xmm2        # xmm2 = ke: North
        movupd   (%r10,%r9),%xmm1    # xmm1 = ke: South
        movapd   %xmm2,%xmm3         # xmm3 = ke: North
        addpd    %xmm1,%xmm2         # xmm2 = ke: North + South
        subpd    %xmm1,%xmm3         # xmm3 = ke: North - South
        addq     $16,16(%rbp)        # ke += 2
#       ---------------------------- #
        mov      -8(%rbp),%r10       # fv(:)
        movupd   (%r10),%xmm4        # xmm4 = fv: North
        movupd   (%r10,%r9),%xmm1    # xmm1 = fv: South
        movapd   %xmm4,%xmm5         # xmm5 = fv: North
        addpd    %xmm1,%xmm4         # xmm4 = fv: North + South
        subpd    %xmm1,%xmm5         # xmm5 = fv: North - South
        addq     $16,-8(%rbp)        # fv += 2
#       ---------------------------- #
        mov      -16(%rbp),%r10      # fu(:)
        movupd   (%r10),%xmm6        # xmm6 = fu: North
        movupd   (%r10,%r9),%xmm1    # xmm1 = fu: South
        movapd   %xmm6,%xmm7         # xmm7 = fu: North
        addpd    %xmm1,%xmm6         # xmm6 = fu: North + South
        subpd    %xmm1,%xmm7         # xmm7 = fu: North - South
        addq     $16,-16(%rbp)       # fu += 2
#       ---------------------------- #
        mov      24(%rbp),%r10       # ut(:)
        movupd   (%r10),%xmm1        # xmm1 = ut: North
        movupd   (%r10,%r9),%xmm9    # xmm9 = ut: South
        shufpd   $1,%xmm1,%xmm1      # real <-> imag
        shufpd   $1,%xmm9,%xmm9      # real <-> imag
        xorpd    %xmm8,%xmm8         # xmm8 = 0.0
        subpd    %xmm1,%xmm8         # xmm8 = ut: -North
        subpd    %xmm9,%xmm8         # xmm8 = ut: -North -South
        subpd    %xmm1,%xmm9         # xmm9 = ut: -North +South
        addq     $16,24(%rbp)        # ut += 2
#       ---------------------------- #
        mov      32(%rbp),%r10       # vt(:)
        movupd   (%r10),%xmm10       # xmm10 = fu: North
        movupd   (%r10,%r9),%xmm1    # xmm1  = fu: South
        movapd   %xmm10,%xmm11       # xmm11 = fu: North
        addpd    %xmm1,%xmm10        # xmm10 = fu: North + South
        subpd    %xmm1,%xmm11        # xmm11 = fu: North - South
        addq     $16,32(%rbp)        # vt += 2
#       ---------------------------- #
        mov      -24(%rbp),%r10      # tn(:)
        movupd   (%r10),%xmm12       # xmm12 = fu: North
        movupd   (%r10,%r9),%xmm1    # xmm1  = fu: South
        movapd   %xmm12,%xmm13       # xmm13 = fu: North
        addpd    %xmm1,%xmm12        # xmm12 = fu: North + South
        subpd    %xmm1,%xmm13        # xmm13 = fu: North - South
        addq     $16,-24(%rbp)       # tn += 2
#       ---------------------------- #
        mov      %ecx,%ebx           
#       ---------------------------- #
LOOP_N_MKTEND:
#       ------- divergence --------- #
        movupd   (%rdi),%xmm15       # d(w)
#       ------- qq * ke ------------ #
        movddup  (%r14),%xmm1        # xmm1 = qq(w)
        mulpd    %xmm2,%xmm1         # qq(w) * ke(symm)
        addpd    %xmm1,%xmm15        # d(w) += qq(w) * ke(s)
#       ------- qe * fv------------- #
        movddup  (%r12),%xmm1        # xmm1 = qe(w)
        mulpd    %xmm5,%xmm1         # qe(w) * fv(anti)
        subpd    %xmm1,%xmm15        # d(w) -= qe(w) * fv(a)
#       ------- qm * fu ------------ #
        movddup  (%r11),%xmm1        # xmm1 = qm(w)
        mulpd    %xmm6,%xmm1         # qm(w) * fu(symm)
        shufpd   $1,%xmm1,%xmm1      # real <-> imag
        addsubpd %xmm1,%xmm15        # d(w) +- qm(w) * fu(s)
        movupd   %xmm15,(%rdi)       # store d(w)
#       -------- vorticity --------- #
        movupd   (%rdx),%xmm15       # z(w)
#       ------- qe * fu------------- #
        movddup  (%r12),%xmm1        # xmm1 = qe(w)
        mulpd    %xmm7,%xmm1         # qe(w) * fu(anti)
        addpd    %xmm1,%xmm15        # z(w) += qe(w) * fu(a)
#       ------- qm * fv ------------ #
        movddup  (%r11),%xmm1        # xmm1 = qm(w)
        mulpd    %xmm4,%xmm1         # qm(w) * fv(symm)
        shufpd   $1,%xmm1,%xmm1      # real <-> imag
        addsubpd %xmm1,%xmm15        # z(w) +- qm(w) * fv(s)
        movupd   %xmm15,(%rdx)       # store z(w)
#       -------- temperature ------- #
        movupd   (%rsi),%xmm15       # t(w)
#       ------- qe * vt------------- #
        movddup  (%r12),%xmm1        # xmm1 = qe(w)
        mulpd    %xmm11,%xmm1        # qe(w) * vt(anti)
        addpd    %xmm1,%xmm15        # t(w) += qe(w) * vt(a)
#       ------- qc * tn------------- #
        movddup  (%r13),%xmm1        # xmm1 = qc(w)
        mulpd    %xmm12,%xmm1        # qc(w) * tn(symm)
        addpd    %xmm1,%xmm15        # t(w) += qc(w) * tn(s)
#       ------- qm * ut ------------ #
        movddup  (%r11),%xmm1        # xmm1 = qm(w)
        mulpd    %xmm8,%xmm1         # qm(w) * ut(symm)
        addsubpd %xmm1,%xmm15        # t(w) +- qm(w) * ut(s)
        movupd   %xmm15,(%rsi)       # store t(w)
#       ---------------------------- #
        add      $16,%rdi            # d += 2
        add      $16,%rdx            # z += 2
        add      $16,%rsi            # t += 2
        add      $8,%r14             # qq++
        add      $8,%r13             # qc++
        add      $8,%r12             # qe++
        add      $8,%r11             # qm++
        dec      %ebx                # --n
        jz       L_DIV_EXIT          
#       ------- divergence --------- #
        movupd   (%rdi),%xmm15       # d(w)
#       ------- qq * ke ------------ #
        movddup  (%r14),%xmm1        # xmm1 = qq(w)
        mulpd    %xmm3,%xmm1         # qq(w) * ke(anti)
        addpd    %xmm1,%xmm15        # d(w) += qq(w) * ke(a)
#       ------- qe * fv------------- #
        movddup  (%r12),%xmm1        # xmm1 = qe(w)
        mulpd    %xmm4,%xmm1         # qe(w) * fv(symm)
        subpd    %xmm1,%xmm15        # d(w) -= qe(w) * fv(s)
#       ------- qm * fu ------------ #
        movddup  (%r11),%xmm1        # xmm4 = qm(w)
        mulpd    %xmm7,%xmm1         # qm(w) * fu(anti)
        shufpd   $1,%xmm1,%xmm1      # real <-> imag
        addsubpd %xmm1,%xmm15        # d(w) +- qm(w) * fu(a)
        movupd   %xmm15,(%rdi)       # store d(w)
#       -------- vorticity --------- #
        movupd   (%rdx),%xmm15       # z(w)
#       ------- qe * fu------------- #
        movddup  (%r12),%xmm1        # xmm1 = qe(w)
        mulpd    %xmm6,%xmm1         # qe(w) * fu(symm)
        addpd    %xmm1,%xmm15        # z(w) += qe(w) * fu(s)
#       ------- qm * fv ------------ #
        movddup  (%r11),%xmm1        # xmm1 = qm(w)
        mulpd    %xmm5,%xmm1         # qm(w) * fv(anti)
        shufpd   $1,%xmm1,%xmm1      # real <-> imag
        addsubpd %xmm1,%xmm15        # z(w) +- qm(w) * fv(a)
        movupd   %xmm15,(%rdx)       # store z(w)
#       -------- temperature ------- #
        movupd   (%rsi),%xmm15       # t(w)
#       ------- qe * vt------------- #
        movddup  (%r12),%xmm1        # xmm1 = qe(w)
        mulpd    %xmm10,%xmm1        # qe(w) * vt(symm)
        addpd    %xmm1,%xmm15        # t(w) += qe(w) * vt(s)
#       ------- qc * tn------------- #
        movddup  (%r13),%xmm1        # xmm1 = qc(w)
        mulpd    %xmm13,%xmm1        # qc(w) * tn(anti)
        addpd    %xmm1,%xmm15        # t(w) += qc(w) * tn(a)
#       ------- qm * ut ------------ #
        movddup  (%r11),%xmm1        # xmm1 = qm(w)
        mulpd    %xmm9,%xmm1         # qm(w) * ut(anti)
        addsubpd %xmm1,%xmm15        # t(w) +- qm(w) * ut(a)
        movupd   %xmm15,(%rsi)       # store t(w)
#       ---------------------------- #
        add      $16,%rdi            # d += 2
        add      $16,%rdx            # z += 2
        add      $16,%rsi            # t += 2
        add      $8,%r14             # qq++
        add      $8,%r13             # qc++
        add      $8,%r12             # qe++
        add      $8,%r11             # qm++
        dec      %ebx                # --n
        jnz      LOOP_N_MKTEND       
#       ---------------------------- #
L_DIV_EXIT:
        dec      %ecx                
        jnz      LOOP_M_MKTEND       # while (m != 0)
LOOP_M_MKTEND_EXIT:
        add      %rax,-24(%rbp)      # tn += 2 * (nlon-ntp1)
        add      %rax,-16(%rbp)      # fu += 2 * (nlon-ntp1)
        add      %rax,-8(%rbp)       # fv += 2 * (nlon-ntp1)
        add      %rax,16(%rbp)       # ke += 2 * (nlon-ntp1)
        add      %rax,24(%rbp)       # ut += 2 * (nlon-ntp1)
        add      %rax,32(%rbp)       # vt += 2 * (nlon-ntp1)
        dec      %r8d                # --l
        jnz      LOOP_LAT_MKTEND     # while (l != 0)
        pop      %rbx                # restore register
        pop      %r10                
        pop      %r11                
        pop      %r12                
        pop      %r13                
        pop      %r14                
        pop      %r15                
        leave                        
        ret                          

# Fast Direct Legendre Transformation
# =======================================
# E. Kirk - 13-Dec-2011
# =======================================
# call uv2dv(u ,v ,d ,z )
#           (di,si,dx,cx)
#
# arg name  par (bp)  xreg
# -------------------------
#   1   u   rdi
#   2   v   rsi
#   3   d   rdx
#   4   z   rcx
#
# ecx: m loop index
# eax: 2 * (nlon-ntp1) * sizeof(float)
# ebx: n loop index
# r15: ntp1
#
.globl _uv2dv_
_uv2dv_:
.globl uv2dv_
uv2dv_:
        push     %r15                
        push     %r14                
        push     %r13                
        push     %r12                
        push     %r11                
        push     %r10                
        push     %rbx                

        mov      LEG_QE(%rip),%r12   
        mov      LEG_QM(%rip),%r11   
        mov      LEG_NHPP(%rip),%r8d 
        mov      LEG_NLON(%rip),%r9d 
        mov      LEG_NTP1(%rip),%r15d 
        mov      LEG_SPM1(%rip),%ebx 
        shl      $1,%ebx             
        xorpd    %xmm0,%xmm0         # xmm0 = 0
L_INIT_UV2DV:
        movups   %xmm0,(%rdx,%rbx,8) # d(:) = 0.0
        movups   %xmm0,(%rcx,%rbx,8) # z(:) = 0.0
        sub      $2,%rbx             
        jns      L_INIT_UV2DV        
        mov      %r9d,%eax           # eax = nlon
        shlq     $3,%r9              # r9  =nlon * 8
        sub      %r15d,%eax          # eax = nlon - ntp1
        shl      $4,%eax             # eax = 16 * (nlon-ntp1)
LOOP_LAT_UV2DV:
        push     %rdx                # save d(:)
        push     %rcx                # save z(:)
        movl     %r15d,%r14d         # reset m = ntp1
LOOP_M_UV2DV:
#       ------- u ------------------ #
        movupd   (%rdi),%xmm6        # xmm6 = u: North
        movupd   (%rdi,%r9),%xmm1    # xmm1 = u: South
        movapd   %xmm6,%xmm7         # xmm7 = u: North
        addpd    %xmm1,%xmm6         # xmm6 = u: North + South
        subpd    %xmm1,%xmm7         # xmm7 = u: North - South
        add      $16,%rdi            # u += 2
#       ------- v ------------------ #
        movupd   (%rsi),%xmm4        # xmm4 = v: North
        movupd   (%rsi,%r9),%xmm1    # xmm1 = v: South
        movapd   %xmm4,%xmm5         # xmm5 = v: North
        addpd    %xmm1,%xmm4         # xmm4 = v: North + South
        subpd    %xmm1,%xmm5         # xmm5 = v: North - South
        add      $16,%rsi            # v += 2
#       ---------------------------- #
        mov      %r14d,%ebx          
#       ---------------------------- #
LOOP_N_UV2DV:
#       ------- symmetric ---------- #
        movupd   (%rdx),%xmm8        # d(w)
        movupd   (%rcx),%xmm9        # z(w)
#       ------- qe * v ------------- #
        movddup  (%r12),%xmm1        # xmm1 = qe(w)
        mulpd    %xmm5,%xmm1         # qe(w) * v(a)
        subpd    %xmm1,%xmm8         # d(w) -= qe(w) * v(a)
#       ------- qm * u ------------- #
        movddup  (%r11),%xmm1        # xmm1 = qm(w)
        mulpd    %xmm6,%xmm1         # qm(w) * u(s)
        shufpd   $1,%xmm1,%xmm1      # real <-> imag
        addsubpd %xmm1,%xmm8         # d(w) +- qm(w) * u(s)
#       ------- qe * u ------------- #
        movddup  (%r12),%xmm1        # xmm1 = qe(w)
        mulpd    %xmm7,%xmm1         # qe(w) * u(a)
        addpd    %xmm1,%xmm9         # z(w) += qe(w) * u(a)
#       ------- qm * v ------------- #
        movddup  (%r11),%xmm1        # xmm1 = qm(w)
        mulpd    %xmm4,%xmm1         # qm(w) * fv(symm)
        shufpd   $1,%xmm1,%xmm1      # real <-> imag
        addsubpd %xmm1,%xmm9         # z(w) +- qm(w) * v(s)
#       ---------------------------- #
        movupd   %xmm8,(%rdx)        # store d(w)
        movupd   %xmm9,(%rcx)        # store z(w)
#       ---------------------------- #
        add      $16,%rdx            # d += 2
        add      $16,%rcx            # z += 2
        add      $8,%r12             # qe++
        add      $8,%r11             # qm++
        dec      %ebx                # --n
        jz       L_UV2DV_SINGLE      
#       ------- antisymmetric ------ #
        movupd   (%rdx),%xmm8        # d(w)
        movupd   (%rcx),%xmm9        # z(w)
#       ------- qe * v ------------- #
        movddup  (%r12),%xmm1        # xmm1 = qe(w)
        mulpd    %xmm4,%xmm1         # qe(w) * v(s)
        subpd    %xmm1,%xmm8         # d(w) -= qe(w) * v(s)
#       ------- qm * u ------------- #
        movddup  (%r11),%xmm1        # xmm4 = qm(w)
        mulpd    %xmm7,%xmm1         # qm(w) * u(a)
        shufpd   $1,%xmm1,%xmm1      # real <-> imag
        addsubpd %xmm1,%xmm8         # d(w) +- qm(w) * u(a)
#       ------- qe * u ------------- #
        movddup  (%r12),%xmm1        # xmm1 = qe(w)
        mulpd    %xmm6,%xmm1         # qe(w) * u(s)
        addpd    %xmm1,%xmm9         # z(w) += qe(w) * u(s)
#       ------- qm * v ------------- #
        movddup  (%r11),%xmm1        # xmm1 = qm(w)
        mulpd    %xmm5,%xmm1         # qm(w) * v(a)
        shufpd   $1,%xmm1,%xmm1      # real <-> imag
        addsubpd %xmm1,%xmm9         # z(w) +- qm(w) * v(a)
#       ---------------------------- #
        movupd   %xmm8,(%rdx)        # store d(w)
        movupd   %xmm9,(%rcx)        # store z(w)
#       ---------------------------- #
        add      $16,%rdx            # d += 2
        add      $16,%rcx            # z += 2
        add      $8,%r12             # qe++
        add      $8,%r11             # qm++
        dec      %ebx                # --n
        jnz      LOOP_N_UV2DV        
#       ---------------------------- #
L_UV2DV_SINGLE:
        dec      %r14d               
        jnz      LOOP_M_UV2DV        # while (m != 0)
LOOP_M_UV2DV_EXIT:
        add      %rax,%rdi           # u += 2 * (nlon-ntp1)
        add      %rax,%rsi           # v += 2 * (nlon-ntp1)
        pop      %rcx                # restore z(:)
        pop      %rdx                # restore d(:)
        dec      %r8d                # --l
        jnz      LOOP_LAT_UV2DV      # while (l != 0)
        pop      %rbx                # restore register
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
# call tend(q ,qn,uq,vq)
#          (di,si,dx,cx)
#
# arg name  par (bp)  xreg
# -------------------------
#   1   q   rdi
#   2   qn  rsi
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
        push     %r15                
        push     %r14                
        push     %r13                
        push     %r12                
        push     %r11                
        push     %r10                
        push     %rbx                

        mov      %rcx,%r14           # vq(:)

        mov      LEG_QC(%rip),%r13   
        mov      LEG_QE(%rip),%r12   
        mov      LEG_QM(%rip),%r11   
        mov      LEG_SPM1(%rip),%ecx 
        mov      LEG_NHPP(%rip),%r8d 
        mov      LEG_NLON(%rip),%r9d 
        mov      LEG_NTP1(%rip),%r15d 
        shl      $1,%ecx             
        xorps    %xmm0,%xmm0         # xmm0 = 0
L_INIT_QTEND:
        movups   %xmm0,(%rdi,%rcx,8) # q(:) = 0.0
        sub      $2,%rcx             
        jns      L_INIT_QTEND        
        movl     %r9d,%eax           # eax = nlon
        shlq     $3,%r9              # r9  =nlon * 8
        sub      %r15d,%eax          # eax = nlon - ntp1
        shl      $4,%eax             # eax = 16 * (nlon-ntp1)
        xor      %r10,%r10           # r10 = 0
LOOP_LAT_QTEND:
        push     %rdi                # save q(:)
        movl     %r15d,%ecx          # reset m = ntp1
LOOP_M_QTEND:
#       -------- uq ---------------- #
        movupd   (%rdx),%xmm1        # xmm1 = uq: North
        movupd   (%rdx,%r9),%xmm9    # xmm9 = uq: Souqh
        shufpd   $1,%xmm1,%xmm1      # real <-> imag
        shufpd   $1,%xmm9,%xmm9      # real <-> imag
        xorpd    %xmm8,%xmm8         # xmm8 = 0.0
        subpd    %xmm1,%xmm8         # xmm8 = uq: -North
        subpd    %xmm9,%xmm8         # xmm8 = uq: -North -Souqh
        subpd    %xmm1,%xmm9         # xmm9 = uq: -North +Souqh
        add      $16,%rdx            # uq += 2
#       -------- vq ---------------- #
        movupd   (%r14),%xmm10       # xmm10 = fu: North
        movupd   (%r14,%r9),%xmm1    # xmm1  = fu: Souqh
        movapd   %xmm10,%xmm11       # xmm11 = fu: North
        addpd    %xmm1,%xmm10        # xmm10 = fu: North + Souqh
        subpd    %xmm1,%xmm11        # xmm11 = fu: North - Souqh
        add      $16,%r14            # vq += 2
#       ------- qn ----------------- #
        movupd   (%rsi),%xmm12       # xmm12 = fu: North
        movupd   (%rsi,%r9),%xmm1    # xmm1  = fu: Souqh
        movapd   %xmm12,%xmm13       # xmm13 = fu: North
        addpd    %xmm1,%xmm12        # xmm12 = fu: North + Souqh
        subpd    %xmm1,%xmm13        # xmm13 = fu: North - Souqh
        add      $16,%rsi            # qn += 2
#       ---------------------------- #
        mov      %ecx,%ebx           
#       ---------------------------- #
LOOP_N_QTEND:
#       -------- humidity ---------- #
        movupd   (%rdi),%xmm15       # q(w)
#       ------- qe * vt------------- #
        movddup  (%r12,%r10,8),%xmm1 # xmm1 = qe(w)
        mulpd    %xmm11,%xmm1        # qe(w) * vt(anti)
        addpd    %xmm1,%xmm15        # q(w) += qe(w) * vt(a)
#       ------- qc * tn------------- #
        movddup  (%r13,%r10,8),%xmm1 # xmm1 = qc(w)
        mulpd    %xmm12,%xmm1        # qc(w) * tn(symm)
        addpd    %xmm1,%xmm15        # q(w) += qc(w) * tn(s)
#       ------- qm * ut ------------ #
        movddup  (%r11,%r10,8),%xmm1 # xmm1 = qm(w)
        mulpd    %xmm8,%xmm1         # qm(w) * ut(symm)
        addsubpd %xmm1,%xmm15        # q(w) +- qm(w) * ut(s)
        movupd   %xmm15,(%rdi)       # store q(w)
#       ---------------------------- #
        add      $16,%rdi            # q += 2
        inc      %r10                # qx index
        dec      %ebx                # --n
        jz       L_QTEND_EXIT        
#       -------- humidity ---------- #
        movupd   (%rdi),%xmm15       # q(w)
#       ------- qe * vt------------- #
        movddup  (%r12,%r10,8),%xmm1 # xmm1 = qe(w)
        mulpd    %xmm10,%xmm1        # qe(w) * vt(symm)
        addpd    %xmm1,%xmm15        # q(w) += qe(w) * vt(s)
#       ------- qc * tn------------- #
        movddup  (%r13,%r10,8),%xmm1 # xmm1 = qc(w)
        mulpd    %xmm13,%xmm1        # qc(w) * tn(anti)
        addpd    %xmm1,%xmm15        # q(w) += qc(w) * tn(a)
#       ------- qm * ut ------------ #
        movddup  (%r11,%r10,8),%xmm1 # xmm1 = qm(w)
        mulpd    %xmm9,%xmm1         # qm(w) * ut(anti)
        addsubpd %xmm1,%xmm15        # q(w) +- qm(w) * ut(a)
        movupd   %xmm15,(%rdi)       # store q(w)
#       ---------------------------- #
        add      $16,%rdi            # q += 2
        inc      %r10                # qx index
        dec      %ebx                # --n
        jnz      LOOP_N_QTEND        
#       ---------------------------- #
L_QTEND_EXIT:
        dec      %ecx                
        jnz      LOOP_M_QTEND        # while (m != 0)
LOOP_M_QTEND_EXIT:
        add      %rax,%rsi           # qn += 2 * (nlon-ntp1)
        add      %rax,%rdx           # uq += 2 * (nlon-ntp1)
        add      %rax,%r14           # vq += 2 * (nlon-ntp1)
        pop      %rdi                # reset q(:)
        dec      %r8d                # --l
        jnz      LOOP_LAT_QTEND      # while (l != 0)
        pop      %rbx                # restore register
        pop      %r10                
        pop      %r11                
        pop      %r12                
        pop      %r13                
        pop      %r14                
        pop      %r15                
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

        mov      (%rdi),%eax         
        mov      %eax,LEG_NTP1(%rip) # truncation + 1
        mov      (%rsi),%eax         
        mov      %eax,LEG_NLON(%rip) # longitudes
        mov      (%rdx),%eax         
        mov      %eax,LEG_NHPP(%rip) # longitudes
        mov      (%rcx),%eax         
        mov      %eax,LEG_SPM1(%rip) # spectral modes - 1
        mov      (%r8),%rax          
        mov      %rax,LEG_PVOR(%rip) # planetary vorticity
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

        mov      %rdi,LEG_QI(%rip)   
        mov      %rsi,LEG_QJ(%rip)   
        mov      %rdx,LEG_QC(%rip)   
        mov      %rcx,LEG_QE(%rip)   
        ret                          

.globl _legpolb_
_legpolb_:
.globl legpolb_
legpolb_:

        mov      %rdi,LEG_QM(%rip)   
        mov      %rsi,LEG_QQ(%rip)   
        mov      %rdx,LEG_QU(%rip)   
        mov      %rcx,LEG_QV(%rip)   
        ret                          

        .data                        
        .balign   8                   
LEG_QI: .quad 0
LEG_QJ: .quad 0
LEG_QC: .quad 0
LEG_QE: .quad 0
LEG_QM: .quad 0
LEG_QQ: .quad 0
LEG_QU: .quad 0
LEG_QV: .quad 0

LEG_PVOR:
        .quad    0                   
LEG_NTP1:
        .long    0                   
LEG_NLON:
        .long    0                   
LEG_NHPP:
        .long    0                   
LEG_SPM1:
        .long    0                   
