program cmpmat
integer, parameter :: matrx  = 6018
integer, parameter :: ien    =   72
integer, parameter :: kb     = 4*ien+4
integer, parameter :: km     = kb+1
real (kind=8) ea(kb,matrx),ta(km,matrx),sa(matrx)
real (kind=8) eb(kb,matrx),tb(km,matrx),sb(matrx)
real (kind=8) ed(kb,matrx),td(km,matrx),sd(matrx)
open(76,file='mat76',form='unformatted')
open(77,file='mat77',form='unformatted')

read (76) ea
read (76) ta
read (76) sa

read (77) eb
read (77) tb
read (77) sb

ed = ea - eb
td = ta - tb
sd = sa - sb

print *,'ed',maxval(ed),minval(ed)
print *,'td',maxval(td),minval(td)
print *,'sd',maxval(sd),minval(sd)

stop
end
