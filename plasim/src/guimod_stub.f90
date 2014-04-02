! Dummy module replacement for guimod.f90
! Use this for environments with no X11

subroutine guistart
return
end subroutine guistart

subroutine guistep_puma
return
end subroutine guistep_puma

subroutine guistep_plasim
return
end subroutine guistep_plasim

subroutine guistop
return
end subroutine guistop

subroutine guihor(yn,f,klev,pm,pa)
character (len=*) :: yn
return
end subroutine guihor

subroutine guigv(yn,f)
character (len=*) :: yn
return
end subroutine guigv

subroutine guigvcol(yname,f,klon)
character (len=*) :: yname

return
end subroutine guigvcol


subroutine guigtcol(f,klon)
return
end subroutine guigtcol


subroutine guid3dcol(yname,f,klon,klev,pm,pa)
character (len=*) :: yname
return
end subroutine guid3dcol

subroutine guips(fp)
return
end subroutine guips

subroutine guiput(yn,f,k1,k2,k3)
character (len=*) :: yn
return
end subroutine guiput

subroutine guihorlsg(yn,f,mask,klev,pm,pa)
character (len=*) :: yn 
return
end subroutine guihorlsg

