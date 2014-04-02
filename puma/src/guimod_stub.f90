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

subroutine guihor(yn,f,pm,pa)
character (len=*) :: yn
return
end subroutine guihor

subroutine guigv(yn,f)
character (len=*) :: yn
return
end subroutine guigv

subroutine guips(fp)
return
end subroutine guips

subroutine guiput(yn,f,k1,k2,k3)
character (len=*) :: yn
return
end subroutine guiput

subroutine guigt(f)
return
end
