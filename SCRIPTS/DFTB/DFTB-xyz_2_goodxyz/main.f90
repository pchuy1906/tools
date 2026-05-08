program main
implicit none

!LAPACK SUBROUTINE
!details search the subroutine
Integer, parameter :: ndim = 3
integer :: M=ndim
integer :: N=ndim
integer :: LDA=ndim
integer :: IPIV(ndim)
integer :: INFO
real(8) :: WORK(ndim)
integer :: LWORK = ndim

real(8), parameter :: pi=3.14159265358979d0
real(8) :: f1
character(len=5) :: atom1, atom2, natom_str
character(len=50) :: XYZfile
integer :: natom, iatom, iatom1, iatom2
real(8) :: box(3,3), invbox(3,3), Vunit
real(8) :: xyz(3), tmpxyz(4)
character(len=5) :: aname
real(8) :: dist
integer :: npair, nframe

real(8) :: xyz_cell(3), xyz1(3), xyz2(3), dxyz(3), tran(3)
integer :: ix, iy, iz, ixyz

integer :: itype, ntype, id_print
character (len=10) :: a_str, b_str, c_str
real(8) :: a, b ,c

write(*,*) "usage: ./a.out arg1_xyzfile a b c"
write(*,*)

CALL getarg(1, XYZfile)
CALL getarg(2, a_str)
CALL getarg(3, b_str)
CALL getarg(4, c_str)

read(a_str , *) a
read(b_str , *) b
read(c_str , *) c

box = 0.d0
box(1,1) = a
box(2,2) = b
box(3,3) = c

invbox = box
Call DGETRF( M, N, invbox, LDA, IPIV, INFO )
CALL DGETRI(N, invbox, N, IPIV, WORK, LWORK, INFO)

write(*,*) "The XYZfile is ", XYZfile
write(*,*) "Read the XYZ file:"

Open(unit=11, file=XYZfile)
Open(unit=22, file="newfile.xyz")

nframe = 0
do while (.true.)
    read(11,*,END=200) natom
    read(11,*) 
    nframe = nframe + 1
    !
    write(22,*) natom
    write(22,*) a, b, c
    do iatom = 1, natom
        read(11,*) aname, xyz(1:3), tmpxyz(1:4)
        xyz = matmul(xyz,invbox)
        do ixyz = 1,3
            xyz(ixyz) = xyz(ixyz)-dble(floor(xyz(ixyz)))
        enddo
        xyz = matmul(xyz,box)
        write(22,'(a4,3f15.6)') aname, xyz(1:3)
    enddo
enddo
200 continue

write(*,*) "the total number of frame is ", nframe

end program


Subroutine smallest_d(O1xyz, O2xyz_ini, O2xyz, dist, invbox, box)
real(8):: O1xyz(3), O2xyz_ini(3), O2xyz(3), invbox(3,3), box(3,3), dist
real(8):: Cxyz(3)
integer:: k4
Cxyz = O1xyz - O2xyz_ini
Cxyz = matmul(Cxyz,invbox)
do k4 = 1,3
  Cxyz(k4) = Cxyz(k4)-dble(nint(Cxyz(k4)))
enddo
Cxyz = matmul(Cxyz,box)
dist = sqrt(dot_product(Cxyz, Cxyz))
O2xyz= O1xyz - Cxyz
End Subroutine smallest_d
