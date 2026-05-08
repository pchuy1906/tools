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
character(len=5) :: a_str, b_str, c_str

character(len=50) :: XYZfile
integer :: natom, iatom, iatom1, iatom2
real(8) :: box(3,3), invbox(3,3), Vunit
real(8), allocatable :: xyz(:,:)
character(len=5), allocatable :: aname(:), atype(:)
real(8) :: dist
integer :: npair, nframe

real(8) :: xyz_cell(3), xyz1(3), xyz2(3), dxyz(3), tran(3)
integer :: ix, iy, iz

integer :: itype, ntype, id_print



write(*,*) "usage: ./a.out arg1_natom arg2_a arg3_b arg4_c arg5_xyzfile"
write(*,*)

CALL getarg(1, natom_str)
CALL getarg(2, a_str)
CALL getarg(3, b_str)
CALL getarg(4, c_str)
CALL getarg(5, XYZfile)

box = 0.0d0
read(a_str , *) box(1,1)
read(b_str , *) box(2,2)
read(c_str , *) box(3,3)

read(natom_str , *) natom

write(*,*) "The XYZfile is ", XYZfile
write(*,*) "Found the number of atom: ", natom 
write(*,*) "read file ntype.dat for the number of atom types"

Open(unit=44, file="ntype.dat")
read(44,*) ntype
allocate(atype(ntype))
do itype = 1, ntype
  read(44,*) atype(itype)
enddo
close(44)

write(*,*) "Read the XYZ file:"

allocate(xyz(natom,3), aname(natom))

Open(unit=11, file=XYZfile)
Open(unit=22, file="file.gen")

nframe = 0
do while (.true.)
  read(11,*,END=200)
  read(11,*) 
  Vunit = box(1,1)*box(2,2)*box(3,3)
  nframe = nframe + 1
  !
  write(22,*) natom
  write(22,*) atype(:)
  do iatom = 1, natom
    read(11,*) aname(iatom), xyz(iatom,1:3)
    do itype = 1, ntype
      if (aname(iatom)==atype(itype)) id_print = itype
    enddo
    write(22,*) iatom, id_print, xyz(iatom,1:3)
  enddo
  write(22,*) 0.0d0, 0.0d0, 0.0d0
  write(22,*) box(1,:)
  write(22,*) box(2,:)
  write(22,*) box(3,:)
  !
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


Subroutine add_histogram_z(Nz, zarray, z0, denz, scale_f)
integer :: Nz, iz
real(8) :: zarray(Nz), denz(Nz), z0, width, scale_f
real(8) :: zmin, zmax, del
zmin = zarray(1)
zmax = zarray(Nz)
del  = (zmax-zmin)/dble(Nz-1)
iz   = int((z0-zmin)/del) + 1
if ( (iz<0) .or. (iz>Nz) ) then
  denz(1) = denz(1) + 0.0
else
  denz(iz) = denz(iz) + 1.0*scale_f
endif
End subroutine add_histogram_z

