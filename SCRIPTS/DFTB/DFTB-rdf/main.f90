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
character(len=10) :: atom1, atom2, natom_str, a_str, b_str, c_str
character(len=50) :: XYZfile
integer :: natom, iatom, iatom1, iatom2
real(8) :: box(3,3), invbox(3,3), Vunit
real(8), allocatable :: xyz(:,:)
character(len=5), allocatable :: aname(:)
real(8) :: dist
integer :: npair, nframe

real(8) :: xyz_cell(3), xyz1(3), xyz2(3), dxyz(3), tran(3)
integer :: ix, iy, iz

integer, parameter :: nr = 1000
real(8) :: rmin, rmax, r(nr), gr(nr)
integer :: ir
rmin = 0.2d0
rmax =10.0d0
do ir = 1, nr
  r(ir) = rmin + (rmax-rmin)/dble(nr-1)*dble(ir-1)
enddo
gr = 0.0d0

write(*,*) "usage: ./a.out arg1_fe arg2_fe natom HISTORY"
write(*,*)

CALL getarg(1, atom1)
CALL getarg(2, atom2)
CALL getarg(3, natom_str)
CALL getarg(4, a_str)
CALL getarg(5, b_str)
CALL getarg(6, c_str)
CALL getarg(7, XYZfile)

read(natom_str , *) natom

box = 0.0d0
read(a_str , *) box(1,1)
read(b_str , *) box(2,2)
read(c_str , *) box(3,3)

write(*,*) "The XYZfile is ", XYZfile
write(*,*) "Found the number of atom: ", natom 
write(*,*) "Calculate RDF between ", atom1, " and ", atom2
write(*,*)
write(*,*) "Read the XYZ file:"

allocate(xyz(natom,3), aname(natom))

Open(unit=11, file=XYZfile)

nframe = 0
do while (.true.)
  read(11,*,END=200)
  read(11,*) 
  Vunit = box(1,1)*box(2,2)*box(3,3)
  nframe = nframe + 1
  !
  invbox = box
  Call DGETRF( M, N, invbox, LDA, IPIV, INFO )
  CALL DGETRI(N, invbox, N, IPIV, WORK, LWORK, INFO)

  !
  do iatom = 1, natom
    read(11,*) aname(iatom), xyz(iatom,1:3)
  enddo
  !
  npair = 0
  do iatom1 = 1, natom-1
    do iatom2 = iatom1+1, natom
      if ( ((aname(iatom1) .eq. atom1) .and. (aname(iatom2) .eq. atom2)) .or. &
          &((aname(iatom2) .eq. atom1) .and. (aname(iatom1) .eq. atom2)) ) then
        Call smallest_d(xyz(iatom1,1:3), xyz(iatom2,1:3), xyz_cell, dist, invbox, box)
        npair = npair + 1
        do ix = -1,1
          do iy = -1,1
            do iz = -1,1
              tran = dble(ix)*box(1,:) + dble(iy)*box(2,:) + dble(iz)*box(3,:)
              xyz1 = xyz(iatom1,1:3)
              xyz2 = xyz_cell + tran
              dxyz = xyz1 - xyz2
              dist = sqrt(dot_product(dxyz, dxyz))
              Call add_histogram_z(nr, r, dist, gr, 1.0d0)
            enddo
          enddo
        enddo
      endif
    enddo
  enddo
enddo
200 continue

write(*,*) "the total number of frame is ", nframe

open(unit=66,file='gr.dat')
do ir = 1, nr-1
  f1 = 4.0d0/3.0d0*pi*(r(ir+1)**3 - r(ir)**3)
  write(66,*) r(ir), gr(ir)/f1/(dble(npair)/Vunit)/dble(nframe)
enddo
close(66)

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

