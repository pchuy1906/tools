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
character(len=5) :: atom1, atom2, natom_str, option_str, rcut_str
integer :: option
real(8) :: rcut
character(len=50) :: XYZfile, TMP_STR
integer :: natom, iatom, iatom1, iatom2
real(8) :: box(3,3), invbox(3,3), Vunit
real(8), allocatable :: xyz(:,:)
character(len=5), allocatable :: aname(:)
real(8) :: dist
integer :: nbond, nframe

real(8) :: xyz_cell(3), xyz1(3), xyz2(3), dxyz(3), tran(3)
integer :: ix, iy, iz

write(*,*) "usage: ./a.out C N  1.5 traj.xyz  cell(1-orthor   2-non-orthor)"
write(*,*)

CALL getarg(1, atom1)
CALL getarg(2, atom2)
CALL getarg(3, rcut_str)
read(rcut_str , *) rcut
CALL getarg(4, XYZfile)
CALL getarg(5, option_str)
read(option_str , *) option

write(*,*) "The XYZfile is ", XYZfile
write(*,*) "Calculate RDF between ", atom1, " and ", atom2
write(*,*)
write(*,*) "Read the XYZ file:"

Open(unit=11, file=XYZfile)
read(11,*) natom
write(*,*) "Found the number of atom: ", natom
rewind (11)

allocate(xyz(natom,3), aname(natom))

open(unit=66,file='bond.dat')

nframe = 0
do while (.true.)
    read(11,*,END=200)
    box = 0.0d0
    if (option==1) then
        read(11,*) box(1,1), box(2,2), box(3,3)
        Vunit = box(1,1)*box(2,2)*box(3,3)
    elseif (option==2) then
        read(11,*) box(1,:), box(2,:), box(3,:)
        !!! NEED TO BE CHANGED TO CALCULATE DET of matrix
        Vunit = box(1,1)*box(2,2)*box(3,3)
    elseif (option==3) then
        read(11,*) TMP_str, box(1,:), box(2,:), box(3,:)
        !!! NEED TO BE CHANGED TO CALCULATE DET of matrix
        Vunit = box(1,1)*box(2,2)*box(3,3)
    else
        write(*,*) "NOT IMPLEMENTED YET!"
    endif
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
    nbond = 0
    do iatom1 = 1, natom-1
      do iatom2 = iatom1+1, natom
        if ( ((aname(iatom1) .eq. atom1) .and. (aname(iatom2) .eq. atom2)) .or. &
            &((aname(iatom2) .eq. atom1) .and. (aname(iatom1) .eq. atom2)) ) then
          Call smallest_d(xyz(iatom1,1:3), xyz(iatom2,1:3), xyz_cell, dist, invbox, box)
          if (dist .le. rcut) then
              nbond = nbond + 1
          endif
        endif
      enddo
    enddo
  
    write(66,*) nframe, nbond

enddo

close(66)

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


