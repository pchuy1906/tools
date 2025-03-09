import numpy as np

import argparse
parser = argparse.ArgumentParser(description='making the larger system')
# Arguments supported by the code.
parser.add_argument("--file_dump", default='dump.lammps', help='file dump in LAMMPS')
parser.add_argument("--Ny", default=20, help='Ny')
parser.add_argument("--Nz", default=80, help='Nz')
parser.add_argument('--idtarget', nargs='+', type=int)
#
args        = parser.parse_args()
file_dump    = args.file_dump
Ny           = args.Ny
Nz           = args.Nz
idtarget     = args.idtarget

f  = open(file_dump ,"rt")
nframe=0
while True:

    # read line "ITEM: TIMESTEP"
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break

    # read step number
    tmp = f.readline().split()
    nstep = int(tmp[0])

    # read line "ITEM: NUMBER OF ATOMS"
    tmp = f.readline()

    # read the number of atoms
    tmp = f.readline().split()
    natom = int(tmp[0])

    # read line "ITEM: BOX BOUNDS pp pp ss"
    tmp = f.readline()

    # read the cell
    tmp = f.readline().split()
    if nframe==0:
        Lx = float(tmp[1])-float(tmp[0])
    tmp = f.readline().split()
    if nframe==0:
        Ly = float(tmp[1])-float(tmp[0])
    tmp = f.readline().split()
    if nframe==0:
        Lz = float(tmp[1])-float(tmp[0])

    # read line "ITEM: ATOMS id type element mass x y z vx vy vz"
    tmp = f.readline()

    xyz  = np.zeros(shape=(natom,3))
    vxyz = np.zeros(shape=(natom,3))
    mass = np.zeros(shape=(natom))
    #asym = np.zeros(shape=(natom))
    asym = np.array([ "C" for k in range(natom)])

    nframe += 1

    for iatom in range(natom):
        tmp = f.readline().split()
        id = int(tmp[0])-1
        asym[id] = tmp[2]
        mass[id] = float(tmp[3])
        TMP_XYZ = [float(tmp[i]) for i in range(4,7)]
        TMP_VXYZ = [float(tmp[i]) for i in range(7,10)]


#        ! read line "141179  1        Cr        52        3.71006 3.74448 5.04227    1.0503 -0.326707 2.75573 -0 0"
#        read(11,*)   id,     TMP_INT, asym(id), mass(id), xyz(id,1:3),               vxyz(id,1:3)
#    enddo
#
#    do itarget = 1, ntarget
#        if (idtarget(itarget) == nstep) then
#            write(*,*) "AAA", nstep
#            call compute_Temp(natom, mass, xyz, vxyz, Lx, Ly, Lz, Ny, Nz, ngrid, Tgrid, mgrid, nstep)
#            write(*,*) "Good"
#        endif
#    enddo
#
#!    if (idtarget(ntarget) < nstep) go to 101
#!
#!    if (nframe==1) then
#!        write(66,*) natom
#!        write(66,'(3f15.8)') Lx, Ly, Lz
#!        do iatom = 1, natom
#!            write(66,'(A,3f15.8)') asym(iatom), xyz(iatom,1:3)
#!        enddo
#!    endif
#!
#enddo
#200 continue
#!101 continue
#!close(66)
#!
print ("the total number of frame is ", nframe)
#
#end program
#
#
#
#Subroutine compute_Temp(natom, mass, xyz, vxyz, Lx, Ly, Lz, Ny, Nz, ngrid, Tgrid, mgrid, nstep)
#real(8), parameter :: lammps2kcalmol = 2390.05736, kb = 0.0019872067, Navo = 0.6022140857
#integer :: natom, Ny, Nz, nstep
#real(8) :: mass(natom), xyz(natom,3), vxyz(natom,3)
#real(8) :: Lx, Ly, Lz
#real(8) :: ngrid(Ny,Nz), Tgrid(Ny,Nz), mgrid(Ny,Nz), vgrid(Ny,Nz)
#real(8) :: T0
#
#real(8) :: dy, dz
#integer :: idy, idz
#real(8) :: v2
#character(len=40) :: str_nstep
#
#integer :: iatom
#
#T0    = 0.d0
#Tgrid = 0.d0
#ngrid = 0.d0
#mgrid(:,:) = 0.d0
#vgrid = 0.d0
#
#dy = Ly/dble(Ny)
#dz = Lz/dble(Nz)
#
#!write(*,*) "Ny=", Ny
#!write(*,*) "Nz=", Nz
#!write(*,*) "Ly=", Ly
#!write(*,*) "Lz=", Lz
#!write(*,*) "dy=", dy
#!write(*,*) "dz=", dz
#
#do iatom = 1, natom
#    idy = int(xyz(iatom,2)/dy)+1
#    idz = int(xyz(iatom,3)/dz)+1
#    !write(*,*) "y/z= ", xyz(iatom,2), xyz(iatom,3)
#    !write(*,*) "iatom/idy/idz= ", iatom, idy, idz
#
#    ngrid(idy,idz) = ngrid(idy,idz) + 1.0
#    v2 = vxyz(iatom,1)*vxyz(iatom,1) + vxyz(iatom,2)*vxyz(iatom,2) + vxyz(iatom,3)*vxyz(iatom,3)
#
#    !Tgrid(idy,idz) = Tgrid(idy,idz) + mass(iatom)*v2
#    mgrid(idy,idz) = mgrid(idy,idz) + mass(iatom)
#    !write(*,*) iatom, idy, idz, mass(iatom)
#    !vgrid(idy,idz) = vgrid(idy,idz) + vxyz(iatom,3)
#    !T0 = T0 + mass(iatom)*v2
#enddo
#
#!write(*,*) "T0 mv2", T0
#!T0 = T0*lammps2kcalmol
#!write(*,*) "T0 mv2 (kcal/mol)", T0
#!T0 = T0/(3.0*kb)
#!write(*,*) "T0 mv2/3*kb (T)", T0
#
#
#
#
#
#!do idy = 1, Ny
#!    do idz = 1, Nz
#!        if (ngrid(idy,idz) > 0.5d0) then
#!            !Tgrid(idy,idz) = (lammps2kcalmol * Tgrid(idy,idz)/ngrid(idy,idz) ) / (3.d0 * kb)
#!            mgrid(idy,idz) = mgrid(idy,idz)/(Navo * Lx*dy*dz)
#!            vgrid(idy,idz) = vgrid(idy,idz)/ngrid(idy,idz)
#!        endif
#!    enddo
#!enddo
#
#write (str_nstep, *) nstep
#str_nstep = adjustl(str_nstep)
#!write(*,*) "contour"// trim(str_nstep)// ".dat"
#!write(*,*) Ny, Nz
#!
#!Open(unit=666, file="contour_T_"// trim(str_nstep)// ".dat")
#!write(666,*) Ny, Nz
#!do idy = 1, Ny
#!    do idz = 1, Nz
#!        write(666,'(3f20.3)') (dble(idy)-1.0d0)*dy, (dble(idz)-1.0d0)*dz, Tgrid(idy,idz)
#!    enddo
#!enddo
#!close(666)
#!
#!
#!Open(unit=777, file="contour_density_"// trim(str_nstep)// ".dat")
#!write(777,*) Ny, Nz
#!do idy = 1, Ny
#!    do idz = 1, Nz
#!        write(777,'(3f15.4)') (dble(idy)-1.0d0)*dy, (dble(idz)-1.0d0)*dz, mgrid(idy,idz)
#!    enddo
#!enddo
#!close(777)
#!
#!
#!Open(unit=888, file="contour_vz_"// trim(str_nstep)// ".dat")
#!write(888,*) Ny, Nz
#!do idy = 1, Ny
#!    do idz = 1, Nz
#!        write(888,'(3f15.4)') (dble(idy)-1.0d0)*dy, (dble(idz)-1.0d0)*dz, vgrid(idy,idz)
#!    enddo
#!enddo
#!close(888)
#!
#!write(*,*) "OK here"
#End Subroutine
#
#
print (asym)
print (mass)
print (TMP_XYZ)
print (TMP_VXYZ)

