c       lmp2gmx.f: Translates a lammps data file to Gromacs .gro
c       and Gromacs .itp forcefield files. User must specify input
c       and output file names in in.lmp2gmx file. 
c
c       Authors: Chris Lorenz chris.lorenz@kcl.ac.uk   &
c                Mohamed Ali al-Badri mohamed.al-badri@kcl.ac.uk
c
c       License: GNU General Public License v3.0


      integer maxatypes,maxbtypes,maxantypes,maxdtypes,maxitypes
      integer maxatoms
      parameter (maxatypes=20,maxbtypes=20,maxantypes=50,maxdtypes=120)
      parameter (maxitypes=25,maxatoms=5000)
      character*5 atype,fftype(maxatypes)
      character*80 fconfig,fgro,fitp,fff
      integer natoms,nbonds,nangles,ndihedrals,nimpropers,natypes
      integer nbtypes,nantypes,ndtypes,nitypes,i,j,j1,j2,j3,j4
      integer atnum(maxatypes)
      double precision xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz,mass(maxatypes)
      double precision pc(2,maxatypes),bc(2,maxbtypes),ac(2,maxantypes)
      double precision dc(4,maxdtypes),dcn(4,maxdtypes),ic(2,maxitypes)
     
      dcn = 0.0d0 
      open(2,file='in.convert_lammps_to_gromacs')
      read(2,*) 
      read(2,*) fconfig
      read(2,*)
      read(2,*) fgro
      read(2,*) 
      read(2,*) fitp
      read(2,*) 
      read(2,*) fff
      close(2)
      
      open(10,file=fconfig)
      open(12,file=fgro)
      open(14,file=fitp)
      open(16,file=fff)
     
      read(10,*) 
      read(10,*) natoms
      read(10,*) nbonds
      read(10,*) nangles
      read(10,*) ndihedrals
      read(10,*) nimpropers
      read(10,*)
      read(10,*) natypes
      read(10,*) nbtypes
      read(10,*) nantypes
      read(10,*) ndtypes
      read(10,*) nitypes
      read(10,*)
      read(10,*) xlo,xhi
      read(10,*) ylo,yhi
      read(10,*) zlo,zhi
      lx = xhi-xlo
      ly = yhi-ylo
      lz = zhi-zlo
      read(10,*) 
      read(10,*)
      read(10,*)
      do i = 1,natypes
         read(10,*) j,mass(j),fftype(j)
      end do
      read(10,*)
      read(10,*)
      read(10,*)
      do i = 1,natypes
         read(10,*) j,pc(1,j),pc(2,j)
      end do
      read(10,*)
      read(10,*)
      read(10,*)
      do i = 1,nbtypes
         read(10,*) j,bc(1,j),bc(2,j)
      end do
      read(10,*) 
      read(10,*)
      read(10,*)
      do i = 1,nantypes
         read(10,*) j,ac(1,j),ac(2,j)
      end do
      read(10,*)
      read(10,*)
      read(10,*) 
      do i = 1,ndtypes
         read(10,*) j,dc(1,j),dc(2,j),dc(3,j),dc(4,j)
      end do
      read(10,*)
      read(10,*)
      read(10,*)
      do i = 1,nitypes
         read(10,*) j,ic(1,j),ic(2,j)
      end do
      read(10,*)
      read(10,*)
      read(10,*)
      write(12,*)
      write(12,'(i6)') natoms
      write(14,'(a16)') '[ moleculetype ]'
      write(14,'(a4, i1)') 'GRAP',3
      write(14,*) 
      write(14,'(a9)') '[ atoms ]'
      do i = 1,natoms
         read(10,*) j,jm,jt,qj,xj,yj,zj
         
         if (int(mass(jt)) .eq. 1) then
            atype = '    H'
            atnum(jt) = 1
         else if (int(mass(jt)) .eq. 12) then
            atype = '    C'
            atnum(jt) = 6
         else if (int(mass(jt)) .eq. 15) then
            atype = '    O'
            atnum(jt) = 8
         else if (int(mass(jt)) .eq. 14) then
            atype = '    N'
            atnum(jt) = 7
         else 
            write(6,*) 'Unknown element: ', mass(jt)
            write(6,*) 'Add to conditional statement and recompile.'
            stop
         end if
         write(14,'(i5,1x,a5,i5,1x,a4,a5,i5,2f10.6)') j,fftype(jt),jm,
     &    'GRAP',atype,j,qj,mass(jt)
         write(12,'(i5,2a5,i5,3f8.3)') jm,'GRAP',atype,j,(xj/10.0d0),
     &    (yj/10.0d0),(zj/10.0d0)
      end do
      write(14,*)
      write(16,'(a12)') '[ defaults ]'
      write(16,'(2(i1,1x),a3,2(1x,f3.1))') 1,3,'yes',0.5,0.5
      write(16,*)
      write(16,'(a13)') '[ atomtypes ]'
      write(16,*)
      do it = 1,natypes
         write(16,'(a5,i5,2f10.6,a3,2(1x,e13.7))') fftype(it),atnum(it)
     &      ,mass(it),0.0d0,' A ',pc(2,it)/10.0d0,pc(1,it)*418.40d0
      end do
      write(16,'(a75)') 'opls_111   OW  8      9.95140    -0.834       A
     &    3.15061e-01  6.36386e-01'
      write(16,'(a75)') 'opls_112   HW  1      4.03200     0.417       A
     &    0.00000e+00  0.00001e+00'
      write(12,'(3f10.5)') lx/10.0d0,ly/10.0d0,lz/10.0d0 
      read(10,*)
      write(14,*)
      read(10,*)
      write(14,'(a9)') '[ bonds ]'
      read(10,*)
      do i = 1,nbonds
         read(10,*) j,jt,j1,j2
c        check that the '1' is the right number for the functional form
         write(14,'(2(i6,1x),i3,2(1x,e13.7))') j1,j2,1,bc(2,jt)/10.0d0,
     &      bc(1,jt)*418.40d0
      end do
      read(10,*)
      write(14,*)
      read(10,*)
      write(14,'(a10)') '[ angles ]'
      read(10,*)
      write(14,*)
      do i = 1,nangles
         read(10,*) j,jt,j1,j2,j3
c        check that '1' is the right number for the functional form
         write(14,'(3(i6,1x),i3,2(1x,e13.7))') j1,j2,j3,1,ac(2,jt),
     &      ac(1,jt)*418.40d0
      end do
      read(10,*)
      write(14,*)
      read(10,*)
      write(14,'(a13)') '[ dihedrals ]'
      read(10,*)
      write(14,*)
      do i = 1,ndihedrals
          read(10,*) j,jt,j1,j2,j3,j4
c         input correct functional form for the opls dihedrals
c         Also I assume that the values for the different force
c         constants are what found in the file for gromacs, check
c         this --->
c
c        Changed to R-B type dihedral:
         dcn(1,jt)= dc(1,jt) + (0.50d0*dc(2,jt)) + dc(3,jt)
         dcn(2,jt)= (-0.50d0*dc(2,jt))
         dcn(3,jt)= (-1.0d0*dc(3,jt))
         dcn(4,jt)= (-2.0d0*dc(4,jt))
         if (j1 .ne. j2 .and. j1 .ne. j3 .and. j1 .ne. j4) then
            write(14,'(4(i6,1x),i3,6(1x,e13.7))') j1,j2,j3,j4,3,
     &        dcn(1,jt)*418.40d0,dcn(2,jt)*418.40d0,dcn(3,jt)*418.40d0,
     &        dcn(4,jt)*418.40d0,0.0d0,0.0d0
         end if
      end do
      read(10,*)
      write(14,*)
      read(10,*)
      write(14,'(a13)') '[ impropers ]'
      read(10,*)
      write(14,*)
      do i = 1,nimpropers
         read(10,*) j,jt,j1,j2,j3,j4
         if (j1 .ne. j2 .and. j1 .ne. j3 .and. j1 .ne. j4) then
           write(14,'(4(i6,1x),i3,2(1x,e13.7))') j1,j2,j3,j4,2,
     &          ic(2,jt),ic(1,jt)*418.40d0
         end if
      end do
      close(10)
      close(14)
      close(12)
      close(16)
      end
