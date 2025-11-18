!-----------------------------------------------------
!> @file xyz2jou.f90
!> @brief Coverts UTM/XYZ file to CUBIT/Trelis journal file.
!>
!>  This program converts a UTM or XYZ file to a CUBIT journal file. The UTM or XYZ
!>  file contains the three columns of X, Y, and Z coordinates, respectively.
!>
!> <!-- @athor Hom Nath Gharti (hgharti_AT_princeton_DOT_edu -->
!>
!> ## Compile:
!>  - in parent folder, type:
!>  make
!>  OR
!>  - in src/ folder, type
!>  gfortran xyz2jou.f90 -o xyz2jou
!>
!> ## Usage:
!>  ./bin/xyz2jou \em input_file [\em Options] \n\n
!>	Example: \n
!>  ./bin/xyz2jou ./input/xyz2jou_example.utm
!>
!> ## Options:
!> - -nx: Use this option if you know the number of points in a line along X axis
!>       . This will speed up the processing. For example, -nx=100. If it is not
!>       defined, nx is automatically determined.
!> - -nskip: Use this option if you want to skip (downsample) certain number of
!>       successive points. This will skip along both X and Y axes. For example,
!>       -nskip=2. [DEFAULT 0].
!-----------------------------------------------------
include 'stringmanip.f90'
program xyz2jou
use strings
implicit none
integer,parameter :: kreal=selected_real_kind(12)
character(len=180) :: cubit_fname,inp_fname,out_fname
character(len=80) :: file_head,ext,strblk,path
integer :: ierr,ip1,ip2,narg,nx,nx0,npx
! total points, total points written
integer :: np,npw
! total curves, total curves written
integer :: nc,ncw
integer :: ipskip,icskip
integer :: i_arg,itmp,istat,nskip
! tolerance for straight-line slope discrepancy
real(kind=kreal),parameter :: eps=0.1_kreal
real(kind=kreal) :: m0,m,x0,x,y0,y,z
logical :: isnx0,setm0,islastc,islastp

print*,'Running...'
!----input and initialisation----
narg=command_argument_count()
if (narg <= 0) then
  print*,'ERROR: no input file!'
  stop
endif
! get input file name
call get_command_argument(1, inp_fname)

isnx0=.false.
nskip=0
if(narg>1)then
  do i_arg=2,narg
    call get_command_argument(i_arg, strblk)
    ! read nx parameter
    call look_int(strblk,'-nx=',itmp,istat)
    if(istat.eq.0)then
      nx0=itmp
      isnx0=.true.
      cycle
    endif

    ! read nskip parameter
    call look_int(strblk,'-nskip=',itmp,istat)
    if(istat.eq.0)then
      nskip=itmp
      cycle
    endif

    ! unrecognized option
    write(*,*)'ERROR: unrecognized option: ',trim(adjustl(strblk))
    stop

  enddo
endif

! display input info
write(*,*)'------------------------------------------'
write(*,*)'input file: ',trim(inp_fname)
if(isnx0)then
  write(*,*)'nx0 defined: YES'
  write(*,*)'nx0: ',nx0
else
  write(*,*)'nx0 defined: NO'
endif
write(*,*)'nskip: ',nskip
write(*,*)'------------------------------------------'

open(10,file=inp_fname,status='old',action='read',iostat=ierr)
if(ierr/=0)then
  print*,'ERROR: cannot open file "',trim(inp_fname),'"!'
  stop
endif
call parse_file(inp_fname,path,file_head,ext)
out_fname=trim(file_head)//'.jou'
open(20,file=out_fname,status='replace',action='write')

read(10,*,iostat=ierr)x,y,z
if(ierr/=0)then
  print*,'ERROR: cannot read file!'
  stop
endif

ip1=1; npx=1
np=1; npw=1
write(20,*)'reset'
write(20,*)'create vertex ',x,y,z
ipskip=0; icskip=0
x0=x; y0=y

nc=0; ncw=0;
setm0=.true.
islastc=.false.
m0=100.0_kreal ! this value is redundant
do
  read(10,*,iostat=ierr)x,y,z
  if(ierr/=0)exit
  np=np+1
  npx=npx+1
  ipskip=ipskip+1
  islastp=.false.
  !write(20,*)'create vertex ',x,y,z
  ! compute slope
  m=(y-y0)/(x-x0)
  if(setm0)then
    !write(20,*)'create vertex ',x,y,z
    m0=m
    setm0=.false.
    npx=1
  else
    if((isnx0.and.npx.eq.nx0) .or. (.not.isnx0.and.abs(m-m0).gt.eps))then
      nx=npx ! point in each curve
      !ip2=np-1
      ip2=npw!-1
      nc=nc+1 ! number of curves
      icskip=icskip+1

      !if(nc.eq.1)allocate(lastpx(2*nx))
      ! do not skip the first curve
      if(nc.eq.1 .or. icskip.gt.nskip)then
        write(20,*)'create curve spline vertex ',ip1,' to ',ip2
        write(20,*)('delete vertex all')
        ncw=ncw+1
        !if(ncw.ge.5)stop
        icskip=0
        islastc=.true.
      endif
      if((isnx0.or.nc>1).and.(nx0.ne.nx))then
        print*,'WARNING: number of Xpoints mismatch!'
        print*,'nx0:',nx0,' nx:',nx
      endif
      x0=x; y0=y
      nx0=nx
      setm0=.true.
      ip1=ip2+1
      if(nc.le.1)print*,'Total points in X direction:',nx0
    endif
    !write(20,*)'create vertex ',x,y,z
  endif
  if((nc.eq.0 .or. icskip.gt.nskip-1) .and. ipskip.gt.nskip)then
    npw=npw+1
    write(20,*)'create vertex ',x,y,z
    ipskip=0
    islastp=.true.
    islastc=.false.
  endif
enddo
close(10)

!! do not skip last point
!if(.not.islastp)then
!  npw=npw+1
!  write(20,*)'create vertex ',x,y,z
!  ipskip=0
!  islastp=.true.
!endif

! last curve
if(.not.islastc)then
  !ip2=np
  ip2=npw
  nc=nc+1 ! number of curves
  ! do not skip last curve
  write(20,*)'create curve spline vertex ',ip1,' to ',ip2
  write(20,*)('delete vertex all')
  ncw=ncw+1
endif
if(nc*nx.ne.np)then
  print*,'WARNING: total number of points mismatch!'
endif
print*,'Total points:',np
print*,'Total points written:',npw
print*,'Total points skipped:',np-npw
print*,'Total curves:',nc
print*,'Total curves written:',ncw
print*,'Total curves skipped:',nc-ncw

write(20,*)'create surface skin curve all'

write(20,*)('delete curve all')
write(20,*)('compress all ')
cubit_fname=trim(file_head)//'_surface.cub'
write(20,*)'save as "',trim(cubit_fname),'" overwrite'
close(20)
print*,'Complete!'
stop
end program xyz2jou
!===============================================================================
