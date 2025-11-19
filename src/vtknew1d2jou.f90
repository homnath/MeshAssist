!------------------------------------------------------------------------------
!> @file vtknew1d2jou.f90
!> @brief Converts NEW VTK 1D file to CUBIT/Trelis journal file.
!>
!>  This program reads ASCII NEW vtk files with unstructured grid of lines and
!>  points only, and removes the redundant lines. The redundant nodes can be
!>  removed within the paraview itself using the 'Clean to Grid' filter.
!>
!> <!-- @athor Hom Nath Gharti (hgharti_AT_princeton_DOT_edu -->
!>
!> ## Compile:
!>  gfortran vtknew1d2jou.f90 -o vtknew1d2jou
!>
!> ## Usage:
!>  vtknew1d2jou \em input_file \n\n
!>  vtknew1d2jou vtknew1d2jou_example.vtk
!>
!------------------------------------------------------------------------------
include 'stringmanip.f90'
program vtknew1d2jou
use strings
implicit none
integer,parameter :: kreal=selected_real_kind(12)
type nodal
  integer :: nelmt
  integer :: elmt(0:4) ! I have considered maximum of 5 elements
  real(kind=kreal) :: x,y,z
end type nodal
type(nodal),allocatable :: node(:)

type elemental
  integer :: nod1,nod2
  integer :: id
end type elemental
type(elemental),allocatable :: elmt(:)

integer :: i,i1,i2,j,narg,np,next_node,next_elmt,pnode,pelmt,e
integer :: temp_nc,nc,dumi,temp_n1,temp_n2,numc,nump,nis,node_is(0:1000),       &
curve(0:30000),ncurve
integer,allocatable :: n1(:),n2(:),new_n1(:),new_n2(:),cid(:),nb(:),nmir(:),  &
op2np(:)
real(kind=kreal),allocatable :: xyz(:,:)
character(len=10) :: dumc
character(len=180) :: outjou_fname,inp_fname,outvtk_fname
character(len=80) :: file_head,ext,path
character(len=256) :: line
integer :: ios
integer,allocatable :: connect(:,:)
! Input and initialisation
narg=command_argument_count()
if (narg <= 0) then
  print*,'ERROR: no input file!'
  stop
endif
! get input file name
call get_command_argument(1, inp_fname)

open(unit=10,file=trim(inp_fname),action='read')
do
  read(10,'(a)',iostat=ios)line
  print*,first_word(line)
  if(first_word(line).eq.'POINTS')then
    read(line,*)dumc,np,dumc
    exit
  endif
enddo
print*,'Number of points:',np
allocate(xyz(0:2,0:np-1),nb(0:np-1),op2np(0:np-1),nmir(0:np-1))

read(10,*)xyz

do
  read(10,'(a)',iostat=ios)line
  if(first_word(line).eq.'CELLS')then
    read(line,*)dumc,temp_nc,dumi
    exit
  endif
enddo
temp_nc=temp_nc-1
print*,'Number of cells:',temp_nc
allocate(n1(0:temp_nc-1),n2(0:temp_nc-1),cid(0:temp_nc-1))

do
  read(10,'(a)',iostat=ios)line
  if(first_word(line).eq.'CONNECTIVITY')exit
enddo
allocate(connect(0:1,0:temp_nc-1))
read(10,*)connect

print*,connect(:,0)
print*,connect(:,temp_nc-2)
! Remove cell with same node
nc=0
do i=0,temp_nc-1
  temp_n1=connect(0,i)
  temp_n2=connect(1,i)
  if (temp_n1/=temp_n2) then
    n1(nc)=temp_n1
    n2(nc)=temp_n2
    nc=nc+1
  end if
end do
close(10)
deallocate(connect)
! Remove double curves
cid=1 ! Initially all cid are 1
numc=nc
do i=0,nc-1
  do j=0,nc-1
    if(i/=j .and. cid(i)/=0 .and. cid(j)/=0)then
      if((n1(i)==n1(j) .and. n2(i)==n2(j)) .or. (n1(i)==n2(j) .and.           &
        n2(i)==n1(j)))then
        cid(j)=0
        numc=numc-1
      end if
    end if
  end do
end do

allocate(new_n1(0:numc-1),new_n2(0:numc-1))
numc=0
do i=0,nc-1
  if(cid(i)/=0)then
    new_n1(numc)=n1(i)
    new_n2(numc)=n2(i)
    numc=numc+1
  end if
end do

cid=1 ! Now all cid become 1

! Number of neighbours of the points
nb=0
do i=0,numc-1
  nb(new_n1(i))=nb(new_n1(i))+1
  nb(new_n2(i))=nb(new_n2(i))+1
end do

! Remove free line elements
nc=numc
do i=0,numc-1
  if(nb(new_n1(i))==1 .or. nb(new_n2(i))==1)then
    cid(i)=0
    nc=nc-1
  end if
end do
deallocate(n1,n2)

allocate(n1(0:nc-1),n2(0:nc-1))
! Renumbering cells
nc=0
do i=0,numc-1
  if(cid(i)/=0)then
    n1(nc)=new_n1(i)
    n2(nc)=new_n2(i)
    nc=nc+1
  end if
end do
deallocate(new_n1,new_n2)

! New number of neibours
nb=0
do i=0,nc-1
  nb(n1(i))=nb(n1(i))+1
  nb(n2(i))=nb(n2(i))+1
end do

! Remove free nodes and renumber
nmir=-1
op2np=-1
nump=0
do i=0,np-1
  if(nb(i)/=0)then
    nmir(nump)=i
    op2np(i)=nump
    nump=nump+1
  end if
end do

print*,'Number of new cells:',nc

! Final nodes and connectivity
allocate(node(0:nump-1),elmt(0:nc-1))
do i=0,nump-1
  node(i)%x=xyz(0,nmir(i))
  node(i)%y=xyz(1,nmir(i))
  node(i)%z=xyz(2,nmir(i))

  node(i)%nelmt=0 ! Initialize
end do
deallocate(xyz)

do i=0,nc-1
  elmt(i)%nod1=op2np(n1(i))
  elmt(i)%nod2=op2np(n2(i))
  elmt(i)%id=0
end do
deallocate(n1,n2)

! Processing for the creation of cubit journal file
! Neigbours of nodes in new system
do i=0,nc-1
  i1=elmt(i)%nod1
  i2=elmt(i)%nod2

  node(i1)%elmt(node(i1)%nelmt)=i
  node(i2)%elmt(node(i2)%nelmt)=i

  node(i1)%nelmt=node(i1)%nelmt+1
  node(i2)%nelmt=node(i2)%nelmt+1
end do

! Count intersecting nodes
nis=0
do i=0,nump-1
  if(node(i)%nelmt>2)then
    node_is(nis)=i
    nis=nis+1
  end if
end do

print*,'Number of intersecting nodes: ',nis

call parse_file(inp_fname,path,file_head,ext)
outvtk_fname=trim(file_head)//'_new.vtk'
outjou_fname=trim(file_head)//'.jou'

open(unit=10,file=trim(outvtk_fname),action='write')
open(unit=20,file=trim(outjou_fname),action='write')

write(10,'(a)')'# vtk DataFile Version 3.0'
write(10,'(a)')'vtk output'
write(10,'(a)')'ASCII'
write(10,'(a)')'DATASET UNSTRUCTURED_GRID'
write(10,'(a,i5,a)')'POINTS ',nump,' float'

do i=0,nump-1
  write(10,'(3(f16.6,1x))')node(i)%x,node(i)%y,node(i)%z
  write(20,'(a,3(f16.6,1x))')'create vertex ',node(i)%x,node(i)%y,node(i)%z
end do
write(10,'(a,1x,i5,1x,i5)')'CELLS',nc,3*nc

do i=0,nc-1
  !if(cid(i)/=0)then
    write(10,'(i2,1x,i5,1x,i5)')2,elmt(i)%nod1,elmt(i)%nod2
  !end if
end do
write(10,'(a,1x,i5)')'CELL_TYPES',nc
write(10,'(i5)')(3,i=0,nc-1)
close(10)

! Create splines for crossing curves
int_node: do i=0,nis-1
  !write(*,'(5(i4,1x))')node_is(i),elmt(node(node_is(i))%elmt(0:2))%id
  int_elmt: do j=0,node(node_is(i))%nelmt-1

    ncurve=0
    curve=-1
    e=node(node_is(i))%elmt(j)

    ! This element is already included in the curve
    if(elmt(e)%id==1)cycle int_elmt
    curve(ncurve)=node_is(i)
    ncurve=ncurve+1

    pnode=node_is(i)
    pelmt=e
    inf_loop: do
      ! Next node
      if(elmt(pelmt)%nod1==pnode)then
        next_node=elmt(pelmt)%nod2
      else if(elmt(pelmt)%nod2==pnode)then
        next_node=elmt(pelmt)%nod1
      else
        print*,'ERROR: Node and element mismatch!',pnode,elmt(pelmt)%nod1,    &
        elmt(pelmt)%nod2
        stop
      end if

      curve(ncurve)=next_node
      ncurve=ncurve+1

      elmt(pelmt)%id=1 ! Now it is included in the curve

      ! Meet another intersection point or starting point
      if(node(next_node)%nelmt>2 .or. next_node==node_is(i))exit inf_loop

      pnode=next_node

      ! Next element
      if(node(next_node)%elmt(0)==pelmt)then
        next_elmt=node(next_node)%elmt(1)
      else if(node(next_node)%elmt(1)==pelmt)then
        next_elmt=node(next_node)%elmt(0)
      else
        print*,'ERROR: illegal next element for a node!'
        stop
      end if
      pelmt=next_elmt

      ! print*,'IntNode: ',pnode,pelmt
      !if(ncurve==6)exit
    end do inf_loop
    !print*,node(24)%elmt(0),node(24)%elmt(1)
    write(20,*)'create curve spline vertex ',curve(0:ncurve-1)+1
  end do int_elmt
end do int_node

! Create curves for sigle loops
do i=0,nc-1
  ncurve=0
  curve=-1

  if(elmt(i)%id==1)cycle

  curve(ncurve)=elmt(i)%nod1
  ncurve=ncurve+1

  pnode=elmt(i)%nod1
  pelmt=i
  do
    ! Next node
    if(elmt(pelmt)%nod1==pnode)then
      next_node=elmt(pelmt)%nod2
    else
      next_node=elmt(pelmt)%nod1
    end if

    curve(ncurve)=next_node
    ncurve=ncurve+1

    elmt(pelmt)%id=1 ! Now it is included in the curve

    ! Meet another intersection point
    if(node(next_node)%nelmt>2 .or. next_node==elmt(i)%nod1)exit

    pnode=next_node

    ! Next element
    if(node(next_node)%elmt(0)==pelmt)then
      next_elmt=node(next_node)%elmt(1)
    else
      next_elmt=node(next_node)%elmt(0)
    end if
    pelmt=next_elmt

    print*,'IntNode: ',ncurve,pnode,pelmt

    if(elmt(next_elmt)%id==1)exit
    !if(ncurve==6)exit
  end do
  write(20,*)'create curve spline vertex ',curve(0:ncurve-1)+1
  print*,i
end do
deallocate(node,elmt)

close(20)
end program vtknew1d2jou
