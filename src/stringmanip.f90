! parse file name and return path, file head, and extension
subroutine parse_file(fname,path,head,ext)
implicit none
character(len=*),intent(in) :: fname
character(len=*),intent(out) :: path,head,ext
integer :: i,ipath,iext,slen

slen=len(fname)

! set default output
path=""
head=""
ext=""

if(len_trim(fname)==0)return

! find dot position
iext=slen+1
do i=slen,1,-1
  if(fname(i:i)==".")then
    iext=i
    exit
  endif
enddo

! find slash position
ipath=0
do i=iext,1,-1
  if(fname(i:i)=="/" .or. fname(i:i)=="\")then
    ipath=i
    exit
  endif
enddo

! determine path, head, and extension
head=fname(ipath+1:iext-1)
path=fname(1:ipath)
ext=fname(iext+1:slen)
return
end subroutine parse_file
!==============================================================================

subroutine look_int(strblk,vartag,ival,istat)
implicit none
character(len=*),intent(in) :: strblk,vartag
integer,intent(out) :: ival,istat
integer :: ios,ncstr,ncvar
ncstr=len(trim(strblk))
ncvar=len(trim(vartag))
istat=-1
if(vartag(1:ncvar) .eq. strblk(1:ncvar))then
  ! variable found
  read(strblk(ncvar+1:ncstr),*,iostat=ios)ival
  if(ios.ne.0)return
  istat=0
endif
return
end subroutine look_int
!==============================================================================
