Module Fileparser
  implicit none
  Contains
  
  Subroutine bytereader(keyword,rd_len,ifile,i_err)
    implicit none
    character(len=*), intent(out) :: keyword  !< Parsed keyword
    integer, intent(out) :: rd_len !< Length of parsed keyword
    integer, intent(in) :: ifile  !< File to read from
    integer, intent(out) :: i_err  !< Error status of reading
    logical :: rd_done,rd_start
    
    rd_done=.false.
    rd_start=.false.
    rd_len=0
    do while(.not.rd_done.and.rd_len<len(keyword))
      rd_len=rd_len+1
      read(ifile,'(a1)',advance='no',end=20,eor=10) keyword(rd_len:rd_len)
      rd_start=rd_start.or.(keyword(rd_len:rd_len)/=" ".and.keyword(rd_len:rd_len)/=":")
      rd_done=rd_start.and.(keyword(rd_len:rd_len)==" ".or.keyword(rd_len:rd_len)==":")
      if(keyword(rd_len:rd_len)==":") keyword(rd_len:rd_len)=""
    end do
    
    i_err=0
    keyword=adjustl(keyword(1:rd_len)//'')
    return
    ! final word
    10  continue
    i_err=10
    keyword=adjustl(keyword(1:rd_len)//'')
    return
    ! end of file
    20  continue
    i_err=20
    keyword=adjustl(keyword(1:rd_len)//'')
    return
    
  End Subroutine bytereader
  
  
  !> Convert lower case characters to upper case
  Subroutine small2caps(str)
    implicit none
    character(len=*),intent(inout):: str  !< string to convert
    
    integer i
    
    do i=1,len(str)
      if(str(i:i)>="a" .and. str(i:i)<= "z") str(i:i)=achar(iachar(str(i:i))-32)
    end do
    
  End Subroutine small2caps
  
  
  !> Convert upper case characters to lower case
  Subroutine caps2small(str)
    implicit none
    character(len=*),intent(inout):: str  !< string to convert
    
    integer i
    
    do i=1,len(str)
      if(str(i:i)>="A" .and. str(i:i)<= "Z") str(i:i)=achar(iachar(str(i:i))+32)
    end do
    
  End Subroutine caps2small

End Module Fileparser
