Subroutine ReadCoefficients
  use Constants
  use Parameters
  use FileParser
  implicit none
  character(len=50)   :: keyword, cache 
  integer             :: rd_len,i_err,i,j,i_errb, pos, Ndim, I_temp
  Real(dp)            :: tmpdp
  logical             :: comment

  open(ifileno,file=CoeffFile)
  
  do
    10 continue
    keyword=""
    call bytereader(keyword,rd_len,ifileno,i_errb)
    call caps2small(keyword)
    comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&
    (scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1.or.&
    (scan(trim(keyword),'!')==1))
    
    if (comment) then
      read(ifileno,*) cache
    else
    
      keyword=trim(keyword)
      
  SELECT CASE (keyword)
  CASE ("alat")
      read(ifileno, '(A)', iostat=i_err) cache
      read(ifileno, *, iostat=i_err)((alat(i,j),i=1,3),j=1,3)
      if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
      
  CASE ("coeffdisp")
      read(ifileno, '(A)', iostat=i_err) cache
      pos = scan(cache, ':')
      cache = trim(cache(pos+1:))
      read(cache, *, iostat=i_err) Ndim
      do i = 1, Ndim
      read(ifileno,*,iostat=i_err) CoeffDisp(i)
      end do
      if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
      
  CASE ("coeffjijsm")
      read(ifileno, '(A)', iostat=i_err) cache
      pos = scan(cache, ':')
      cache = trim(cache(pos+1:))
      read(cache, *, iostat=i_err) Ndim
      do i = 1, Ndim
      read(ifileno,*,iostat=i_err) CoeffJijsm(i)
      end do
      if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
      
  CASE ("coeffjijhf")
      read(ifileno, '(A)', iostat=i_err) cache
      pos = scan(cache, ':')
      cache = trim(cache(pos+1:))
      read(cache, *, iostat=i_err) Ndim
      do i = 1, Ndim
      read(ifileno,*,iostat=i_err) CoeffJijhf(i)
      end do
      if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

  CASE ("coeffjijsmhf")
      read(ifileno, '(A)', iostat=i_err) cache
      pos = scan(cache, ':')
      cache = trim(cache(pos+1:))
      read(cache, *, iostat=i_err) Ndim
      do i = 1, Ndim
      read(ifileno,*,iostat=i_err) CoeffJijsmhf(i)
      end do
      if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
      
  CASE ("coeffu")
      read(ifileno, '(A)', iostat=i_err) cache
      pos = scan(cache, ':')
      cache = trim(cache(pos+1:))
      read(cache, *, iostat=i_err) Ndim
      do i = 1, Ndim
      read(ifileno,*,iostat=i_err) Coeffu(i)
      end do
      if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
      
  CASE ("coeffdispu")
      read(ifileno, '(A)', iostat=i_err) cache
      pos = scan(cache, ':')
      cache = trim(cache(pos+1:))
      read(cache, *, iostat=i_err) Ndim
      do i = 1, Ndim
      read(ifileno,*,iostat=i_err) CoeffDispu(i)
      end do
      if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
      
  CASE ("coeffepsdisp")
      read(ifileno, '(A)', iostat=i_err) cache
      pos = scan(cache, ':')
      cache = trim(cache(pos+1:))
      read(cache, *, iostat=i_err) Ndim
      do i = 1, Ndim
      read(ifileno,*,iostat=i_err) CoeffEpsDisp(i)
      end do
      if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
      
  CASE ("coeffepsu")
      read(ifileno, '(A)', iostat=i_err) cache
      pos = scan(cache, ':')
      cache = trim(cache(pos+1:))
      read(cache, *, iostat=i_err) Ndim
      do i = 1, Ndim
      read(ifileno,*,iostat=i_err) CoeffEpsu(i)
      end do
      if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
      
  CASE ("coeffeps")
      read(ifileno, '(A)', iostat=i_err) cache
      pos = scan(cache, ':')
      cache = trim(cache(pos+1:))
      read(cache, *, iostat=i_err) Ndim
      do i = 1, Ndim
      read(ifileno,*,iostat=i_err) CoeffEps(i)
      end do
      if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

  CASE ("dwq")
      read(ifileno, '(A)', iostat=i_err) cache
      read(ifileno,*,iostat=i_err) ((DWq(j,i), j=1,3), i=1,3)
      if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

   CASE ("efield")
       read(ifileno, '(A)', iostat=i_err) cache
       do i = 1, 4
       read(ifileno,*,iostat=i_err) EAmp(i), EPhi(i), GateField(i)
       end do
       if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
 
   CASE ("mass")
       read(ifileno, '(A)', iostat=i_err) cache
       pos = scan(cache, ':')
       cache = trim(cache(pos+1:))
       read(cache, *, iostat=i_err) Ndim
       do i = 1, Ndim
       read(ifileno,*,iostat=i_err) mass(i)
       end do
       do i = 1, Ndim
       mass(i) = mass(i)*mpme
       end do
       if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

   CASE ("nosemass")
       read(ifileno, '(A)', iostat=i_err) cache
       pos = scan(cache, ':')
       cache = trim(cache(pos+1:))
       read(cache, *, iostat=i_err) Ndim
       do i = 1, Ndim
       read(ifileno,*,iostat=i_err) nosemass(i)
       end do
       if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

   CASE ("efermi")
       read(ifileno, '(A)', iostat=i_err) cache
       pos = scan(cache, ':')
       cache = trim(cache(pos+1:))
       read(cache, *, iostat=i_err) Ndim
       Allocate(Efermi(Ndim))
       do i = 1, Ndim
       read(ifileno,*,iostat=i_err) Efermi(i)
       end do
       if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

   CASE ("contour")
       read(ifileno, '(A)', iostat=i_err) cache
       pos = scan(cache, ':')
       cache = trim(cache(pos+1:))
       read(cache, *, iostat=i_err) Ndim
       read(ifileno,*,iostat=i_err) ContourMin, ContourMax, ContourHeight
       read(ifileno,*,iostat=i_err) ContourNPoints(:)
       if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

   CASE ("wanniersite")
       read(ifileno, '(A)', iostat=i_err) cache
       pos = scan(cache, ':')
       cache = trim(cache(pos+1:))
       read(cache, *, iostat=i_err) NumWannSites
       Allocate(SiteOrbInfo(3,NumWannSites))
       do i = 1, NumWannSites
       Read(ifileno,*,iostat=i_err) SiteOrbInfo(1,i), SiteOrbInfo(2,i), SiteOrbInfo(3,i)
       end do
       do i = 1, 2
         do j = 1, NumWannSites
           If(SiteOrbInfo(3,j).gt.0) then
             JijSites(i) = j
             SiteOrbInfo(3,j) = SiteOrbInfo(3,j) - 1
             exit
           end if
         end do
       end do
       if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

   CASE ("poscar")
       read(ifileno, '(A)', iostat=i_err) cache
       pos = scan(cache, ':')
       cache = trim(cache(pos+1:))
       Read(cache, *, iostat=i_err) unitcell%nion
       Allocate(unitcell%site(3,unitcell%nion))
       Read(ifileno, *, iostat=i_err)
       Read(ifileno, *, iostat=i_err) tmpdp
       Read(ifileno, *, iostat=i_err) ((unitcell%latt_a(i,j), j=1,3), i=1,3)
       Read(ifileno, *, iostat=i_err)
       Read(ifileno, *, iostat=i_err)
       Read(ifileno, *, iostat=i_err)
       unitcell%latt_a = tmpdp*unitcell%latt_a
       Read(ifileno, *, iostat=i_err) ((unitcell%site(j,i), j=1,3), i=1,unitcell%nion)
       if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
     
  CASE DEFAULT
      if(len(trim(keyword))>0) then
        read(ifileno,*)
      end if
  END SELECT
    end if
    
    if (i_errb==20) goto 20
    if (i_errb==10) goto 10
    
  end do
  
  20 continue
  close(ifileno)
  
End Subroutine ReadCoefficients
