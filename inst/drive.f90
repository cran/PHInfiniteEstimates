      program drive
implicit none
      integer nr,nc
      integer,dimension(:),allocatable:: nvec,resp
      real,dimension(:,:),allocatable:: dm,xx
      logical,dimension(:),allocatable:: conditionon
      integer(kind=8),allocatable,dimension(:)::newcnt
      integer ii,jj,nret,nm
      character(len=20) fmt
      logical verbose
      verbose=.false.
      nr=8
      nc=4
      allocate(nvec(nr),resp(nr),dm(nr,nc),conditionon(nc))
      do jj=1,nc-1
         conditionon(jj)=.true.
      end do
      conditionon(nc)=.false.
      open(30,file="input")
      fmt='( 4(f4.1,x),i2,x,i2)'
      do ii=1,nr
         read(30,*) (dm(ii,jj),jj=1,nc),nvec(ii),resp(ii)
         write(6,fmt) (dm(ii,jj),jj=1,nc),nvec(ii),resp(ii)
      end do
      close(30)
      nret=1
      allocate(xx(nret,nc),newcnt(nret))
      call network(nr,nc,dm,nvec,resp,conditionon,verbose,nret,nm,xx,newcnt)
!     write(6,*) "Finished first network run; nret",nret
      deallocate(xx,newcnt)
!     write(6,*) "Finished network deallocation; nret",nret
      allocate(xx(nret,nc),newcnt(nret))
      call network(nr,nc,dm,nvec,resp,conditionon,verbose,nret,nm,xx,newcnt)
      deallocate(nvec,resp,dm,conditionon)
      do ii=1,nret
         write(6,*) (xx(ii,jj),jj=1,nm),newcnt(ii)
      end do 
      deallocate(xx,newcnt)
      end
