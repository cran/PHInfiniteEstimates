      subroutine network(nr,nc,dm,nvec,resp,conditionon,verbose,nret,nm,xx,newcnt)
implicit none
      logical,save::flip
      integer nr,nc,nret
      real dm(nr,nc)
      integer nvec(nr),resp(nr)
      logical conditionon(nr),verbose
      integer,save:: na,nb
      integer ii,jj,nu,nm
!     integer kk
      real,dimension(:),allocatable::sst
      real,dimension(:,:),allocatable,save::bnds
      real xx(nret,nc)
      real,dimension(:,:),allocatable,save::pa,pb
      integer(kind=8),dimension(:),allocatable,save::cnta,cntb
      integer(kind=8) newcnt(nret)
      if(nret.eq.1) then
         allocate(sst(nc),bnds(nc,2))
         nm=0
         do ii=1,nc
            if(.not.conditionon(ii)) nm=nm+1
         end do
!        write(6,*) "dm"
!        do ii=1,nr
!           write(6,*) (dm(ii,jj),jj=1,nc)
!        end do
         do jj=1,nc
            sst(jj)=0.0
            do ii=1,nr
               sst(jj)=sst(jj)+dm(ii,jj)*resp(ii)
            end do 
         end do
         na=1
         nu=1
!        write(6,*) "Allocation A, na",na," nc",nc
         allocate(pa(na,nc),cnta(na))
         do jj=1,nc
            pa(1,jj)=0.0
         end do
         cnta(1)=1
         flip=.false.
         do ii=1,nr
            call dobnds(nr,nc,dm,bnds,ii,nvec)
            if(verbose) write(6,*) "ii",ii, " of ", nr
            if(.not.flip) then
               nb=nu*(nvec(ii)+1)
!              write(6,*) "Allocation B, nb",nb," nc",nc,"flip",flip
               allocate(pb(nb,nc),cntb(nb))
               call expand(pa,pb,cnta,cntb,na,nb,nc,dm,nvec,nr,ii)
!              write(6,*) "Exited expand first branch"
!              write(6,*) "cnta",cnta
!              write(6,*) "Deallocate A"
               deallocate(cnta,pa)
!              write(6,*) "Deallocation done"
!              write(6,*) "Entering cntdup in branch a, nb",nb
               call cntdup(pb,nb,nc,cntb,nu,bnds,sst,conditionon)
!              write(6,*) "New cntb",(cntb(jj),jj=1,nb)
!              do jj=1,nb
!                 if(cntb(jj).gt.0) write(6,*) "pb",(pb(jj,kk),kk=1,nc)
!              end do
            else
               na=nu*(nvec(ii)+1)
!              write(6,*) "Allocation A, na",na," nc",nc,"flip",flip
               allocate(pa(na,nc),cnta(na))
               call expand(pb,pa,cntb,cnta,nb,na,nc,dm,nvec,nr,ii)
!              write(6,*) "Exited expand second branch"
!              write(6,*) "cntb",cntb
!              write(6,*) "Deallocate B"
               deallocate(cntb,pb)
!              write(6,*) "Deallocation done"
!              write(6,*) "Entering cntdup in branch b, na",na
               call cntdup(pa,na,nc,cnta,nu,bnds,sst,conditionon)
!              write(6,*) "New cnta",(cnta(jj),jj=1,na)
!              do jj=1,na
!                 if(cnta(jj).gt.0) write(6,*) "pa",(pa(jj,kk),kk=1,nc)
!              end do
            endif
            flip=.not.flip
!           write(6,*) "flip",flip,"nu",nu
         end do
!        write(6,*) "After the end of the loop,flip,nu,nm",flip,nu,nm
         nret=nu
         deallocate(sst,bnds)
      else
!        write(6,*) "In stage 2, flip,nret,nm",flip,nret,nm
!        allocate(xx(nret,nm),newcnt(nret))
         if(flip) then
            call copytox(pb,cntb,conditionon,nb,nc,xx,newcnt,nret,nm)
            deallocate(pb,cntb)
         end if
         if(.not.flip) then
!           write(6,*) "Before copytox na",na
!           write(6,*) "Before copytox cnta(na)",cnta(na)
            call copytox(pa,cnta,conditionon,na,nc,xx,newcnt,nret,nm)
            deallocate(pa,cnta)
         end if
      end if
      return
      end
