!############################################################
      subroutine cntdup(xx,nr,nc,cnt,nu,bnds,sst,conditionon)
implicit none
      integer nr,nc
      integer(kind=8) cnt(nr)
      logical conditionon(nc)
      real xx(nr,nc),bnds(nc,2),sst(nc)
      integer ii,jj,kk,nu
      logical match
!     write(6,*) "Entering cntdup nr,nc",nr,nc
!     do ii=1,2
!        write(6,*) "Bounds",(bnds(jj,ii),jj=1,nc)
!     end do
!     write(6,*) "Sst",(sst(jj),jj=1,nc)
      do ii=1,nr
         do kk=1,nc
            if(conditionon(kk)) then
!              write(6,*) "Checking bounds ii",ii,"kk",kk,"xx(ii,kk)",xx(ii,kk)
               if((xx(ii,kk)+bnds(kk,1)).lt.sst(kk)) then
!                 write(6,*) "Bound check 1"
                  cnt(ii)=0
               endif
               if((xx(ii,kk)-bnds(kk,2)).gt.sst(kk)) then
!                 write(6,*) "Bound check 2"
                  cnt(ii)=0
               endif
            endif
         end do
      end do
      do kk=1,nr-1
         if(cnt(kk).gt.0) then
            do jj=kk+1,nr
               if(cnt(jj).gt.0) then
                  match=.true.
                  do ii=1,nc
                     if(xx(jj,ii).ne.xx(kk,ii)) match=.false.
                  end do!ii
!                 write(6,*) "kk,jj,match",kk,jj,match
                  if(match) then
                     cnt(kk)=cnt(kk)+cnt(jj)
                     cnt(jj)=0
                  end if
               end if
            end do
         end if
      end do
      nu=0
      do kk=1,nr
         if(cnt(kk)>0) nu=nu+1
      end do
!     write(6,*) "About to write xx, nr=",nr," nu=",nu," nc=",nc
!     do ii=1,nr
!        if(cnt(ii).gt.0) write(6,*) "ii",ii,"xx", (xx(ii,jj),jj=1,nc),cnt(ii)
!     end do
      return
      end
