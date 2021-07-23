!###########################################################
      subroutine copytox(pp,cnt,conditionon,nr,nc,xx,newcnt,nu,nm)
implicit none
      integer nr,nc,nu,nm
      integer(kind=8) cnt(nr),newcnt(nu)
      real pp(nr,nc),xx(nu,nm)
      integer kk,ii,jj,ll
      integer(kind=8) zero
      logical conditionon(nc)
!     logical test
!     write(6,*) "Entered copytox, nu ",nu,"cnt(nr)",cnt(nr)
!     write(6,*) "In copytox, nu ",nu,"cnt",(cnt(kk),kk=1,nr)
      kk=0
      zero=0
      do ii=1,nr
!        test=cnt(ii).gt.zero
         if(cnt(ii).gt.zero) then
            kk=kk+1
!           write(6,*) "ii,kk,nu",ii,kk,nu,"cnt(ii)",cnt(ii)
            newcnt(kk)=cnt(ii)
            ll=0
            do jj=1,nc
               if(.not.conditionon(jj)) then
                  ll=ll+1
                  xx(kk,ll)=pp(ii,jj)
               end if
            end do
         end if
      end do
      return
      end
