      subroutine dobnds(nr,nc,dm,bnds,cur,nvec)
implicit none
      integer nr,nc,cur
      integer nvec(nr)
      real dm(nr,nc), bnds(nc,2)
!     real mx(2)
      integer jj,kk
!     write(6,*) "In dobnds, cur",cur
      do kk=1,nc
         bnds(kk,1)=0.0
         bnds(kk,2)=0.0
         if(cur.lt.nr) then
!           mx(1)=dm(cur+1,kk)
!           mx(2)=dm(cur+1,kk)
!           write(6,*) "mx",mx,"kk",kk
!           do jj=cur+1,nr
!              mx(1)=max(mx(1),dm(jj,kk))
!              mx(2)=min(mx(2),dm(jj,kk))
!           end do
            do jj=cur+1,nr
! Column 1 of bnds is max, column 2 is min.
               bnds(kk,2)=bnds(kk,2)+min(dm(jj,kk),0.0)*nvec(jj)
               bnds(kk,1)=bnds(kk,1)+max(dm(jj,kk),0.0)*nvec(jj)
            end do
         end if
      end do
      return
      end
