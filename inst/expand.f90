!#############################################################
      subroutine expand(pa,pb,cnta,cntb,na,nb,nc,dm,nvec,nobs,ll)
implicit none
      integer na,nb,nc,nobs
      integer(kind=8) cnta(na),cntb(nb),bc
      real pa(na,nc),pb(nb,nc),dm(nobs,nc)
      integer nvec(nobs)
      integer ii,jj,kk,ll,mm,ee
!     write(6,*) "Enter expand, na",na,"ll",ll,"nvec(ll)",nvec(ll),"added row",(dm(ll,jj),jj=1,nc)
!     if(nb.ne.(na*(1+nvec(ll)))) write(6,*) "Error"
!     do ii=1,na
!        write(6,*) "pa", (pa(ii,jj),jj=1,nc),cnta(ii)
!     end do
      mm=0
      do ii=1,na
         if(cnta(ii).gt.0) then
            mm=mm+1
            bc=1
            do jj=0,nvec(ll)
               ee=(mm-1)*(nvec(ll)+1)+jj+1
!              write(6,*) "Writing to row ",ee
               do kk=1,nc
                  pb(ee,kk)=pa(ii,kk)+jj*dm(ll,kk)
               end do 
               cntb(ee)=cnta(ii)*bc
               bc=(bc*(nvec(ll)-jj))/(jj+1)
            end do 
         end if
      end do
!     do ii=1,nb
!        write(6,*) "ii of nb",ii,nb
!        if(cntb(ii).gt.0) write(6,*) "ii of nb",ii,nb,"pb", (pb(ii,jj),jj=1,nc),cntb(ii)
!     end do
!     write(6,*) "About to leave expand"
      return
      end
