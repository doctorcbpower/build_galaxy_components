module nrutils_modules
  integer(kind=4) :: iseed
contains
  ! Calculation of integrals...
  SUBROUTINE qromb(func,a,b,ss)
    implicit none
    INTEGER JMAX,JMAXP,K,KM
    REAL*8 a,b,func,ss,EPS
    EXTERNAL func
    PARAMETER (EPS=1.e-6, JMAX=30, JMAXP=JMAX+1, K=5, KM=K-1)
    INTEGER j
    REAL*8 dss,h(JMAXP),s(JMAXP)
    h(1)=1.
    do 11 j=1,JMAX
       call trapzd(func,a,b,s(j),j)
       if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
       endif
       s(j+1)=s(j)
       h(j+1)=0.25*h(j)
11  continue
    pause 'too many steps in qromb'
  END SUBROUTINE qromb

  SUBROUTINE qromo(func,a,b,ss,choose)
    implicit none    
    INTEGER JMAX,JMAXP,K,KM
    REAL*8 a,b,func,ss,EPS
    EXTERNAL func,choose
    PARAMETER (EPS=1.e-6, JMAX=30, JMAXP=JMAX+1, K=5, KM=K-1)
    INTEGER j
    REAL*8 dss,h(JMAXP),s(JMAXP)
    h(1)=1.
    do 11 j=1,JMAX
       call choose(func,a,b,s(j),j)
       if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
       endif
       s(j+1)=s(j)
       h(j+1)=h(j)/9.
11  continue
    pause 'too many steps in qromo'
  END SUBROUTINE qromo
  
  SUBROUTINE qqromb(func,a,b,c,ss)
    implicit none
    INTEGER JMAX,JMAXP,K,KM
    REAL*8 a,b,c,func,ss,EPS
    EXTERNAL func
    PARAMETER (EPS=1.d-6, JMAX=30, JMAXP=JMAX+1, K=5, KM=K-1)
    INTEGER j
    REAL*8 dss,h(JMAXP),s(JMAXP)
    
    h(1)=1.
    do 11 j=1,JMAX
       call trapzd3(func,a,b,c,s(j),j)
       if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
       endif
       s(j+1)=s(j)
       h(j+1)=0.25*h(j)
11  continue
    pause 'too many steps in qromb'
  END SUBROUTINE qqromb
      
  SUBROUTINE mqromb(func,a,b,xp,yp,np,ss,mode)
    implicit none
    integer mode
    integer np
    real*8 xp(np),yp(np)
    INTEGER JMAX,JMAXP,K,KM
    REAL*8 a,b,func,ss,EPS
    EXTERNAL func
    PARAMETER (EPS=1.d-6, JMAX=30, JMAXP=JMAX+1, K=5, KM=K-1)
    INTEGER j
    REAL*8 dss,h(JMAXP),s(JMAXP)

    h(1)=1.
    do 11 j=1,JMAX
       if(mode.eq.0) then
          call trapzd2(func,a,b,xp,yp,np,s(j),j)
       else
          call midpnt2(func,a,b,xp,yp,np,s(j),j)       
       end if
       if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
       endif
       s(j+1)=s(j)
       if(mode.eq.0) then       
          h(j+1)=0.25*h(j)
       else
          h(j+1)=h(j)/9
       endif
11  continue
    pause 'too many steps in qromb'
  END SUBROUTINE mqromb

  SUBROUTINE polint(xa,ya,n,x,y,dy)
    implicit none
    INTEGER n,NMAX
    REAL*8 dy,x,y,xa(n),ya(n)
    PARAMETER (NMAX=10)
    INTEGER i,m,ns
    REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
    ns=1
    dif=abs(x-xa(1))
    do 11 i=1,n
       dift=abs(x-xa(i))
       if (dift.lt.dif) then
          ns=i
          dif=dift
       endif
       c(i)=ya(i)
       d(i)=ya(i)
11  continue
    y=ya(ns)
    ns=ns-1
    do 13 m=1,n-1
       do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12     continue
       if (2*ns.lt.n-m)then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       endif
       y=y+dy
13  continue
    return
  END SUBROUTINE polint

  SUBROUTINE trapzd(func,a,b,s,n)
    implicit none
    INTEGER n
    REAL*8 a,b,s,func
    EXTERNAL func
    INTEGER it,j
    REAL*8 del,sum,tnm,x
    if (n.eq.1) then
       s=0.5*(b-a)*(func(a)+func(b))
    else
       it=2**(n-2)
       tnm=it
       del=(b-a)/tnm
       x=a+0.5*del
       sum=0.
       do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11     continue
       s=0.5*(s+(b-a)*sum/tnm)
    endif
    return
  END SUBROUTINE trapzd

  SUBROUTINE trapzd2(func,a,b,xp,yp,np,s,n)
    implicit none
    integer np
    real*8 xp(np),yp(np)
    INTEGER n
    REAL*8 a,b,s,func
    EXTERNAL func
    INTEGER it,j
    REAL*8 del,sum,tnm,x
    if (n.eq.1) then
       s=0.5*(b-a)*(func(a,xp,yp,np)+func(b,xp,yp,np))
    else
       it=2**(n-2)
       tnm=it
       del=(b-a)/tnm
       x=a+0.5*del
       sum=0.
       do 11 j=1,it
          sum=sum+func(x,xp,yp,np)
          x=x+del
11     continue
       s=0.5*(s+(b-a)*sum/tnm)
    endif
    return
  END SUBROUTINE trapzd2

  SUBROUTINE trapzd3(func,a,b,c,s,n)
    implicit none
    INTEGER n
    REAL*8 a,b,c,s,func
    EXTERNAL func
    INTEGER it,j
    REAL*8 del,sum,tnm,x
    if (n.eq.1) then
       s=0.5*(b-a)*(func(a,c)+func(b,c))
    else
       it=2**(n-2)
       tnm=it
       del=(b-a)/tnm
       x=a+0.5*del
       sum=0.
       do 11 j=1,it
          sum=sum+func(x,c)
          x=x+del
11     continue
       s=0.5*(s+(b-a)*sum/tnm)
    endif
    return
  END SUBROUTINE trapzd3
  
  SUBROUTINE midpnt(func,a,b,s,n)
    implicit none
    INTEGER n
    REAL*8 a,b,s,func
    EXTERNAL func
    INTEGER it,j
    REAL*8 ddel,del,sum,tnm,x
    if (n.eq.1) then
       s=(b-a)*func(0.5*(a+b))
    else
       it=3**(n-2)
       tnm=it
       del=(b-a)/(3.*tnm)
       ddel=del+del
       x=a+0.5*del
       sum=0.
       do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11     continue
       s=(s+(b-a)*sum/tnm)/3.
    endif
    return
  END SUBROUTINE midpnt

  SUBROUTINE midpnt2(func,a,b,xp,yp,np,s,n)
    implicit none
    integer np
    real*8 xp(np),yp(np)
    INTEGER n
    REAL*8 a,b,s,func
    EXTERNAL func
    INTEGER it,j
    REAL*8 ddel,del,sum,tnm,x
    if (n.eq.1) then
       s=(b-a)*func(0.5*(a+b),xp,yp,np)
    else
       it=3**(n-2)
       tnm=it
       del=(b-a)/(3.*tnm)
       ddel=del+del
       x=a+0.5*del
       sum=0.
       do 11 j=1,it
          sum=sum+func(x,xp,yp,np)
          x=x+ddel
          sum=sum+func(x,xp,yp,np)
          x=x+del
11     continue
       s=(s+(b-a)*sum/tnm)/3.
    endif
    return
  END SUBROUTINE midpnt2

  SUBROUTINE midinf(funk,aa,bb,s,n)
    INTEGER n
    REAL*8 aa,bb,s,funk
    EXTERNAL funk
    INTEGER it,j
    REAL*8 a,b,ddel,del,sum,tnm,func,x
    func(x)=funk(1./x)/x**2
    b=1./aa
    a=1./bb
    if (n.eq.1) then
       s=(b-a)*func(0.5*(a+b))
    else
       it=3**(n-2)
       tnm=it
       del=(b-a)/(3.*tnm)
       ddel=del+del
       x=a+0.5*del
       sum=0.
       do j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
       end do
       s=(s+(b-a)*sum/tnm)/3.
    endif
    return
  end SUBROUTINE midinf
    
  ! Calculation of Bessel functions...
  
  REAL*8 FUNCTION bessj0(x)
    implicit none
    REAL*8 x
    REAL*8 ax,xx,z
    REAL*8 p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,&
         & s1,s2,s3,s4,s5,s6,y
    SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,&
         & s5,s6
    DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,&
         & -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,&
         & .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
    DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,&
         & 651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2,&
         & s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0,&
         & 59272.64853d0,267.8532712d0,1.d0/
    if(abs(x).lt.8.)then
       y=x**2
       bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*&
            & (s4+y*(s5+y*s6)))))
    else
       ax=abs(x)
       z=8./ax
       y=z**2
       xx=ax-.785398164
       bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*&
            & p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
    endif
    return
  END FUNCTION bessj0

  REAL*8 FUNCTION bessj1(x)
    REAL*8 ax,x,xx,z
    REAL*8 p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,&
         &s1,s2,s3,s4,s5,s6,y
    SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,&
         &s5,s6
    DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,&
         &242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2,&
         &s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0,&
         &99447.43394d0,376.9991397d0,1.d0/
    DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,&
         &.2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,&
         &-.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
    if(abs(x).lt.8.)then
       y=x**2
       bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+&
            &y*(s4+y*(s5+y*s6)))))
    else
       ax=abs(x)
       z=8./ax
       y=z**2
       xx=ax-2.356194491
       bessj1=sqrt(.636619772d0/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*&
            &p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.d0,x)
    endif
    return
  END FUNCTION bessj1

  REAL*8 FUNCTION bessi0(x)
    implicit none
    REAL*8 x
    REAL*8 ax
    REAL*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
    SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
    DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,&
         & 1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
    DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,&
         & 0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,0.2635537d-1,&
         & -0.1647633d-1,0.392377d-2/
    if (abs(x).lt.3.75) then
       y=(x/3.75)**2
       bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
    else
       ax=abs(x)
       y=3.75/ax
       bessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*&
            & (q7+y*(q8+y*q9))))))))
    endif
    return
  END FUNCTION bessi0

  REAL*8 FUNCTION bessk0(x)
    implicit none
    REAL*8 x
    REAL*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
    SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
    DATA p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,&
         & 0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
    DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,&
         & -0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/
    if (x.le.2.0) then
       y=x*x/4.0
       bessk0=(-log(x/2.0)*bessi0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*&
            & (p6+y*p7))))))
    else
       y=(2.0/x)
       bessk0=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*&
            & q7))))))
    endif
    return
  END FUNCTION bessk0

  REAL*8 FUNCTION bessi1(x)
    implicit none
    REAL*8 x
    REAL*8 ax
    REAL*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
    SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
    DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,&
         & 0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
    DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,&
         & -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,-0.2895312d-1,&
         & 0.1787654d-1,-0.420059d-2/
    if (abs(x).lt.3.75) then
       y=(x/3.75)**2
       bessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
    else
       ax=abs(x)
       y=3.75/ax
       bessi1=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*&
            &(q7+y*(q8+y*q9))))))))
       if(x.lt.0.)bessi1=-bessi1
    endif
    return
  END FUNCTION bessi1

   REAL*8 FUNCTION bessk1(x)
     implicit none
     REAL*8 x
     REAL*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
     SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
     DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0,&
          & -0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
     DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,-0.3655620d-1,&
          & 0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/
     if (x.le.2.0) then
        y=x*x/4.0
        bessk1=(log(x/2.0)*bessi1(x))+(1.0/x)*(p1+y*(p2+y*(p3+y*(p4+y*&
             & (p5+y*(p6+y*p7))))))
     else
        y=2.0/x
        bessk1=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*&
             & q7))))))
     endif
     return
   END FUNCTION bessk1

   ! Generate random numbers...

   REAL*8 FUNCTION ran3(idum)
     implicit none
     INTEGER idum
     INTEGER MBIG,MSEED,MZ
     REAL*8 FAC
     PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
     INTEGER i,iff,ii,inext,inextp,k
     INTEGER mj,mk,ma(55)
     SAVE iff,inext,inextp,ma
     DATA iff /0/
     if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
           ii=mod(21*i,55)
           ma(ii)=mk
           mk=mj-mk
           if(mk.lt.MZ)mk=mk+MBIG
           mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
   END FUNCTION ran3

   ! Solve ODE
   
   SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,&
        &     rkqs,xx,yy,n,phiz,fintrp)
!     use structure
     implicit none
     
     INTEGER nbad,nok,nvar,MAXSTP,KMAXX,NMAX
     PARAMETER (MAXSTP=1000,NMAX=2,KMAXX=2000)
     
     REAL*8 eps,h1,hmin,x1,x2,ystart(nvar),TINY
     parameter(TINY=1.d-30)
     EXTERNAL derivs,rkqs
     
     INTEGER i,kmax,kount,nstp
     REAL*8 dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),&
          &     yp(NMAX,KMAXX),yscal(NMAX)
     COMMON /path/ kmax,kount,dxsav,xp,yp
     
     integer n
     real*8 xx(*),yy(*)
     
     real*8 phiz,fintrp
     external phiz,fintrp
     real*8 sum
     
     kmax=KMAXX
     
     x=x1
     h=sign(h1,x2-x1)
     nok=0
     nbad=0
     kount=0
     
     do i=1,nvar
        y(i)=ystart(i)
     end do
     
     if (kmax.gt.0) xsav=x-2.*dxsav
     
     do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        
        do i=1,nvar
           yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
        end do
        
        if(kmax.gt.0)then
           
           if(abs(x-xsav).gt.abs(dxsav)) then
              
              if(kount.lt.kmax-1)then
                 kount=kount+1
                 xp(kount)=x
                 do i=1,nvar
                    yp(i,kount)=y(i)
                 end do
                 xsav=x
              endif            ! if(kount.lt.kmax-1) then
              
           endif               ! if(abs(x-xsav).gt.abs(dxsav)) then
        endif                  ! if(kmax.gt.0) then
        
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
        
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        
        if(hdid.eq.h)then
           nok=nok+1
        else
           nbad=nbad+1
        endif
        
        if((x-x2)*(x2-x1).ge.0.)then
           do 14 i=1,nvar
              ystart(i)=y(i)
14         continue
           if(kmax.ne.0)then
              kount=kount+1
              xp(kount)=x
              do 15 i=1,nvar
                 yp(i,kount)=y(i)
15            continue
           endif

           n=kount
           do i=1,n
              xx(i)=xp(i)
              yy(i)=yp(1,i)
           end do
           
           return
        endif
        if(abs(hnext).lt.hmin) pause 'stepsize smaller than minimum in odeint'
        h=hnext
16      continue
        pause 'too many steps in odeint'
        return
      END SUBROUTINE odeint
      
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
!        implicit none
        INTEGER n,NMAX
        REAL*8 eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
        PARAMETER (NMAX=50)
        INTEGER i
        external derivs
        REAL*8 errmax,h,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,PSHRNK,&
             & ERRCON
        PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)
        h=htry
1       call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
        errmax=0.
        do 11 i=1,n
           errmax=max(errmax,abs(yerr(i)/yscal(i)))
11      continue

        errmax=errmax/eps
        if(errmax.gt.1.)then
           h=SAFETY*h*(errmax**PSHRNK)
           if(h.lt.0.1*h)then
              h=.1*h
           endif
           xnew=x+h
           if(xnew.eq.x)pause 'stepsize underflow in rkqs'
           goto 1
        else
           if(errmax.gt.ERRCON)then
              hnext=SAFETY*h*(errmax**PGROW)
           else
              hnext=5.*h
           endif
           hdid=h
           x=x+h
           do 12 i=1,n
              y(i)=ytemp(i)
12         continue
        return
     endif
  END SUBROUTINE rkqs

  SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
!    implicit none
    INTEGER n,NMAX
    REAL*8 h,x,dydx(n),y(n),yerr(n),yout(n)
    PARAMETER (NMAX=50)
    external derivs
    INTEGER i
    REAL*8 ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),&
         & ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,&
         & B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
    PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,&
         & B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,&
         & B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,&
         & B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,&
         & C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,&
         & DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,&
         & DC6=C6-.25)
    do 11 i=1,n
       ytemp(i)=y(i)+B21*h*dydx(i)
11  continue
    call derivs(x+A2*h,ytemp,ak2)
    do 12 i=1,n
       ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12  continue
    call derivs(x+A3*h,ytemp,ak3)
    do 13 i=1,n
       ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13  continue
    call derivs(x+A4*h,ytemp,ak4)
    do 14 i=1,n
       ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14  continue
    call derivs(x+A5*h,ytemp,ak5)
    do 15 i=1,n
       ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+&
            & B65*ak5(i))
15  continue
    call derivs(x+A6*h,ytemp,ak6)
    do 16 i=1,n
       yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16  continue
    do 17 i=1,n
       yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
17  continue
    return
  END SUBROUTINE rkck
      
  ! Sorting routine
  
  SUBROUTINE indexx(n,arr,indx)
    implicit none
    INTEGER n,indx(n),M,NSTACK
    REAL arr(n)
    PARAMETER (M=7,NSTACK=50)
    INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
    REAL a
    do 11 j=1,n
       indx(j)=j
11  continue
    jstack=0
    l=1
    ir=n
1   if(ir-l.lt.M)then
       do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
             if(arr(indx(i)).le.a)goto 2
             indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13     continue
       if(jstack.eq.0)return
       ir=istack(jstack)
       l=istack(jstack-1)
       jstack=jstack-2
    else
       k=(l+ir)/2
       itemp=indx(k)
       indx(k)=indx(l+1)
       indx(l+1)=itemp
       if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
       endif
       if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
       endif
       if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
       endif
       i=l+1
       j=ir
       indxt=indx(l)
       a=arr(indxt)
3      continue
       i=i+1
       if(arr(indx(i)).lt.a)goto 3
4      continue
       j=j-1
       if(arr(indx(j)).gt.a)goto 4
       if(j.lt.i)goto 5
       itemp=indx(i)
       indx(i)=indx(j)
       indx(j)=itemp
       goto 3
5      indx(l)=indx(j)
       indx(j)=indxt
       jstack=jstack+2
       if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
       if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
       else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
       endif
    endif
    goto 1
  END SUBROUTINE indexx

  REAL*8 FUNCTION expint(n,x)  
    implicit none
    
    INTEGER n,MAXIT  
    REAL*8 x,EPS,FPMIN,EULER  
    PARAMETER (MAXIT=100,EPS=1.e-7,FPMIN=1.e-30,EULER=.5772156649)  
    INTEGER i,ii,nm1  
    REAL a,b,c,d,del,fact,h,psi  
    nm1=n-1  
    if(n.lt.0.or.x.lt.0..or.(x.eq.0..and.(n.eq.0.or.n.eq.1)))then  
       pause 'bad arguments in expint'  
    else if(n.eq.0)then  
       expint=exp(-x)/x  
    else if(x.eq.0.)then  
       expint=1./nm1  
    else if(x.gt.1.)then  
       b=x+n  
       c=1./FPMIN  
       d=1./b  
       h=d  
       do 11 i=1,MAXIT  
          a=-i*(nm1+i)  
          b=b+2.  
          d=1./(a*d+b)  
          c=b+a/c  
          del=c*d  
          h=h*del  
          if(abs(del-1.).lt.EPS)then  
             expint=h*exp(-x)  
             return  
          endif
11     continue  
          pause 'continued fraction failed in expint'  
       else  
          if(nm1.ne.0)then  
             expint=1./nm1  
          else  
             expint=-log(x)-EULER  
          endif
          fact=1.  
          do 13 i=1,MAXIT  
             fact=-fact*x/i  
             if(i.ne.nm1)then  
                del=-fact/(i-nm1)  
             else  
                psi=-EULER  
                do 12 ii=1,nm1  
                   psi=psi+1./ii  
12              continue  
                   del=fact*(-log(x)+psi)  
                endif  
                expint=expint+del  
                if(abs(del).lt.abs(expint)*EPS) return  
13         continue  
	        pause 'series failed in expint'  
      endif
      return  
   END FUNCTION expint
  
end module nrutils_modules
