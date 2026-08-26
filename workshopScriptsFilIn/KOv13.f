c      KO see Springer book by Wilkins, M.
c              "Computer simulations of Dynamic Phenomena"
       parameter(jj = 3600)       !total number of nodes allowed
                                 !this number should be greater than the number
                                 !of nodes that you want so that nodes can be added
                                 !if fracture occurs
       parameter(jmat=50)        !number of materials
       parameter(nmax=4)        !number of placeholders for time: t(0) and t(1) current, t(2) and t(3) new
       integer n, j,imat,ibc(0:jj), icontact,idebug,jdebug,iskip,ipoint
       integer j1(jj),ieostemp(jmat),istrtemp(jmat)
       integer        ieos(2,0:jj),  istr(0:jj),jjj,jzz,ndebug,i,amr
       integer icompact(0:nmax,jj)

       real*8  m(0:jj), r(0:nmax,0:jj), U(0:nmax,0:jj)
       real*8  t(0:nmax),pfrac(0:jj),pvoid
       real*8  phi(0:nmax,0:jj), sigmar(0:nmax,0:jj)
       real*8  sigmao(0:nmax,0:jj),Temp(0:nmax,0:jj)
       real*8  beta(0:nmax,0:jj), P(0:nmax,0:jj), q(0:nmax,0:jj)
       real*8  s1(0:nmax,0:jj),s2(0:nmax,0:jj), s3(0:nmax,0:jj)
       real*8  rho(0:jj), V(0:nmax,0:jj), entropy(0:nmax,0:jj)
       real*8  epsi1(0:nmax,0:jj),epsi2(0:nmax,0:jj),epdt1(0:nmax,0:jj)
       real*8  epsip(0:nmax,3,0:jj)
       real*8  epstrain(0:nmax,0:jj),eqstress(0:nmax,0:jj),eqs
       real*8  eestrain(0:nmax,0:jj),epsilmax(0:jj)
       real*8  e00,e0,e01,e02,e03,e04
       real*8  E(0:nmax,0:jj), K(0:nmax,0:jj)
       real*8  Y(0:jj),Y0(0:jj),mu(0:jj),mu0(0:jj),Ytry
       real*8  deltaZ, Vdot ,dt_min, dr_min,delt_temp
       real*8  deltar,  qbar, ipstrain
       real*8  deltexact
       real*8  a, b, Co, CL, d ,a_min,b_min,rho_min
       real*8  EOS(0:jmat,13),init(jmat,16),bc(-9:9,5),LL(jj),xstart(jj)
       real*8  str(0:jmat,10)
       real*8 SGLbeta(0:jj),SGLn(0:jj),SGLb(0:jj)
       real*8  JOAJO(0:jj),JOBJO(0:jj),JOCJO(0:jj)
       real*8  JOMJO(0:jj),JONJO(0:jj),JOTJO(0:jj),JOPOS(0:jj)
       real*8  xx,xa,xb,deltat_0,deltat,delt,tstop
c       real*8  U_0,P_0,V_0,E_0,rho_0,a_0,up
c       real*8  U_1,P_1,V_1,E_1,a,
       real *8 zero
       real*8  aa,ww,rr,bb,dti0,dti1,Usave
       real*8  phiv,phij,betav,betaj,dthalf
       real*8  qtotal,mvtotal,ketotal,ietotal,etotal,ke3total
       real*8  bs1,bs2,bs3,bs4,bs5,bs6,bs7,bs8,bs9,bs10,bs11,bs12
       real*8  tskip,dtskip
       real*8  rho_local,stemp,ctemp,gtemp,v0,v00,vv
c       real*8  alpha
c      real*8  pe,ps,ap,ae
       real*8  k1,k2,k3,gamma0
       real*8  Us,up,PH,EH,TH,strain,P0,T0
       real*8  En2j1,diffE
       character*20 name
       character*256 text
c
c     Initalize some key varables
       Co        = 2.0d0  !artifical viscoity constants, default=2
       CL        = 1.0d0  !artifical viscoity constants, default=1
       pvoid     = 0.d0   !pressure in the void
       iskip     = 1000    !number of iterations to skip to echo to screen, must be >= 1
       icontact  = 0      !contact flag this should not be changed but i think it could
                          ! be elimated if the code were cleaned up.
       ipfrac    = 0    ! this will be set to 1 if any pfrac is non-zero
       idebug    = 0      ! turn on debug
       ndebug    = 1      ! what time steup to you want to print to screen
       jdebug    = 2      ! what node do you want to look at when debugging
       d         = 1.d0   !defines geometry (1-1D planar, 2-1D cylinderical-never tested), see willkins
       deltat_0  = 0.1d-7
c
       tskip     = 0.d0   ! this is the value at which data starts getting written to file
       dtskip    = 0.01  ! This is the amount of time skipped between writes to file
       amr       = 1      ! this is an adaptive mesh refinement (1 is on - 0 is off)
c---------------------------------------------------------------------
c     zero all variables
c
c      print *,'Initializing variables ...'
       zero = 0.d0
       do n = 0,nmax
          do j = 0,jj
             ibc(j)      = 9  ! these are non-used nodes
             U(n,j)      = zero
             r(n,j)      = zero
             t(n)        = deltat_0*dfloat(n)
             deltat      = deltat_0
             phi(n,j)    = zero
             sigmar(n,j) = zero
             sigmao(n,j) = zero
             beta(n,j)   = zero
             q(n,j)      = zero
             s1(n,j)     = zero
             s2(n,j)     = zero
             s3(n,j)     = zero
             epsi1(n,j)  = zero
             epsi2(n,j)  = zero
             epsip(n,1,j)= zero
             epsip(n,2,j)= zero
             epsip(n,3,j)= zero
             epstrain(n,j)=zero
             K(n,j)      = zero
             Y(j)        = zero
             deltaZ      = zero
             pfrac(j)    = zero
             Temp(n,j)   = zero
             entropy(n,j)= zero
            enddo
       enddo
c
      do j=1,jmat
       do i=1,7
        eos(j,i)= 0.d0
       enddo
       do i=1,4
        init(j,i)= 0.d0
       enddo
      enddo
      do j=1,jj
       do jzz=1,2
        ieos(jzz,j)= 0
        istr(    j)= 0 
       enddo
      enddo
c---------------------------------------------------------------------
c     Read input Data file
c
c      print *,'running'
c      read(*,*) foo
      imat = 0
      jsum = 0
      iReadInput = 1
      name='ko.dat'
      open (unit=33, file=name, form='formatted', status='unknown')
      open (unit=31,file='ko.in',form='formatted',status='unknown')
      open (unit=35, file='ke.dat', form='formatted', status='unknown')
      print *,'opened file'
c      read(*,*) foo
      Read(31,'(A240)') text
      print *,text
      Read(31,'(A240)') text
      print *,text
c      read(*,*) foo
      do while(iReadInput .eq. 1) ! Loop through ko.in and read material input data 
 5    imat=imat+1
      print *,'reading',imat
c      read(*,*) foo
      Read(31,'(A240)') text
      print *,Text
c      print *,'hello world'
c             print *,'testing text'
c             read(*,*) foo
      if(text(1:7) .eq. '       ') then
       iReadInput = 0
      else
       print *,'there is a material insert to be read.'
       ieos(1,jj-(imat-1)) = imat

*       Read(text,'(I4,I5,I6,2e9.5,5e7.3,3e7.4,10e7.4,3e5.2,7e10.4)') !g77
       Read(text,*)
     & ieostemp(imat),istrtemp(imat),
     & j1(imat),LL(imat),xstart(imat),
     & init(imat,1),init(imat,2),init(imat,3),
     & init(imat,4),init(imat,5),
     & eos(imat,1),eos(imat,2),eos(imat,3),eos(imat,4),
     & eos(imat,5),eos(imat,6),eos(imat,7),eos(imat,8),eos(imat,9),
     & eos(imat,10),eos(imat,11),eos(imat,12),eos(imat,13),
     &  str(imat,1),str(imat,2), str(imat,3),str(imat,4),
     &  str(imat,5),str(imat,6), str(imat,7),str(imat,8),
     &  str(imat,9),str(imat,10)
       print *,'Input line ',imat,' as follows:'
       print *,
     & ieostemp(imat),istrtemp(imat),
     & j1(imat),LL(imat),xstart(imat),
     & init(imat,1),init(imat,2),init(imat,3),
     & init(imat,4),init(imat,5),
     & eos(imat,1),eos(imat,2),eos(imat,3),eos(imat,4),
     & eos(imat,5),eos(imat,6),eos(imat,7),eos(imat,8),eos(imat,9),
     & eos(imat,10),eos(imat,11),eos(imat,12),eos(imat,13),
     &  str(imat,1),str(imat,2), str(imat,3),str(imat,4),
     &  str(imat,5),str(imat,6), str(imat,7),str(imat,8),
     &  str(imat,9),str(imat,10)
       write(*,'(I4,I5,I6,7e9.3,3e9.4,10e9.4,3e7.2,7e10.4)') !g77
     & ieostemp(imat),istrtemp(imat),
     & j1(imat),LL(imat),xstart(imat),
     & init(imat,1),init(imat,2),init(imat,3),
     & init(imat,4),init(imat,5),
     & eos(imat,1),eos(imat,2),eos(imat,3),eos(imat,4),
     & eos(imat,5),eos(imat,6),eos(imat,7),eos(imat,8),eos(imat,9),
     & eos(imat,10),eos(imat,11),eos(imat,12),eos(imat,13),
     &  str(imat,1),str(imat,2), str(imat,3),str(imat,4),
     &  str(imat,5),str(imat,6), str(imat,7),str(imat,8),
     &  str(imat,9),str(imat,10)
       if ( int(j1(imat)/2) .ne. float(j1(imat))/2.) then
        print *,'Error: nodes specificed for a material must be even'
        read(*,*) foo
       endif
       jsum = jsum + j1(imat)
       print *,sum,j1(imat)
       if (jsum .ge. jj) then
        print *,'Error: Either increase jj in the code and recompile 
     &               OR reduce the # nodes in kov3.in'
        read(*,*) foo
       endif
c      read(*,*) foo
       endif
       enddo
       Print *,'End material read from ko.in: ',imat-1,' Materials read'
       imat = imat-1
c      read(*,*) foo
c
c---------------------------------------------------------------------
c     Input boundary conditions
c
      Read(31,'(A240)') text
      print *,text
c      read(*,*) foo
      Read(31,'(A240)') text
      print *,text
c      read(*,*) foo
      Read(31,'(A240)') text
      print *,text
      Read(text,'(I7,5e7.3)') ! g77
c      Read(text,*) ! g90
     & ibc(0),bc(-1,1),bc(-1,2),bc(-1,3),bc(-1,4),bc(-1,5)
      print *,'read Boundary Conditions...',ibc(0)
c      read(*,*) foo
      Read(31,'(A150)') text
      print *,text
      Read(text,'(I7,5e7.3)')  ! g77
c      Read(text,*)           ! f90
     & ibc(jj),bc( 1,1),bc( 1,2),bc( 1,3),bc( 1,4),bc( 1,5) !note that ibc is stored at jjj just temporarily
      print *,'read Boundary Conditions...',ibc(jj)
c
      Read(31,'(A113)') text
      print *,text
      Read(31,'(A113)') text
      print *,text
      Read(31,'(A113)') text
      print *,text
      Read(text,*) tstop,dtskip ! this is the time to run to and the time skip between writes to ko.dat
c      Read(text,*) tstop  ! f90
      close(31)
c---------------------------------------------------------------------
c      Distritazation LOOP
c
      print *,'distrize ..'
c      read(*,*) foo
        r(0,1) = xstart(1)
        r(1,1) = xstart(1)
      do jjj=1,imat
c       print *,'discretize ..',jjj,LL(jjj),xstart(1),ipoint
       deltar = LL(jjj)/dfloat(j1(jjj)/2)
       if (jjj .eq. 1) then  ! this if assigns the ipoint to the first node of the material
                             ! and checks to see it materials are initally in contact.
                             ! for 1D mesoscale simulations there could be a gap between materials
        ipoint = 0           ! material pointer
       elseif (abs(r(0,ipoint+j1(jjj-1))-xstart(jjj)) .lt. 1.e-5) then
c       no gap between materials
        r(0,ipoint+j1(jjj-1)) = xstart(jjj)
        r(1,ipoint+j1(jjj-1)) = xstart(jjj)
        ibc(ipoint+j1(jjj-1)) = 0 !this is a regular node with j-1 and j+1 neighbor
        ipoint                = ipoint+j1(jjj-1)
       elseif ( r(0,ipoint+j1(jjj-1)) .lt. xstart(jjj)) then
c       inital gap between materials
        ipoint      = ipoint+j1(jjj-1)+2 ! notice there are two nodes gap
           ! nodes so that the even/odd or whole/half step counting
           ! remains consistent 
        ibc(ipoint  ) = -2  ! here the gap between materials is set to BC = 2 and -2
        ieos(2,ipoint-1)=0  ! eos type: 1-MG, etc.
        istr(  ipoint-1)=0
        ibc(ipoint-2) =  2
       else
       print *,'Input Geometry Error!'
       print *,'Materials initially overlapping'
       read(*,*) foo
       endif
c
       r(0,ipoint)  = xstart(jjj)
       r(1,ipoint)  = xstart(jjj)
       ibc(1)       = 0
       U(0,ipoint)  = init(jjj,2)
       U(1,ipoint)  = init(jjj,2)
c
       do j=ipoint+2,ipoint+j1(jjj),2      !node definitions
          ibc(j)       = 0      !this defines the node as a centeral difference
          r(0,j)       = r(0,j-2) + deltar
          r(1,j)       = r(1,j-2) + deltar
          U(0,j)       = init(jjj,2)
          U(1,j)       = init(jjj,2)
       enddo
c
       do j=ipoint+1,ipoint+j1(jjj)-1,2    !cell definitions initial conditions
          ibc(j)       = 0      !this defines the cell as a centeral difference
          r(0,j)       = 0.5d0*(r(0,j+1)+r(0,j-1))
          r(1,j)       = 0.5d0*(r(1,j+1)+r(1,j-1))
          P(0,j)       = init(jjj,1)
          P(1,j)       = init(jjj,1)
          rho(j)       = init(jjj,3)  !rho is not updated in time therefore it is always rho0
          E(0,j)       = init(jjj,4)  ! this gets overwritten below
          E(1,j)       = init(jjj,4)
          ieos(1,j)    = jjj            ! material number
          ieos(2,j)    = ieostemp(jjj)  ! type of equation of statep
          istr(  j)    = istrtemp(jjj)  !info was stored there (RHS) just temp.
          if (istr(j) .eq. 3) then ! true for Johnson-Cook
               Y(j)    = str(jjj,4) ! AJO from ko.in
          else
               Y(j)    = eos(jjj,5) ! yield from ko.in
          endif
c          print *,'Yield',j,jjj,istr(j),y(j),str(jjj,4),eos(jjj,5)
          mu(j)        = eos(jjj,6)
          pfrac(j-1)   = eos(jjj,7)
          pfrac(j)     = eos(jjj,7)
          pfrac(j+1)   = eos(jjj,7)
          if (eos(jjj,7) .ne.  0.0d0) ipfrac =1 
          SGLbeta(j)   = str(jjj,1) ! Stienberg Gun plastic strain coefficient
          SGLn(j)      = str(jjj,2) ! Stienberg Gun plastic strain exponent
          SGLb(j)      = str(jjj,3) ! Stienberg Gun plastic strain-rate coefficient
          SGLbeta(j+1) = str(jjj,1) ! Stienberg Gun plastic strain coefficient
          SGLn(j+1)    = str(jjj,2) ! Stienberg Gun plastic strain exponent
          SGLb(j+1)    = str(jjj,3) ! Stienberg Gun plastic strain-rate coefficient
          JOAJO(j)     = str(jjj,4) ! Johnson-Cook plastic strain coefficient
          JOBJO(j)     = str(jjj,5) ! Johnson-Cook plastic strain exponent
          JOCJO(j)     = str(jjj,6) ! Johnson-Cook strain-rate coefficient
          JOMJO(j)     = str(jjj,7) ! Johnson-Cook plastic strain coefficient
          JONJO(j)     = str(jjj,8) ! Johnson-Cook plastic strain exponent
          JOTJO(j)     = str(jjj,9) !  Johnson-Cook coefficient
          JOPos(j)     = str(jjj,10) !  Johnson-Cook plastic strain coefficient
          JOAJO(j+1)   = str(jjj,4) !  Johnson-Cook plastic strain coefficient
          JOBJO(j+1)   = str(jjj,5) ! Johnson-Cook plastic strain exponent
          JOCJO(j+1)   = str(jjj,6) ! Johnson-Cook plastic strain-rate coefficient
          JOMJO(j+1)   = str(jjj,7) ! Johnson-Cook plastic strain coefficient
          JONJO(j+1)   = str(jjj,8) ! Johnson-Cook plastic strain coefficient
          JOTJO(j+1)   = str(jjj,9) ! Johnson-Cook plastic strain exponent
          JOPos(j+1)   = str(jjj,10) ! Johnson-Cook plastic strain-rate coefficient
c       print *,JOAJO(j+1),JOBJO(j+1),JOCJO(j+1)
          enddo
          ieos(1,jj-(jjj-1))     = 0   !this assigns values to the last cell center
          ieos(2,jj-(jjj-1))     = 0
          ibc(j1(jjj)+ipoint)    = 2
          ibc(j1(jjj)+1+ipoint)  = 9
      enddo  ! jjj loop through imat
          ibc(j1(imat)+ipoint)    = ibc(jj)  !this assigns values to the last node
c          print *,'ipoint',j1(imat)+ipoint,ibc(jj)
          ibc(jj)=9
c
c      print *,'Assign mass to nodes'
c
        do j=0,jj-2,2
             if (ibc(j+1) .eq. 0) then
              m(j+1) = rho(j+1)*((r(0,j+2)**d- r(0,j)**d)/d)
             endif
        enddo
c        read(*,*) foo
c---------------------------------------------------------------------
c set initial specific volume
      do j=0,jj-2,2
       if (ibc(j+1) .eq. 0) then
        V(0,j+1)=rho(j+1)*((r(0,j+2)**d-r(0,j)**d)/d)/m(j+1)
        V(1,j+1)=rho(j+1)*((r(0,j+2)**d-r(0,j)**d)/d)/m(j+1)
c
c        E(0,j+1)=P(0,j+1)/(v(0,j+1)*eos(ieos(1,j+1),4)-1)
c        E(1,j+1)=P(1,j+1)/(v(1,j+1)*eos(ieos(1,j+1),4)-1)
        Temp(0,j+1) = 300.d0 ! this is set on pg 61 of wilkins start of last paragraph
        Temp(1,j+1) = 300.d0
        E(0,j+1)    = 0.d0 !                  -eos(ieos(1,j+1),8)*Temp(0,j+1)
        E(1,j+1)    = 0.d0 !                  -eos(ieos(1,j+1),8)*Temp(1,j+1)
        P(0,j+1)    =-rho(j+1)*eos(ieos(1,j+1),4)
     &                        *eos(ieos(1,j+1),8)*Temp(0,j+1)
        P(1,j+1)    =-rho(j+1)*eos(ieos(1,j+1),4)
     &                        *eos(ieos(1,j+1),8)*Temp(0,j+1)
        P(0,j+1)    = 0.d0
        P(1,j+1)    = 0.d0
c
        entropy(0,j+1) = 6.8d-5 ! i got this from a VT website EES thing
        entropy(1,j+1) = 6.8d-5 ! units mbar-cc/K/g
c        print *,'Initial Conditions'
c        print *,'Pres [Mbar]',P(0,j+1)
c        print *,'Engr [MBar-cc/cc] and Cv',E(0,j+1),eos(ieos(i,j+1),8)
c        print *,'Volm [cc]',((r(0,j+2)**d-r(0,j)**d)/d)
c        print *,'temp*[K]',Temp(0,j+1),rho(j+1),m(j+1),v(0,j+1)
c        print *,'entropy',entropy(0,j+1),eos(ieos(1,j+1),8)
c        read(*,*) foo
c        print *,'here1',j
c        icompact(0,j-2)  = 0 ! compaction flag for use with snow plow type models initilized to zero
c      print *,'here2',j
c        icompact(1,j-2)  = 0 ! compaction flag for use with snow plow type models initilized to zero
c      print *,'here3',j
c        icompact(2,j-2)  = 0 ! compaction flag for use with snow plow type models initilized to zero
c        print *,'i am here'
c	read(*,*) foo
       endif
      enddo
c      read(*,*)
c---------------------------------------------------------------------
c     Output Initial Conditions to screen
c
      n=1
c      write(*,98)
c     &'n','rad','rad','Vel','Vel'
c     &,'pres','Vol'
c      do j=0,jj-2
c       write(*,'(I4,0I2,10f15.7)')
c       write(*,*)
c     &  j,  !ieos(1,j+1),ieos(2,j+1),
c     &  ibc(j),ibc(j+1),ibc(j+2),
c     &  r(n,j),r(n,j+1),r(n,j+2),U(n,j),P(n,j+1),
c     &  v(n,j+1),m(j+1),rho(j+1),Y(j+1)
c       enddo
c
       print *,'Start Main Loop,  goto tstop = ',tstop
c      read(*,*) foo
c
 98   format(2A5,11x,A9,A9,A9,A10,a8)
 99   format(1I5,11x,9e9.2,A10,A8)
c
c*****************************************************************************
c******************************** MAIN LOOP **********************************
c*****************************************************************************
c
       ncount = -1
c       if (t(n) .le. tstop) then
 1      ncount = ncount + 2 ! when this use to count to an integer, ncount was the max integer
        n     =  1    ! this use to be the time step counter when all n data was stored
                      ! n must remain n=1 because the time data is no longer stored

      if( m(1) .ne. m(3) ) then
c this can happen by roundoff error. adjust the number of nodes slightly
       print *,'Masses do not match',n,m(1),m(3)
       print *,rho(j+1),r(0,j+2),r(0,j),d,r(0,j+2)-r(0,j)
       read(*,*) foo
      endif
c
cccccccccccccccc   Contact Check   ccccccccccccccccccccccc
c
      if (1 .eq. 1) then
       deltmin = 1.d0 ! find smallest time to contact, initially set this to something large
 3     do j=0,jj,2
      if (ibc(j) .ne. 9) then
        r(n+2,j) =   r(n,j) + U(n+1,j)*deltat    !deltat is full time step between t(n) and t(n+2)
        r(n+1,j) = ( r(n,j) + r(n+2,j) ) / 2.d0
      endif
      enddo
c  This loop is looking for the contact with the smallest time step
      do j=2,jj-2,2
      if ( ibc(j) .ne. 9 .and. r(n+2,j) .lt. r(n+2,j-2) .and.
     &     ibc(j-1) .eq. 9) then ! ibc(j)=9 is for void
          icontact = 2 !this is a flag to let the contact connection
                       !know that the nodes should be shifted and joined 
          deltExact = (r(n,j)-r(n,j-2))/(U(n+1,j-2)-U(n+1,j))/2.d0 !  calculate 
                       !  time step that exact makes node r(n,j-2) land on r(n,j)
          if (deltexact .le. deltmin) then
            deltmin  = deltexact
            jcontact = j ! jcontact is the node to be removed
          endif
c          print *,'deltat,deltexact',deltat,deltmin
          print *,'Contact at ncount=', ncount,' and j=',j
c          print *,'velocity',U(n+1,j),U(n+1,j-2)
c          print *,'r(n+2,j-2),r(n+2,j),r(n+2,j+2)',n,j,
c     &             r(n  ,j-2),r(n  ,j),r(n  ,j+2),
c     &             r(n+1,j-2),r(n+1,j),r(n+1,j+2),
c     &             r(n+2,j-2),r(n+2,j),r(n+2,j+2)

          deltat = deltmin
          goto 3  ! go back, advance the node locations and see if there 
                  ! is still a node with overlap.  There can still be
                  ! overlap if there are other nodes that were initially
                  ! closer together than the two found
      endif
      enddo
c
cccccc     Begin the shifting of nodes due to void remove from contact
c
       if (icontact .eq. 2 ) then
       Usave = (m(jcontact+1)*U(n+1,jcontact)+m(jcontact-3)
     &             *U(n+1,jcontact-2))/(m(jcontact-3)+m(jcontact+1))
c
c          do jjj = 0,jj !notice this is shifting every node.
c           print *,'pre-shift',jjj,ibc(jjj),ieos(1,jjj),ieos(2,jjj),
c     &                          r(0,jjj),r(1,jjj),r(2,jjj),r(3,jjj)
c          enddo
c
c          print *,'Start contact shift',j,jj
          do jjj = jcontact-2,jj-6 !notice this is shifting every node.
              ibc(    jjj  )  =    ibc(    jjj+2)
             ieos(1,  jjj  )  =   ieos(1,  jjj+2)
             ieos(2,  jjj  )  =   ieos(2,  jjj+2)
             SGLbeta(jjj)     = SGLbeta(jjj+2) ! Stienberg Gun plastic strain coefficient
             SGLn(jjj)        = SGLn(jjj+2) ! Stienberg Gun plastic strain exponent
             SGLb(jjj)        = SGLb(jjj+2) ! Stienberg Gun plastic strain-rate coefficient
             JOAJO(jjj)       = JOAJO(jjj+2) !  Johnson-Cook plastic strain coefficient
             JOBJO(jjj)       = JOBJO(jjj+2) ! Johnson-Cook plastic strain exponent
             JOCJO(jjj)       = JOCJO(jjj+2) ! Johnson-Cook plastic strain-rate coefficient
             JOMJO(jjj)       = JOMJO(jjj+2) ! Johnson-Cook plastic strain coefficient
             JONJO(jjj)       = JONJO(jjj+2) ! Johnson-Cook plastic strain coefficient
             JOTJO(jjj)       = JOTJO(jjj+2) ! Johnson-Cook plastic strain exponent
             JOPos(jjj)       = JOPos(jjj+2) ! Johnson-Cook plastic strain-rate coefficient
c
                m( jjj  ) =      m(jjj+2)
              rho( jjj  ) =    rho(jjj+2)
                 Y(jjj  ) =      Y(jjj+2)
             pfrac(jjj  ) =  pfrac(jjj+2)  !node value
             pfrac(jjj  ) =  pfrac(jjj+2)  !Cell value
c
                r(0 ,jjj )  =     r(0 ,jjj+2)
                U(0 ,jjj )  =     U(0 ,jjj+2)
              phi(0 ,jjj ) =    phi(0 ,jjj+2)
             beta(0 ,jjj ) =   beta(0 ,jjj+2)
           sigmar(0 ,jjj ) = sigmar(0 ,jjj+2)
           sigmao(0 ,jjj ) = sigmao(0 ,jjj+2)
                V(0 ,jjj ) =      V(0 ,jjj+2)
               s1(0 ,jjj ) =     s1(0 ,jjj+2)
               s2(0 ,jjj ) =     s2(0 ,jjj+2)
               s3(0 ,jjj ) =     s3(0 ,jjj+2)
                E(0 ,jjj ) =      E(0 ,jjj+2)
c
                r(1 ,jjj ) =      r(1 ,jjj+2)
                U(1 ,jjj ) =      U(1 ,jjj+2)
              phi(1 ,jjj ) =    phi(1 ,jjj+2)
             beta(1 ,jjj ) =   beta(1 ,jjj+2)
           sigmar(1 ,jjj ) = sigmar(1 ,jjj+2)
           sigmao(1 ,jjj ) = sigmao(1 ,jjj+2)
                V(1 ,jjj ) =      V(1 ,jjj+2)
               s1(1 ,jjj ) =     s1(1 ,jjj+2)
               s2(1 ,jjj ) =     s2(1 ,jjj+2)
               s3(1 ,jjj ) =     s3(1 ,jjj+2)
                E(1 ,jjj ) =      E(1 ,jjj+2)
          if (jjj .eq. jj) then
            print *,r(0,jjj),r(0,jjj)
            read (*,*)
          endif
          enddo
                U(n+1 ,jcontact  -2) = Usave !setting jcontact node values
             pfrac(    jcontact  -2) = 1.d-2
                     ibc(jcontact-2) = 0
                     ibc(jcontact-1) = 0
          icontact = 0  ! this flag was on in order to force the shift
c          do jjj = 0,jj !notice this is shifting every node.
c           print *,'post shift',jjj,ibc(jjj),ieos(1,jjj),ieos(2,jjj),
c     &                          r(0,jjj),r(1,jjj),r(2,jjj),r(3,jjj)
c          enddo
          goto 3  ! this loops through again to make sure there is not a
          ! second contact that occures at the same instant in time
       endif
       endif
cccccccccccccccccccccc End of contact ccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccc  Advance Time step  cccccccccccccccccccccccccccccc
c
c     advance time step, for first calculation timestep is set by deltat_0
c     after that the time step is calculated after the end of the current time step, deltat
c
      if (icontact .eq. 2) then      ! if there was contact time step adjusted to perfectly
         t(n+1) = t(n)   + deltat/2.d0    ! have nodes touch at the end of the time step
         t(n+2) = t(n+1) + deltat  !*2.0d0
         if  (dti1 - deltat .lt. 0. .or. t(n+2)-t(n+1) .lt. 0.) then
          print *,'Time step Error from contact'
          print *, dti1,deltat,dti1 - deltat , t(n+2)-t(n+1)
 10       read(*,*) foo
          go to 10
         endif
      else
         t(n+1) = t(n) + deltat/2.0d0
         t(n+2) = t(n) + deltat  !*2.0d0
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     2.    Conservation of Momumtum  - start j
c
c  Boundary Conditions
c   -/+ 1 Free End, stress free
c   -/+ 2 and 3 Part of contact check
c   -/+ 4 Fixed end
c
c System is failing in this do loop line 450	

      do j=0,jj-1,2
      if (ibc(j+1) .eq. 0) then !this is true for central difference cells
       sigmar(n,j+1) = (-(P(n,j+1)+q(n-1,j+1))+s1(n,j+1))
       sigmao(n,j+1) = (-(P(n,j+1)+q(n-1,j+1))+s2(n,j+1))
       endif
      enddo
c	
      do j=0,jj,2
c
       delt = t(n+1)-t(n-1)    !this is different than deltat, deltat is the true timestep
      if (ibc(j) .ne. 9) then

       if (ibc(j) .eq. -1) then
        phi(n,j) = (rho(j+1)*((r(n,j+2)-r(n,j))/V(n,j+1)))/2.d0
        beta(n,j) = (sigmar(n,j+1)-sigmao(n,j+1))*V(n,j+1)
     &             /(r(n,j+1)*rho(j+1))
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))
     &  *( sigmar(n,j+1)+bc(-1,1))+delt*beta(n,j)*(d-1.d0)
c       print *,'BC thing',bc(-1,1)
       elseif (ibc(j) .eq. 1) then
        phi(n,j) =  rho(j-1)*((r(n,j)-r(n,j-2))/V(n,j-1))/2.d0
        beta(n,j) = (sigmar(n,j-1)-sigmao(n,j-1))*V(n,j-1)
     &             /(r(n,j-1)*rho(j-1))
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))
     &  *(-bc( 1,1)-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0)
c       U(n+1,j) = 0.d0
c
c this is the new way which treats the voids like inner or outer bc's
       elseif (ibc(j) .eq. -2 .or. ibc(j) .eq. -3) then
        phi(n,j) = (rho(j+1)*((r(n,j+2)-r(n,j))/V(n,j+1)))/2.d0
        beta(n,j) = (sigmar(n,j+1)-sigmao(n,j+1))*V(n,j+1)
     &             /(r(n,j+1)*rho(j+1))
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))
     &  *(sigmar(n,j+1)+pvoid)+delt*beta(n,j)*(d-1.d0)
       elseif (ibc(j) .eq. 2) then
        phi(n,j) =  rho(j-1)*((r(n,j)-r(n,j-2))/V(n,j-1))/2.d0
        beta(n,j) = (sigmar(n,j-1)-sigmao(n,j-1))*V(n,j-1)
     &             /(r(n,j-1)*rho(j-1))
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))
     &  *(pvoid-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0)
       elseif (ibc(j) .eq. 0) then
        phi(n,j) = (0.5d0)*(rho(j+1)*((r(n,j+2)-r(n,j))/V(n,j+1))+
     &   rho(j-1)*((r(n,j)-r(n,j-2))/V(n,j-1)))
        beta(n,j)=( (sigmar(n,j+1)-sigmao(n,j+1))*(V(n,j+1)/rho(j+1))/
     &  (0.5d0*(r(n,j+2)+r(n,j  ))) +
     &            (sigmar(n,j-1)-sigmao(n,j-1))*(V(n,j-1)/rho(j-1))/
     &  (0.5d0*(r(n,j  )+r(n,j-2))) )/2.d0
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))
     &  *(sigmar(n,j+1)-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0)
       elseif (ibc(j) .eq. -4) then
        U(n+1,j) = 0.d0
c       if ( U(n+1,j) .lt. 1.d-5) U(n+1,j) = 0.d0
       elseif (ibc(j) .eq. 4) then
        U(n+1,j) = 0.d0   ! rigid wall
c       if ( U(n+1,j) .lt. 1.d-5) U(n+1,j) = 0.d0
       else
        print *,'Momentum ERROR!!!'
        print *,j,ibc(j)
        read(*,*) foo
       endif
c
      endif
      enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     3.  Position  - within j
c
      do j=0,jj,2
      if (ibc(j) .ne. 9) then
        r(n+2,j) =  r(n,j)+U(n+1,j)*deltat  !*2.d0
        r(n+1,j) = (r(n,j)+r(n+2,j))/2.d0
      endif
      enddo
c debug block
        if (idebug .eq. 1 .and. ncount .ge. ndebug) then
        j=jdebug
        print *,'j=',j,'n=',n,'ncount=',ncount
        print *,'time',t(n-1),t(n),t(n+1),t(n+2)
        print *,'deltat',deltat,delt
        print *,'vel',U(n,j),U(n+1,j)
        print *,'r',r(n,j-2),r(n,j),r(n,j+2)
        print *,'V',V(n,j+1),V(n,j),V(n,j-1)
        print *,'rho',rho(j-1),rho(j),rho(j+1)
        print *,'sigmar',sigmar(n,j+1),sigmar(n,j-1)
        print *,'sigmao',sigmao(n,j-1),sigmao(n,j+1)
        print *,'4',phi(n,j),beta(n,j),U(n-1,j)
        print *,'5',(deltat/phi(n,j))*(sigmar(n,j+1)-0.d0)
        read(*,*) foo
        endif
c
      if (1 .eq. 1) then
       do j=2,jj-2,2
        if ( ibc(j) .ne. 9 .and. r(n+2,j) .lt. r(n+2,j-2) ) then
          print *,'Contact Check Failed',ncount,j
          print *,'Try lowering time step'
          print *,'r(n+2,j-2),u(n+1,j-2)',r(n+2,j-2),U(n+1,j-2)
          print *,'r(n+2,j  ),u(n+1,j  )',r(n+2,j)  ,U(n+1,j)
          print *,'deltat',deltat
          print *,bc(-1,1),bc(1,1)
          read(*,*) foo
       endif
      enddo
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     3.  Relative VOLUME   - within j
c
c      deltat     =  t(n+2)-t(n)
c      print *,'line 626'
c	read(*,*) foo
      do j=0,jj-2,2
      if (ibc(j+1) .eq. 0) then
c
        r(n+1,j+1) = ( r(n+1,j)+r(n+1,j+2) ) /2.d0
        r(n+2,j+1) = ( r(n+2,j)+r(n+2,j+2) ) /2.d0
c
        V(n+2,j+1)=rho(j+1)*((r(n+2,j+2)**d-r(n+2,j)**d)/d)/m(j+1)
        V(n+1,j+1)=rho(j+1)*((r(n+1,j+2)**d-r(n+1,j)**d)/d)/m(j+1)
c
        if ( V(n+2,j+1) .eq. 0.) then
        print *,'------------- Zero volume error! #1-----------------'
        print *,'volume',n,j+1,r(n+2,j+2),r(n+2,j),V(n+2,j+1)
        read(*,*) foo
        endif
        if ( V(n+1,j+1) .eq. 0.) then
        print *,'------------- Zero volume error! #2-----------------'
        print *,'volume',n,j+1,r(n+1,j+2),r(n+1,j),V(n+1,j+1)
        read(*,*) foo
        endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     4. VELOCITY STRAINS  - finish j
c
         epsi1(n+1,j+1) = (U(n+1,j+2)-U(n+1,j))/(r(n+1,j+2)-r(n+1,j))
         epdt1(n+1,j+1) = (epsi1(n+1,j+1)-epsi1(n-1,j+1))/deltat !this is 1st order accurate
        if (d .eq. 1.) then
         epsi2(n+1,j+1) = 0.d0  ! this enforces plane strain
        else
         epsi2(n+1,j+1) = (U(n+1,j+2)+U(n+1,j))/(r(n+1,j+2)+r(n+1,j))
        endif
      endif

      enddo  ! end j loop for volume
c--------------------  AMR Feature  -------------------------------

        if (idebug .eq. 1 .and. n .ge. ndebug) then
        j=jdebug
        print *,'r',j,r(n+2,j+1),m(j+1)
        print *,'epsil',j,V(n+1,j+1),V(n+2,j+1),
     &                    epsi1(n+1,j+1),epsi2(n+1,j+1)
        read(*,*) foo
        endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     5a. STRESSES  - start j
c
c      print *,'line 679'
      do j=0,jj-2,2
      if (ibc(j+1) .eq. 0) then
c       deltat = t(n+2)-t(n)
       if (ieos(2,j+1) .ne. 3) then  
        s1(n+2,j+1) = s1(n,j+1) + 2.d0*mu(j+1)
     &             * ( epsi1(n+1,j+1)*deltat
     &                 - (V(n+2,j+1)-V(n,j+1))/(3.d0*V(n+1,j+1)) )
        s2(n+2,j+1) = s2(n,j+1) + 2.d0*mu(j+1)
     &             * ( epsi2(n+1,j+1)*deltat
     &                 - (V(n+2,j+1)-V(n,j+1))/(3.d0*V(n+1,j+1)) )
        s3(n+2,j+1) = -(s1(n+2,j+1)+s2(n+2,j+1))
c
        s1(n+1,j+1) = ( s1(n+2,j+1) + s1(n,j+1) )/2.d0
        s2(n+1,j+1) = ( s2(n+2,j+1) + s2(n,j+1) )/2.d0
        s3(n+1,j+1) = ( s3(n+2,j+1) + s3(n,j+1) )/2.d0
       else     !Gamma law ideal gas Newtonian Stress
        s1(n+2,j+1) = 4.d0*1.8d-10*epsi1(n+1,j+1)/3.d0
     &             -1.387e-6* (V(n+2,j+1)-V(n,j+1))/V(n+1,j+1) 
c        s1(n+2,j+1) = s1(n,j+1)+4.d0*1.8d-10*epdt1(n+1,j+1)*deltat/3.d0
c     &             -1.387e-6* (V(n+2,j+1)-V(n,j+1))/(3.d0*V(n+1,j+1)) 
c
c        print *,epdt1(n+1,j+1),epsi1(n+1,j+1),epsi1(n-1,j+1)
c
        s1(n+1,j+1) = ( s1(n+2,j+1) + s1(n,j+1) )/2.d0
      endif ! ieos check
c
      endif
      enddo
        if (idebug .eq. 1 .and. n .ge. ndebug) then
      j=jdebug
      print *,'1',j,eos(ieos(1,j+1),6),V(n+2,j+1),V(n,j+1),V(n+1,j+1)
      print *,'s',s1(n+2,j+1),s2(n+2,j+1),s3(n+2,j+1),epsi1(n+1,j+1)
        endif
c      if (idebug .eq. 1 .and. n .ge. ndebug) read(*,*) foo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     6. VON MISES YIELD CONDITION - start j
c
      do j=0,jj-2,2
      if (ibc(j+1) .eq. 0) then
c
c This is yield condition
       k(n+2,j+1) =(s1(n+2,j+1)**2 + s2(n+2,j+1)**2 + s3(n+2,j+1)**2)
     &           -(2.d0/3.d0)*Y(j+1)**2
c       print *,'Yield check',Y(j+1),k(n+2,j+1)
c
       if (k(n+2,j+1) .gt. 0.) then   ! the material yielded!
c
c       Print *,'material yielded',j,k(n+2,j+1)
c        read(*,*) foo
c        print *,s1(n+2,j+1),s2(n+2,j+1),s3(n+2,j+1)
c         xxn = dsqrt(s1(n  ,j+1)**2 + s2(n  ,j+1)**2 + s3(n  ,j+1)**2)
c         xxn = dsqrt(2.d0/3.d0)*Y(j+1)/xxn  ! starting time step equivalent  stress
c
c VON MISES
c   xx is the scaling factor
c   epsip these are the strain rates calculated above
c   eqsn   is the equivalent plastic strain at the current time step n
c   eqsn2s is the equivalent plastic strain at intermediate future time step (n+1)*
c   eqsn2  is the equivalent plastic strain at the next time step (n+1)   
c   eestrain is the total elastic strain calculated from equivalent stress
c   epstrain is the total plastic strain
c
      if (istr(j+1) .eq. 1) then ! this is Von Mises
         ipstrain = 0.d0  !incremental plastic strain
         xx  = dsqrt(s1(n+2,j+1)**2 + s2(n+2,j+1)**2 + s3(n+2,j+1)**2) ! scale factor eq 3.17
         xx  = Y(j+1)/(dsqrt(3.d0/2.d0)*xx)
c
         eqsn  =dsqrt(3.d0/2.d0)*     
     &          dsqrt(s1(n  ,j+1)**2 + s2(n  ,j+1)**2 + s3(n  ,j+1)**2) !incremental equivalent plastic stress at n
         eqsn2s=dsqrt(3.d0/2.d0)*     
     &          dsqrt(s1(n+2,j+1)**2 + s2(n+2,j+1)**2 + s3(n+2,j+1)**2) !incremental equivalent plastic stress at (n+1)*
c
         s1(n+2,j+1) =   xx*s1(n+2,j+1)  !by multiplying by xx we put the stress back on yield surface
         s2(n+2,j+1) =   xx*s2(n+2,j+1)
         s3(n+2,j+1) =   xx*s3(n+2,j+1)
c
         s1(n+1,j+1) = ( s1(n+2,j+1) + s1(n,j+1) )/2.d0
         s2(n+1,j+1) = ( s2(n+2,j+1) + s2(n,j+1) )/2.d0
         s3(n+1,j+1) = ( s3(n+2,j+1) + s3(n,j+1) )/2.d0
c
         epsip(n+2,1,j+1)=epsip(n,1,j+1)+(1.d0/xx-1.d0)*s1(n+2,j+1) ! incremental plastic strain in x-dir
     &                    /(2.d0*mu(j+1))
         epsip(n+2,2,j+1)=epsip(n,2,j+1)+(1.d0/xx-1.d0)*s2(n+2,j+1) ! incremental plastic strain in y-dir
     &                    /(2.d0*mu(j+1))
         epsip(n+2,3,j+1)=epsip(n,3,j+1)+(1.d0/xx-1.d0)*s3(n+2,j+1) ! incremental plastic strain in z-dir
     &                    /(2.d0*mu(j+1))
         eqsn2 =dsqrt(3.d0/2.d0)*     
     &          dsqrt(s1(n+2,j+1)**2 + s2(n+2,j+1)**2 + s3(n+2,j+1)**2) !incremental equivalent plastic stress at (n+1)
         eqstress(n+2,j+1)=eqsn2   ! equivalent plastic stress
         ipstrain = (1.d0/xx-1.d0)*eqsn2/(3.d0*mu(j+1)) !incremental plastic strain
c
        eqsn  = dsqrt(3.d0/2.d0)*
     &          dsqrt(s1(n  ,j+1)**2 + s2(n  ,j+1)**2 + s3(n  ,j+1)**2) !incremental equivalent plastic stress at n
        eqsn2 = dsqrt(3.d0/2.d0)*
     &          dsqrt(s1(n+2,j+1)**2 + s2(n+2,j+1)**2 + s3(n+2,j+1)**2) !incremental equivalent plastic stress at (n+1)
        epstrain(n+2,j+1)=epstrain(n,j+1) + ipstrain
c
       if (eestrain(n+2,j+1)+epstrain(n+2,j+1) .gt. epsilmax(j+1) ) then
         epsilmax(j+1) = eestrain(n+2,j+1)+epstrain(n+2,j+1)
       endif
c
c==========================================================================
       elseif (istr(j+1) .eq. 2) then ! this is SteinbergGund from Wilkins book
c
         ipstrain = 0.d0  !incremental plastic strain
         jconv    = -1    !counts how many iterations it takes to converge the plastic deformation 
c         Ynew     = Y(j+1)
c bisection method to converge yield and strain 
         Ytry     = Y(j+1) !+ 0.003
c
 201      jconv = jconv + 1
          k(n+2,j+1) =(s1(n+2,j+1)**2 + s2(n+2,j+1)**2 + s3(n+2,j+1)**2)
     &           -(2.d0/3.d0)*Ytry**2
          xx  = dsqrt(s1(n+2,j+1)**2 + s2(n+2,j+1)**2 + s3(n+2,j+1)**2) ! scale factor eq 3.17
          xx  = Ytry/(dsqrt(3.d0/2.d0)*xx)
c
          eqsn2 =dsqrt(3.d0/2.d0) * dsqrt((xx*s1(n+2,j+1))**2 
     &           + (xx*s2(n+2,j+1))**2 + (xx*s3(n+2,j+1))**2) !equivalent stress at (n+1), on yield surface
          ipstrain = (1.d0/xx-1.d0)*eqsn2/(3.d0*mu(j+1)) !incremental plastic strain
c
c   Steinberg-Guinan model to increase yield strength, Y(j+1)
      mu(j+1)=
     &     eos(ieos(1,j+1),6) ! this is shear strength in ko.in
     &     *(1.d0 + SGLb(j+1)*P(n,j+1)*V(n,j+1)**(1.d0/3.d0))
c     &      *(1.d0 +SGLb(j+1)* (ipstrain/deltat) ) 
      Ynew =
c     &       Y(j+1)
     &      eos(ieos(1,j+1),5)  ! this is yield strength in ko.in
     &  *(1.d0+SGLbeta(j+1)*(epstrain(n,j+1)+ipstrain))**SGLn(j+1) 
     &        *(1.d0 + SGLb(j+1)*P(n,j+1)*V(n,j+1)**(1.d0/3.d0))
c
         if (j .eq. 1266) then
           print *,'jconv',j,jconv,Ynew,Ytry
         endif
         if (dabs(Ytry- Ynew) .lt. 1d-4) then
          else
c            if (j+1 .eq. 1266) then
c             print *,'jconv',jconv,Ynew,Ytry
c            endif
            Ytry = Ytry + (Ynew-Ytry)*1.0d-4
c            print *,'jconv',jconv,Ynew,Ytry,
c     &                          SGLbeta(j+1),SGLn(j+1),SGLb(j+1)
            goto 201

         endif
c         enddo ! this is the convergence loop for JC yield 
c
c now that the yield surface and strains are converged, update the deviatoric stresses
c
        Y(j+1)      = Ytry
        epstrain(n+2,j+1)=epstrain(n,j+1) + ipstrain
        s1(n+2,j+1) =   xx*s1(n+2,j+1)  !by multiplying by xx we put the stress back on yield surface
        s2(n+2,j+1) =   xx*s2(n+2,j+1)
        s3(n+2,j+1) =   xx*s3(n+2,j+1)
c
        s1(n+1,j+1) = ( s1(n+2,j+1) + s1(n,j+1) )/2.d0
        s2(n+1,j+1) = ( s2(n+2,j+1) + s2(n,j+1) )/2.d0
        s3(n+1,j+1) = ( s3(n+2,j+1) + s3(n,j+1) )/2.d0
c 
         epsip(n+2,1,j+1)=epsip(n,1,j+1)+(1.d0/xx-1.d0)*xx*s1(n+2,j+1) ! incremental plastic strain in x-dir
     &                    /(2.d0*mu(j+1))
         epsip(n+2,2,j+1)=epsip(n,2,j+1)+(1.d0/xx-1.d0)*xx*s2(n+2,j+1) ! incremental plastic strain in y-dir
     &                    /(2.d0*mu(j+1))
         epsip(n+2,3,j+1)=epsip(n,3,j+1)+(1.d0/xx-1.d0)*xx*s3(n+2,j+1) ! incremental plastic strain in z-dir
     &                    /(2.d0*mu(j+1))
         eqstress(n+2,j+1)=eqsn2   ! equivalent plastic stress
c        print *,j+1,Y(j+1)
c
c===============================================================================
       elseif (istr(j+1) .eq. 3) then ! this is Johnson-Cook
c
         ipstrain = 0.d0  !incremental plastic strain
         jconv    = -1    !counts how many iterations it takes to converge the plastic deformation 
c         Ynew     = Y(j+1)
c bisection method 
         Ytry     = Y(j+1) !+ 0.003
c
c         do while ( dabs( Ynew - Ytry ) .gt. 1.d-6) 
c Ytry1
 301      jconv = jconv + 1
          k(n+2,j+1) =(s1(n+2,j+1)**2 + s2(n+2,j+1)**2 + s3(n+2,j+1)**2)
     &           -(2.d0/3.d0)*Ytry**2
          xx  = dsqrt(s1(n+2,j+1)**2 + s2(n+2,j+1)**2 + s3(n+2,j+1)**2) ! scale factor eq 3.17
          xx  = Ytry/(dsqrt(3.d0/2.d0)*xx)
c
          eqsn2    = dsqrt(3.d0/2.d0) * dsqrt((xx*s1(n+2,j+1))**2 
     &           + (xx*s2(n+2,j+1))**2 + (xx*s3(n+2,j+1))**2) ! equivalent stress at (n+1), on yield surface
          ipstrain = (1.d0/xx-1.d0)*eqsn2/(3.d0*mu(j+1))      ! incremental plastic strain
c
c   Johnson-Cook model to increase yield strength, Y(j+1)
c
      mu(j+1)=
c     &      mu(j+1)
     &     eos(ieos(1,j+1),6) ! this is shear strength in ko.in
     &      *(1.d0 +JOBJO(j+1)* (ipstrain/deltat) ) ! this did not seem to make much differemce 
      Ynew =
     &     Y(j+1)*0.d0 + JOAJO(j+1)
     &     + JOBJO(j+1)*(epstrain(n,j+1)*1.d0+ipstrain*1.d0)**JONJO(j+1)
     &          *(1.d0 + JOCJO(j+1)*(ipstrain/deltat) ) 
c     &        *(1.d0 - T-Tr/TM-Tr)  
c
c       if (j .eq. 340) then
c        print *,'Yield change',jconv,y(j+1),JOAJO(j+1),
c     &  epstrain(n,j+1)+ipstrain
c       endif
c        if (ncount .le. 10000 ) then !.and. j+1 .ge. 1266) then 
c         print *,'JO parameters: ',j+1,Y(j+1),JOAJO(j+1),JOBJO(j+1),
c     &   JOCJO(j+1),JOMJO(j+1),JONJO(j+1),JOTJO(j+1),JOPOS(j+1),
c     &   epstrain(n,j+1),ipstrain
c        endif
          if (dabs(Ytry- Ynew) .lt. 1d-4) then
          else
            Ytry = Ytry + (Ynew-Ytry)*1.0d-4
c        print *,'jconv',jconv,Ynew,Ytry,JOBJO(j+1),ipstrain,JONJO(j+1)  
            goto 301
          endif
c         enddo ! this is the convergence loop for JC yield 
c
c now that the yield surface and strains are converged, update the deviatoric stresses
c       
        Y(j+1)      = Ytry
c        print *,jconv,j+1,y(j+1),JOBJO(j+1),
c     &  epstrain(n,j+1),ipstrain,JONJO(j+1)
c
        epstrain(n+2,j+1)=epstrain(n,j+1) + ipstrain
c        print *,'ytry',j,jconv,Y(j+1)
c       print *,'Sum Stress',dsqrt(3.d0/2.d0)*
c     & dsqrt(s1(n+2,j+1)**2+s2(n+2,j+1)**2+s3(n+2,j+1)**2),
c     & Y(j+1),eqsn2, dsqrt(3.d0/2.d0)*
c     & dsqrt(s1(n+2,j+1)**2+s2(n+2,j+1)**2+s3(n+2,j+1)**2)-
c     & Y(j+1)-eqsn2 
         s1(n+2,j+1) =   xx*s1(n+2,j+1)  !by multiplying by xx we put the stress back on yield surface
         s2(n+2,j+1) =   xx*s2(n+2,j+1)
         s3(n+2,j+1) =   xx*s3(n+2,j+1)
c
         s1(n+1,j+1) = ( s1(n+2,j+1) + s1(n,j+1) )/2.d0
         s2(n+1,j+1) = ( s2(n+2,j+1) + s2(n,j+1) )/2.d0
         s3(n+1,j+1) = ( s3(n+2,j+1) + s3(n,j+1) )/2.d0
c 
         epsip(n+2,1,j+1)=epsip(n,1,j+1)+(1.d0/xx-1.d0)*xx*s1(n+2,j+1) ! incremental plastic strain in x-dir
     &                    /(2.d0*mu(j+1))
         epsip(n+2,2,j+1)=epsip(n,2,j+1)+(1.d0/xx-1.d0)*xx*s2(n+2,j+1) ! incremental plastic strain in y-dir
     &                    /(2.d0*mu(j+1))
         epsip(n+2,3,j+1)=epsip(n,3,j+1)+(1.d0/xx-1.d0)*xx*s3(n+2,j+1) ! incremental plastic strain in z-dir
     &                    /(2.d0*mu(j+1))
         eqstress(n+2,j+1)=eqsn2   ! equivalent plastic stress
c==========================================================================
        else
         print *,'Constiutive input error in istr'
         read (*,*) foo
c
c==========================================================================
       endif ! this endif is for VonMises vs. JC or SGL
c
c
c  Accumulation of Elastic Strain, which happens before yield
c
        eestrain(n+2,j+1)=eestrain(n,j+1) +
     &                     (eqsn2-eqsn) /(3.d0*mu(j+1)) ! total strain, see equation after 3.26c
c
        endif !(ibc(j+1) .eq. 0) then
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     7. ARTIFICIAL VISCOSITY - within j
c
       if (1 .eq. 1) then
c       print *,j,u(n+1,j+2),u(n+1,j)
       if( u(n+1,j+2) .lt. u(n+1,j) .and.
     &     V(n+2,j+1)-V(n,j+1) .lt. 0. ) then
        xx = rho(j+1)/V(n+1,j+1)
        if (p(n,j+1) .lt. 0.) then
c         print *,'Negative Pressure!'
         a  = dsqrt(-P(n,j+1)/xx)
        else
         a  = dsqrt( P(n,j+1)/xx)
        endif
        if (dabs(P(n,j+1)) .lt. 0.002d0) a = 0.d0
c       print *,'Sound Speed',j,P(n,j+1),xx,a
        q(n+1,j+1) = Co*Co*xx*    (U(n+1,j+2)-U(n+1,j))**2
     &             + CL*a *xx*dabs(U(n+1,j+2)-U(n+1,j))
       else
        q(n+1,j+1)=0.d0
       endif
      endif
      endif
      enddo !j counter
        if (idebug .eq. 1 .and. n .ge. ndebug) then
         j=jdebug
       print *,'Q',j,q(n+1,j+1),xx,(U(n+1,j+2)-U(n+1,j)),
     &    a,rho(j+1),v(n,j+1),P(n,j+1)
        endif
             if(idebug .eq.1 .and. n .ge. ndebug) read(*,*) foo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     8.  ENERGY  - start j
c     ieos(2,imat)=1 - Mie Grunisen
c     ieos(2,imat)=2 - Gamma law ideal gas
c     ieos(2,imat)=3 - Gamma law ideal gas with Newtonian Stress and heat transfer
c     ieos(2,imat)=4 - snow plow model with KO inputs for Hugoniot
c     ieos(2,imat)=5 - snow plow model with anamolous hugoniot
c     ieos(2,imat)=6 - P-alpha Model
c
      do j=0,jj-2,2
      if (ibc(j+1) .eq. 0) then
c
       if (ieos(2,j+1) .eq. 1) then                    !Mie Grunisen
        gamma0 = eos(ieos(1,j+1),4)
        k1  = rho(j+1)*eos(ieos(1,j+1),1)**2
        k2  =  k1*(2.d0*eos(ieos(1,j+1),2)-gamma0/2.d0)
        k3  =  k1*(3.d0*eos(ieos(1,j+1),2)-gamma0)*eos(ieos(1,j+1),2)
        qbar = (q(n+1,j+1)+q(n-1,j+1))/2.d0
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)
     &             +(d-1.d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*deltat
        strain = 1.d0 - V(n+2,j+1)
        xb = gamma0 !gamma
        if(strain .lt. 0.d0) then
         xa = k1*strain + k3*strain**3
        else
         xa = k1*strain + k2*strain**2 + k3*strain**3
        endif
        E(n+2,j+1)=
     &  ( E(n,j+1)-((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1))
     &            + deltaZ )
     &   / (1.d0 + xb*(V(n+2,j+1)-V(n,j+1))/2.d0)
        E(n+1,j+1)= ( E(n+2,j+1) + E(n,j+1) )/2.d0
c      if (n .eq. 33) then
c        write (*,'(a2,I3,9e9.2)')
c     &   'E ',j,E(n+2,j+1),E(n+2,j+1),xa,P(n,j+1),qbar,deltaZ
c     &   ,V(n+2,j+1),V(n,j+1)
cc      read(*,*) foo
c      endif
c
       elseif (ieos(2,j+1) .eq. 2) then                !Gamma law ideal gas
        qbar   = (q(n+1,j+1)+q(n-1,j+1))/2.d0
        delt   = t(n+1)-t(n)
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)
     &              +(d-1.d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*delt
c
        xa   = 0.d0
        xb   = (eos(ieos(1,j+1),4)-1.d0)/V(n+2,j+1)  !gamma
        E(n+2,j+1)= ( E(n,j+1)
     &   -((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1)) + deltaZ )/
     &  (1.d0+xb*(V(n+2,j+1)-V(n,j+1))/2.d0)
        E(n+1,j+1)= ( E(n+2,j+1) + E(n,j+1) )/2.d0
c      if (j .eq. 0 ) then
c        if (n .lt.2. .or. n .eq. 201)    E(n+2,j+1)=0.d0
c        xx=V(n+2,j+1)-V(n,j+1)
c        print *,E(n,j+1)   /(1.*(1.+xb*(V(n+2,j+1)-V(n,j+1))/2.))
c        print *,xa*xx      /(2.*(1.+xb*(V(n+2,j+1)-V(n,j+1))/2.))
c        print *,P(n,j+1)*xx/(2.*(1.+xb*(V(n+2,j+1)-V(n,j+1))/2.))
c        print *,qbar*xx    /(1.*(1.+xb*(V(n+2,j+1)-V(n,j+1))/2.))
c        print *,deltaZ  /(1.*(1.+xb*(V(n+2,j+1)-V(n,j+1))/2.))
c      endif
c
       elseif (ieos(2,j+1) .eq. 3) then                !Gamma law ideal gas     ???
        qbar = (q(n+1,j+1)+q(n-1,j+1))/2.d0
        delt   = t(n+1)-t(n)
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)
     &              +(d-1.0d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*delt
c
        xx   = rho(j+1)/V(n+2,j+1)
        xa   = 0.d0
        xb   = (eos(ieos(1,j+1),4)-1.d0)/V(n+2,j+1)
        E(n+2,j+1)= ( E(n,j+1)
     &   -((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1))+ deltaZ )/
     &  (1.d0+xb*(V(n+2,j+1)-V(n,j+1))/2.d0)
c
       elseif (ieos(2,j+1) .eq. 4) then                !Snow Plow
        qbar = (q(n+1,j+1)+q(n-1,j+1))/2.d0
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)
     &             +(d-1.d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*deltat
        xb = eos(ieos(1,j+1),4) !gamma
        v0 = eos(ieos(1,j+1),12)
        xx = 1.d0 - (V(n+2,j+1)/rho(j+1))/V0   !this makes the compression relative to the compacted material hugoniot
        if ( v(n+2,j+1)/rho(j+1) .le. v0) icompact(n+2,j+1) = 1
        if(     icompact(n+2,j+1) .eq. 1 .and. xx .lt. 0.d0) then
         xa = eos(ieos(1,j+1),1)*xx+eos(ieos(1,j+1),3)*xx**3
        elseif (icompact(n+2,j+1) .eq. 1 .and. xx .ge. 0.d0) then
         xa = eos(ieos(1,j+1),1)*xx+eos(ieos(1,j+1),2)*xx**2
     &       +eos(ieos(1,j+1),3)*xx**3
        elseif( icompact(n+2,j+1) .eq. 0) then
         xa = 0.d0   !snow plow model - this zeros out pressure until threshold is reached
        else
         print *,'compaction error in energy'
         print *,icompact(n+2,j+1),xx
         read(*,*) foo
        endif
        E(n+2,j+1)=
     &  ( E(n,j+1)-((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1))
     &            + deltaZ )
     &   / (1.d0 + xb*(V(n+2,j+1)-V(n,j+1))/2.d0)
        E(n+1,j+1)= ( E(n+2,j+1) + E(n,j+1) )/2.d0
c
       elseif (ieos(2,j+1) .eq. 5) then                !Snow Plow with anaomolus hugoniot
        qbar = (q(n+1,j+1)+q(n-1,j+1))/2.d0
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)
     &             +(d-1.d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*deltat
        stemp = 0.1048d0      !slope
        ctemp = 0.5124d0      !bulk sound speed
        gtemp = 0.9d0         !gamma
        v0 = eos(ieos(1,j+1),8)
        V00= 1.d0/rho(j+1)
        vv = V(n+2,j+1)/rho(j+1)    !v_local
        if ( v(n+2,j+1)/rho(j+1) .le. v0) icompact(n+2,j+1) = 1
        if (icompact(n+2,j+1) .eq. 1) then  !icompact is the flag, once compact always compact
         xa = ((2.d0*vv-gtemp*(v0 -vv))*ctemp*ctemp*(v0-vv))     !porous hugoniot, see meyers pg 141
     &       /((2.d0*vv-gtemp*(v00-vv))*((v0-stemp*(v0-vv))**2) )
         xa = dabs(xa) !not sure if i need this but it was add so that re-denstended materials don't have negative pressure
        else
         xa = 0.d0
        endif
        E(n+2,j+1)=
     &  ( E(n,j+1)-((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1))
     &            + deltaZ )
     &   / (1.d0 + xb*(V(n+2,j+1)-V(n,j+1))/2.d0)
        E(n+1,j+1)= ( E(n+2,j+1) + E(n,j+1) )/2.d0
c
       elseif (ieos(2,j+1) .eq. 6) then                 !p-alpha
        qbar = (q(n+1,j+1)+q(n-1,j+1))/2.d0
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)
     &             +(d-1.d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*deltat
        xx = 1.d0 - V(n+2,j+1)
        xb = eos(ieos(1,j+1),4) !gamma
        if(xx .lt. 0.) then
        xa = eos(ieos(1,j+1),1)*xx+0.d0*xx**2+eos(ieos(1,j+1),3)*xx**3
        else
        xa = eos(ieos(1,j+1),1)*xx+eos(ieos(1,j+1),2)*xx**2
     &      +eos(ieos(1,j+1),3)*xx**3
        endif
        E(n+2,j+1)=
     &  ( E(n,j+1)-((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1))
     &            + deltaZ )
     &   / (1.d0 + xb*(V(n+2,j+1)-V(n,j+1))/2.d0)
        E(n+1,j+1)= ( E(n+2,j+1) + E(n,j+1) )/2.d0
c
c
       else
        print *,'EOS error: n,j,ibc(j+1),ieos(2,j+1),'
        print *,n,j,ibc(j-1),ieos(2,j-1),ibc(j+1),ieos(2,j+1)
        read(*,*) foo
       endif
      endif
c
      enddo
c
      if(idebug .eq.1  .and. n .ge. ndebug) then
        j=jdebug
        print *,'Energy',j,E(n+2,j+1),E(n,j+1),xx,
     & (2.d0*(1.d0+xb*(V(n+2,j+1)-V(n,j+1))/2.d0))
      endif
      if(idebug .eq.1  .and. n .gt. ndebug) read(*,*) foo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     5b PRESSURE - start j
c
c      print *,'line 935'

      do j = 0,jj-2,2
      if (ibc(j+1) .eq. 0) then
       if (ieos(2,j+1) .eq. 1) then   !     Mie Gruneisen EOS as formulated in Wilkins
        gamma0 = eos(ieos(1,j+1),4)
        k1  = rho(j+1)*eos(ieos(1,j+1),1)**2
        k2  =  k1*(2.d0*eos(ieos(1,j+1),2)-gamma0/2.d0)
        k3  =  k1*(3.d0*eos(ieos(1,j+1),2)-gamma0)*eos(ieos(1,j+1),2)
c        if (j .eq. 2840) then
c        print *,eos(ieos(1,j+1),1),eos(ieos(1,j+1),2),eos(ieos(1,j+1),3)
c        read(*,*) foo
c        endif
        strain      = 1.d0 - V(n+1,j+1)
        if (strain .lt. 0.d0) then
         P(n+1,j+1) = k1*strain + k3*strain**3 + gamma0*E(n+1,j+1)
        else
         P(n+1,j+1) = k1*strain + k2*strain**2 + k3*strain**3
     &                                         + gamma0*E(n+1,j+1)
        endif
        strain      = 1.d0 - V(n+2,j+1)
        if (strain .lt. 0.d0) then
         P(n+2,j+1) = k1*strain + k3*strain**3 + gamma0*E(n+2,j+1)
        else
         P(n+2,j+1) = k1*strain + k2*strain**2 + k3*strain**3
     &                                         + gamma0*E(n+2,j+1)
        endif
c        print *,'Mie Grun EOS 1',P(n+1,j+1),P(n+2,j+1)
c        read(*,*) foo
c       if (icontact .eq. 1) then
c       if (P(n+2,j+1) .lt. 0.) then
c       print *,'Neg Pres',n+1,j+1,xx,v(n+1,j+1),E(n+1,j+1),p(n+1,j+1)
c       print *,'Neg Pres',n+2,j+1,xx,v(n+2,j+1),E(n+2,j+1),p(n+2,j+1)
cc       read(*,*) foo
c       endif
c       endif
       elseif (ieos(2,j+1) .eq. 2) then            ! Gamma Law (perfect gas)
         P(n+2,j+1) = (eos(ieos(1,j+1),4)-1.d0)*E(n+2,j+1)/V(n+2,j+1)
         P(n+1,j+1) = (P(n,j+1)+P(n+2,j+1))/2.d0
       elseif (ieos(2,j+1) .eq. 3) then            ! Gamma Law
         P(n+2,j+1) = (eos(ieos(1,j+1),4)-1.d0)*E(n+2,j+1)/V(n+2,j+1)
         P(n+1,j+1) = (P(n,j+1)+P(n+2,j+1))/2.d0
c        xa          = (u(n+2,j)+u(n+2,j+2))/2.d0
c        xx          = rho(j+1)/V(n+2,j+1)
c         P(n+2,j+1) = P(n,j+1)
c     &     +(eos(ieos(1,j+1),4)+1.d0)*xx*xa*xa/2.d0
c       print *,'pres',rho_0,gamma+1,xa,xx
       elseif (ieos(2,j+1) .eq. 4) then             !Snow Plow
        v0  = eos(ieos(1,j+1),12)
        V00 = 1.d0/rho(j+1)
        vv  = v(n+1,j+1)/rho(j+1)
        if (v0 .eq. 0.) then
         print *,'v0 can not equal zero'
         print *,j,ieos(1,j+1),eos(ieos(1,j+1),12)
         read(*,*) foo
        endif
c
        if (vv .le. V0) icompact(n+1,j+1) = 1
         xx         = 1.d0 - (V(n+1,j+1)/rho(j+1))/V0   !this makes the compression relative to the compacted material hugoniot
        if (icompact(n+1,j+1) .eq. 1 .and. xx .lt. 0.) then
          P(n+1,j+1) = eos(ieos(1,j+1),1)*xx+
     &                 eos(ieos(1,j+1),3)*xx**3+
     &                 eos(ieos(1,j+1),4)*E(n+1,j+1)
        elseif(icompact(n+1,j+1) .eq. 1 .and. xx .ge. 0.d0) then
          P(n+1,j+1) = eos(ieos(1,j+1),1)*xx+
     &                 eos(ieos(1,j+1),2)*xx**2+
     &                 eos(ieos(1,j+1),3)*xx**3+
     &                 eos(ieos(1,j+1),4)*E(n+1,j+1)
        elseif (icompact(n+1,j+1) .eq. 0) then
          P(n+1,j+1) = 0.d0*eos(ieos(1,j+1),4)*E(n+1,j+1)
        else
         Print *,'pressure compaction error'
         read(*,*) foo
        endif
c        if ( icompact(n+1,j+1) .eq. 1) then
c        print *,t(n),j+1,icompact(n+1,j+1),vv,p(n+1,j+1)
c        print *,eos(ieos(1,j+1),1)*xx,eos(ieos(1,j+1),3)*xx**3
c     &         ,eos(ieos(1,j+1),4)*E(n+1,j+1)
c        print *,xx,eos(ieos(1,j+1),1),eos(ieos(1,j+1),3)
c     &         ,eos(ieos(1,j+1),4)
c        read(*,*) foo
c        else
c        print *,t(n),j+1,icompact(n+1,j+1),vv,p(n+1,j+1)
c        endif
c
        if (v(n+2,j+1)/rho(j+1) .le. V0) icompact(n+2,j+1) = 1
         xx         = 1.d0 - (V(n+2,j+1)/rho(j+1))/V0     !this makes the compression relative to the compacted material hugoniot
         if (icompact(n+2,j+1) .eq. 1 .and. xx .lt. 0.) then
          P(n+2,j+1) = eos(ieos(1,j+1),1)*xx+
     &                 eos(ieos(1,j+1),3)*xx**3+
     &                 eos(ieos(1,j+1),4)*E(n+2,j+1)
         elseif (icompact(n+2,j+1) .eq. 1 .and. xx .ge. 0.) then
          P(n+2,j+1) = eos(ieos(1,j+1),1)*xx+
     &                 eos(ieos(1,j+1),2)*xx**2+
     &                 eos(ieos(1,j+1),3)*xx**3+
     &                 eos(ieos(1,j+1),4)*E(n+2,j+1)
        elseif (icompact(n+2,j+1) .eq. 0) then  !v .gt. v0
          P(n+2,j+1) = 0.d0*eos(ieos(1,j+1),4)*E(n+2,j+1)
        else
         Print *,'pressure compaction error'
         read(*,*) foo
        endif
c
       elseif (ieos(2,j+1) .eq. 5) then   ! Snow Plow with anaomolous Hugoniot
         V0 = eos(ieos(1,j+1),8)
         V00= 1.d0/rho(j+1)

         vv = V(n+1,j+1)/rho(j+1)    !v_local
         if (VV .le. V0) icompact(n+1,j+1) = 1
         if (icompact(n+1,j+1) .eq. 1) then
          xa = ((2.d0*vv-gtemp*(v0 -vv))*ctemp*ctemp*(v0-vv))    !porous hugoniot, see meyers pg 141
     &        /((2.d0*vv-gtemp*(v00-vv))*((v0-stemp*(v0-vv))**2) )
          if (xa .lt. 0.d0) then
           print *,'Compression less than v0 on anamolous hugoniot, n+1'
           xa=dabs(xa)
          endif
         else
          xa = 0.d0
         endif

         P(n+1,j+1) = xa + eos(ieos(1,j+1),4)*E(n+1,j+1)

         vv = V(n+2,j+1)/rho(j+1)    !v_local
         if (VV .le. V0) icompact(n+2,j+1) = 1
         if (icompact(n+2,j+1) .eq. 1) then
         xa = ((2.d0*vv-gtemp*(v0 -vv))*ctemp*ctemp*(v0-vv))     !porous hugoniot, see meyers pg 141
     &       /((2.d0*vv-gtemp*(v00-vv))*((v0-stemp*(v0-vv))**2) )
          if (xa .lt. 0.d0) then
           print *,'Compression less than v0 on anamolous hugoniot, n+2'
           xa=dabs(xa)
          endif
         else
          xa = 0.d0
         endif
c
       elseif (ieos(2,j+1) .eq. 6) then  !p-alpha model
c
c
       elseif  (ieos(2,j+1) .eq. 7) then  ! Mie Gruneisen CTH like formulation
        P0 = 0.d0
        E0 = 0.d0
        T0 = 293.d0
        strain        = 1.d0 - V(n+1,j+1)
        rho_local = rho(j+1)/V(n+1,j+1)
        up = (U(n+1,j)+U(n+1,j+2))/2.d0
        Us = eos(ieos(1,j+1),1)           ! local shock speed
     &     + eos(ieos(1,j+1),2)*up
     &     + eos(ieos(1,j+1),3)/eos(ieos(1,j+1),1)*up**2
        PH = P0 + rho_local*Us*up
        EH = E0 + up*up/2.d0
        TH = dexp( eos(ieos(1,j+1),4)*strain)*T0
!       E  = EH + eos(ieos(1,j+1),8)*(T-TH)
!       P(n+2,j+1) = PH + eos(ieos(1,j+1),4)*rho_local*(
c         if (Us*up .ne. Us*Us*strain) then
c          print *,'EOS does not jive'
c          read(*,*) foo
c         endif
        xx         = 1.d0 - V(n+2,j+1)
        P(n+2,j+1) = eos(ieos(1,j+1),1)*xx+
     &               eos(ieos(1,j+1),2)*xx**2+
     &               eos(ieos(1,j+1),3)*xx**3+
     &               eos(ieos(1,j+1),4)*E(n+2,j+1)
       else    !this is ieos(2,j+1) else
      Print *,'eos error in setup file for node'
     &         ,j+1,ieos(1,j+1),ieos(2,j+1)
        read(*,*) foo
       endif  ! IEOS
      endif !ibc
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This is the Temperature and Entropy Calculation
c      print *,'line 1104'
      do j=0,jj-2,2
c	print *,'line 1106',j,nc
c	print *,'temp',temp(n+2c,j+1)
c	print *,'E   ',e(n+2,j+1)
c	print *,'rho ',rho(j+1)
c	print *,'ieos',ieos(1,j+1)
c	print *,'eos ',eos(ieos(1,j+1),8)
c      
       xx  = 1.d0 - V(n+2,j+1)
       e00 = -300.d0*eos(ieos(1,j+1),8)
       e01 =   e00 * eos(ieos(1,j+1),4)
       e02 = (eos(ieos(1,j+1),1)**2+eos(ieos(1,j+1),4)**2.d0*e00) /2.d0
       e03 = (4.d0*eos(ieos(1,j+1),2)*eos(ieos(1,j+1),1)**2 +
     &             eos(ieos(1,j+1),4)**3*e00) / 6.d0
       e04 = (18.d0*eos(ieos(1,j+1),2)**2*eos(ieos(1,j+1),1)**2 - 
     &         2.d0*eos(ieos(1,j+1),4)*eos(ieos(1,j+1),2)*
     &              eos(ieos(1,j+1),1)**2 +
     &         eos(ieos(1,j+1),4)**4*e00 ) / 24.d0
       e0 = e00 + e01*xx +  e02*xx**2 + e03*xx**3 + e04*xx**4
       p0 = rho(j+1)*(e01+2.d0*e02*xx+3.d0*e03*xx**2+4.d0*e04*xx**3)
       Temp(n+2,j+1)=(E(n+2,j+1)/rho(j+1)-e0)/eos(ieos(1,j+1),8)
c       print *,'e0',Temp(n+2,j+1),e0,E(n+2,j+1)/rho(j+1)
c       print *,'p0',p0,v(n+2,j+1),rho(j+1)
c       print *,'T0',Temp(n+2,j+1)
c     &             ,(E(n+2,j+1)/rho(j+1))/eos(ieos(1,j+1),8)
c     &             ,-e0/eos(ieos(1,j+1),8)
c       if (rho(j+1) .gt. 0.d0) then
c         print *,Temp(n+2,j+1),e0,E(n+2,j+1)/rho(j+1),xx
c       endif
c       print  *,'t',eos(ieos(1,j+1),1),eos(ieos(1,j+1),2),
c     & eos(ieos(1,j+1),3),eos(ieos(1,j+1),4),
c     & eos(ieos(1,j+1),5),eos(ieos(1,j+1),6)

c     &  ,eos(ieos(1,j+1),8)
c      print *,'line 1109',j,n
       entropy(n+2,j+1)=entropy(n,j+1)+
     &  s1(n+2,j+1)/Temp(n+2,j+1)
c       entropy(n+2,j+1)=entropy(n,j+1)+
c       print *,entropy(n+2,j+1),entropy(n,j+1),
c     &  s1(n+2,j+1),epsi1(n+2,j+1),Temp(n+2,j+1)
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This is the Pressure and Energy convergence check
c     print *,'line 1118'
      do j=0,jj-2,2
       if (ieos(2,j+1) .eq. 1) then                    !Mie Grunisen
        gamma0 = eos(ieos(1,j+1),4)
        k1  = rho(j+1)*eos(ieos(1,j+1),1)**2
        k2  =  k1*(2.d0*eos(ieos(1,j+1),2)-gamma0/2.d0)
        k3  =  k1*(3.d0*eos(ieos(1,j+1),2)-gamma0)*eos(ieos(1,j+1),2)
        qbar = (q(n+1,j+1)+q(n-1,j+1))/2.d0
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)
     &             +(d-1.d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*deltat
        strain = 1.d0 - V(n+2,j+1)
        xb = gamma0 !gamma
        if(strain .lt. 0.d0) then
        xa = k1*strain + k3*strain**3
        else
        xa = k1*strain + k2*strain**2 + k3*strain**3
        endif
        En2j1=
     &  ( E(n,j+1)-((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1))
     &            + deltaZ )
     &   / (1.d0 + xb*(V(n+2,j+1)-V(n,j+1))/2.d0)
         diffE=    E(n+2,j+1)-En2j1
        if (dabs(diffE) .gt. 0.d-6 ) then
         print *,'Energy not converged',En2j1,E(n+2,j+1)
c         if (
        endif
       elseif (  ieos(2,j+1) .eq. 0) then
       !this is for not used nodes
       else
c        print *,'Non-MG EOS',j+1,ieos(2,j+1)
c        read(*,*) foo
       endif
      enddo
cccccccccccccccccccccccccc   Spall   ccccccccccccccccccccccccccccccccccccccccccc
c     Check for tensile fracture, ie does the pressure/stress exceed pfrac at j?
c     check for physical separation
c
      if (1 .eq. ipfrac) then
      do j=0,jj-2,2
       if (-p(n+2,j+1)+s1(n+2,j+1) .gt. pfrac(j)) then
c        print *,j,n,-p(n+2,j+1)+s1(n+2,j+1),pfrac(j)
        print *,'Tensile Soall Fracture at node: ',j
        print *,'Spall at position ',r(n,j),'cm at time ',t(n),'us'
c        print *,'Before Separating:'
c        print *,'The node at:'
c        print *,'left: r(j-2 ,n+2)=',r(n+2,j-2),' U(n-1)=',U(n+1,j-2)
c        print *,'J:    r(j   ,n+2)=',r(n+2,j)  ,' U(n-1)=',U(n+1,j+0)
c        print *,'right:r(j+2 ,n+2)=',r(n+2,j+2),' U(n-1)=',U(n+1,j+2)
c        print *,'      r(j+4 ,n+2)=',r(n+2,j+4),' U(n-1)=',U(n+1,j+4)
c       print *,'Frac',-p(n+2,j+1)+s1(n+2,j+1),pfrac(j)
        do jjj = jj-4,j+2,-2
         do nz = n,n+2
              ibc(    jjj-0)  =    ibc(    jjj-2)
              ibc(    jjj-1)  =    ibc(    jjj-3)
             ieos(1,  jjj-1)  =   ieos(1,  jjj-3)
             ieos(2,  jjj-1)  =   ieos(2,  jjj-3)
             ieos(1,  jjj-0)  =   ieos(1,  jjj-2)
             ieos(2,  jjj-0)  =   ieos(2,  jjj-2)
                m(    jjj-1)  =      m(    jjj-3)
              rho(    jjj-1)  =    rho(    jjj-3)

                r(nz ,jjj-0)  =     r(nz ,jjj-2)
                r(nz ,jjj-1)  =     r(nz ,jjj-3)
                U(nz ,jjj-0)  =     U(nz ,jjj-2)
              phi(nz ,jjj-0) =    phi(nz ,jjj-2)
             beta(nz ,jjj-0) =   beta(nz ,jjj-2)
           sigmar(nz ,jjj-1) = sigmar(nz ,jjj-3)
           sigmao(nz ,jjj-1) = sigmao(nz ,jjj-3)
                V(nz ,jjj-1) =      V(nz ,jjj-3)
               s1(nz ,jjj-1) =     s1(nz ,jjj-3)
               s2(nz ,jjj-1) =     s2(nz ,jjj-3)
               s3(nz ,jjj-1) =     s3(nz ,jjj-3)
                q(nz ,jjj-1) =      q(nz ,jjj-3)
            epsi1(nz ,jjj-1) =  epsi1(nz ,jjj-3)
            epsi2(nz ,jjj-1) =  epsi2(nz ,jjj-3)
                E(nz ,jjj-1) =      E(nz ,jjj-3)
                Y(    jjj-1) =      Y(    jjj-3)
                k(nz ,jjj-1) =      k(nz ,jjj-3)
            pfrac(    jjj-0) =  pfrac(    jjj-2)  !node value
            pfrac(    jjj-1) =  pfrac(    jjj-3)  !Cell value
           enddo
          enddo
c
          do nz = n,n+2  !void cell center
                 U(nz,j+1) = zero
                 U(nz,j+1) = zero
              phi(nz ,j+1) = zero
             beta(nz ,j+1) = zero
           sigmar(nz ,j+1) = zero
           sigmao(nz ,j+1) = zero
                V(nz ,j+1) = zero
               s1(nz ,j+1) = zero
               s2(nz ,j+1) = zero
               s3(nz ,j+1) = zero
                E(nz ,j+1) = zero
                Y(    j+1) = zero
                q(nz ,j+1) = zero
            pfrac(    j+1) = zero
            pfrac(    j+1) = zero
          enddo
                ibc(j)   = 2 !outer
                ibc(j+1) = 9 !void
                ibc(j+2) =-2 !inner
              P(n+1,j+1) = pvoid
              P(n+2,j+1) = pvoid
              P(n+1,j-1) = pvoid
              P(n+2,j-1) = pvoid
              P(n+1,j+3) = pvoid
              P(n+2,j+3) = pvoid
              u(n+1,j+2) = u(n+1,j+4)
              r(n+1,j+1) = r(n+1,j)
              r(n+2,j+1) = r(n+2,j)
              pfrac(j+2) =  pfrac(j+4)  !node value
              pfrac(j+3) =  pfrac(j+4)  !node value
c        print *,'After Separating:',n,j
c        print *,'r(j-2)=',r(n+2,j-2),'U= ',u(n+1,j-2),'P=',P(n+2,j-2)
c        print *,'r(j-1)=',r(n+2,j-1),'U= ',u(n+1,j-1),'P=',P(n+2,j-1)
c        print *,'r(j)  =',r(n+2,j  ),'U= ',u(n+1,j-0),'P=',P(n+2,j+0)
c        print *,'r(j+1)=',r(n+2,j+1),'U= ',u(n+1,j+1),'P=',P(n+2,j+1)
c        print *,'r(j+2)=',r(n+2,j+2),'U= ',u(n+1,j+2),'P=',P(n+2,j+2)
c        print *,'r(j+3)=',r(n+2,j+3),'U= ',u(n+1,j+3),'P=',P(n+2,j+3)
c        print *,'r(j+4)=',r(n+2,j+4),'U= ',u(n+1,j+4),'P=',P(n+2,j+4)
c               read(*,*) foo
c
c      do jzz=1,jj-2,1
cc       write(*,'(I4,5I2,10f7.3)')
c       write(*,'(I4,5I2,3f7.3,2e11.3,5f7.3)')
c     &  jzz,ieos(1,jzz),ieos(2,jzz),
c     &  ibc(jzz-1),ibc(jzz),ibc(jzz+1),
c     &  r(n+2,jzz-1),r(n+2,jzz),r(n+2,jzz+1),
c     &  U(n+1,jzz-1),U(n+1,jzz),U(n+1,jzz+1),
c     &  p(n+2,jzz),Temp(n+2,jzz)
cc     &  ,v(n+2,jzz+1),m(jzz+1),rho(jzz+1),y(jzz+1),pfrac(jzz+2)
c      enddo
c      read(*,*) foo

      ipfrac = 0
      goto 1999
      endif ! this is the "it fractured" if


        enddo ! this is the fracture check
 1999 endif ! this is the fracture on if
c
c************************   Spall End    **************************
ccccccccccccccccccccccccc   TIME-STEP    cccccccccccccccccccccccccc
c
c      print *,'line 1258'
       dt_min =  1.d9
       dr_min   =  r(n+2,0+2)-r(n+2,0  )  !this is just to get us started
c
       jzz    = -1
      do j = 0, jj-2, 2 ! loop through to find the fastest wave to set min time step
      if (ibc(j+1) .eq. 0) then
c       deltat =     t(n+2)    -t(n)
       Vdot   =   ( V(n+2,j+1)-V(n  ,j+1) )/deltat
       deltar = abs(r(n+2,j+2)-r(n+2,j  ) )
c
       b      = 8.d0*(Co**2+CL)*deltar*(Vdot/V(n+1,j+1))
       if (Vdot/V(n+1,j+1) .ge. 0.) b = 0.d0
c
       rho_local = rho(j+1)/V(n+2,j+1)
c       print *,'rho_local',j,rho_local,rho(j+1),1.d0/v(n+2,j+1)
       if (rho_local .gt. 1.e10) then
       print *,'Error: density too high'
       print *,j,ncount,t(n)
       print *,rho(j+1),v(n+2,j+1)
       print *,j,r(n+2,j+2),r(n+2,j)
       print *,V(n+2,j+1),V(n,j+1)
       print *,P(n+1,j),p(n+1,j+1),P(n+2,j),p(n+2,j+1)
       read(*,*) foo
       endif
c       print *,'rho local',rho_local
       a  = sqrt(abs(P(n+2,j+1))/ rho_local ) !I checked units
c       print *,'ts',j,a,b,deltar,Vdot,V(n+1,j+1)
       if ( (a**2+b**2) .ne. 0.d0 ) then
        delt_temp = (2.d0/3.d0)*(deltar/sqrt(a**2+b**2))
        if (delt_temp .lt. dt_min) then
         jzz    = j
         dt_min = delt_temp
         a_min = a
         b_min = b
         rho_min=rho_local
c        print *,'rr',rho_min,deltat,jzz,a_min,b_min,p(n+2,j+1)
c        read(*,*) foo
        endif
       endif
       endif !ibc(j+1)
      enddo !j loop

      dt_min=dt_min/20.d0    !Here is where you can arbitrarly lower time step
                            !the demomenator is suppose to be 1.d0

c      if (ncount .eq. 1) then
c              print *,dtmin
c              read(*,*)
c      endif  
      delt = t(n+1)-t(n-1)
      if (dt_min .gt. 1.1d0*delt ) then
        dt_min = 1.1d0 * delt
      endif
      deltat = (dt_min+delt)/2.d0 ! deltat is the full time step between t(n) t(n+2)
c
c      t(n+3) = t(n+2) + deltat           ! this is now advanced along with all the
c      t(n+4) = t(n+2) + deltat + deltat  ! other state variables
c
      if (idebug .eq. 1) then
      if (icontact .eq. 1) then
      Print *,'End of first contact'
             if(idebug .eq.1 )read(*,*) foo
      endif
      endif
c        if (deltat .lt. 0.1d-5) then  ! this is just larry swabies suggestion for stability
c        print *,rho_min,deltat,jzz,a_min,b_min
c        endif

c*******************OUTPUT SOLUTION to Screen**********************
c

      if (idebug .eq. 1 .and. ncount .gt. ndebug) then
      do j = 1,jj-2,2
        write(*,'(2I4,11f7.4)')
     & n+1,j,t(n),r(n+1,j-1),r(n+1,j+1),U(n+1,j-1),U(n+1,j+1)
     &,P(n+2,j),q(n+2,j),e(n+2,j),V(n+2,j)
      enddo
      read(*,*) foo
      endif
c
      jzz = ncount
      if (int((jzz-1)/iskip) .eq. dfloat(jzz-1)/dfloat(iskip) .or.
     &                t(n+2) .ge. tstop) then
c      if (n .eq. 1 .or. t(n) .ge. tskip ) then
      qtotal  = 0.d0  ! total artifical viscosity
      mvtotal = 0.d0  ! total momentum
      etotal  = 0.d0  ! total energy
      ketotal = 0.d0  ! total kinetic energy
      ietotal = 0.d0  ! total internal energy
      ke3total= 0.d0  ! total kinetic energy
c
      do j = 0,jj-2,2
       if (ibc(j+1) .eq. 0) then
       qtotal =qtotal  + q(n+2,j)  
       mvtotal=mvtotal + m(j+1)*(  u(n+1,j) + u(n+1,j+2) )/2.d0
       ketotal=ketotal + m(j+1)*(((u(n+1,j) + u(n+1,j+2))/2.d0)**2)/2.d0
       ietotal =ietotal  +(E(n+2,j+1)-E(n,j+1))/deltat
     &         + P(n+1,j+1)*(V(n+2,j+1)-V(n,j+1))/deltat
       etotal=ketotal-ietotal
       if (ieos(1,j-1) .eq. 2) then
       ke3total=ke3total + m(j+1)*(((u(n+1,j) + u(n+1,j+2))/2.d0)**2)/2.d0
       endif
       endif
      enddo
      if (ncount .eq. 1) then
      write(*,'(A8,5X,A8,8X,A8,8X,A8,8X,A8,8X,A8)')
     &   ' ncount ','Time ','mvtotal ','KEtotal ','ietotal ','etotal  '
      endif
        write(*,'(I8,6e15.6)')
     & ncount,t(n+1),mvtotal,KEtotal,ietotal,etotal,r(n,30)
c       ioffset = 10
c       do j=98-ioffset,98,2
c         write(*,99)
c     & n+2,j+1,t(n),r(n+1,j),r(n+1,j+2),U(n+1,j),U(n+1,j+2)
c     &  ,P(n+2,j+1),e(n+2,j+1),q(n+1,j+1)
c       enddo
c       do j=98,98+ioffset,2
c         write(*,99)
c     & n+2,j+1,t(n),r(n+1,j),r(n+1,j+2),U(n+1,j),U(n+1,j+2)
c     &  ,P(n+2,j+1),e(n+2,j+1),q(n+1,j+1)
c       enddo
c      print *, ' '
c
c      skip = (jj/2)/256
c      skip = nn/skip
c      do j=ioffset,nn,skip
c      write(*,98)
c     &'n','Cell','dt','rad','rad','Vel','Vel'
c     &,'pres','Vol'
c      print *,' '
c      enddo
c      write(*,'(2I4,6e9.2,6f7.4)')
c     & n,j,deltat,r(n,j-1),r(n,j+1),U(n,j-1),U(n,j+1),P(n,j)
c     & ,V(n,j),rho(j)
      endif
c
c ******************* Write solution to File  *************************
c       write(*,*) 'Enter file name'
c       read(*,*) foo name
c       name='al.txt'
c       open (unit=33, file=name, form='formatted', status='unknown')
c       do n=1,nn,2
        if (t(n) .eq. 0. .or. t(n) .ge. tskip ) then
          do j=0,jj-4,2
         bs1 =     V(n  ,j+1) !rho(j+1)
         bs2 =     U(n+1,j  )
         bs3 =     P(n  ,j+1)
         bs4 =     E(n  ,j+1)
         bs5 =     q(n+1,j+1)
         bs6 =    s1(n  ,j+1)
         bs7 =    s2(n  ,j+1)
         bs8 = epsi1(n-1,j+1)
         bs9 =  Temp(n  ,j+1)
 
         if (dabs(bs1)    .lt. 1.d-99)   bs1 = 0.d0  !this is because if the expoent is
         if (dabs(U(n+1,j  ))    .lt. 1.d-99)   bs2 = 0.d0  ! larger than 99 it will not write
         if (dabs(bs3)    .lt. 1.d-99)   bs3 = 0.d0  ! to the file correctly, so I just
         if (dabs(bs4)    .lt. 1.d-99)   bs4 = 0.d0  ! zero it out
         if (dabs(bs5)    .lt. 1.d-99)   bs5 = 0.d0
         if (dabs(bs6)    .lt. 1.d-99)   bs6 = 0.d0
         if (dabs(bs7)    .lt. 1.d-99)   bs7 = 0.d0
         if (dabs(bs8)    .lt. 1.d-99)   bs8 = 0.d0
         bs10 = -bs3 + bs6
         bs11 = -bs3 + bs7
         bs12= rho(j+1)/V(n+1,j+1)

c      epsi1(n+1,j+1) = (U(n+1,j+2)-U(n+1,j))/(r(n+1,j+2)-r(n+1,j))
c      epsi2(n+1,j+1)
c       if (j .eq. 3140) print *,'u',U(n+1,j)
c            write(33,'(3I6,2I3,30d26.18)')  ! "D" works with matlab not octave
            write(33,'(3I6,2I3,30e26.18)')  ! works with matlab and octave
     &     ncount,    j+1,       jj,  ibc(j),ibc(j+1),
     &     t(n)  , r(n,j), r(n,j+1),r(n,j+2),     bs1,
     &     bs2   ,    bs3, bs12    ,     bs4,     bs5,
     &     bs6   ,    bs7, k(n,j+1),ketotal,      bs9,
     &     bs10  ,    bs11, epsi1(n+1,j+1),epstrain(n,j+1),  
     &     epsip(n,1,j+1),epsip(n,2,j+1),epsip(n,3,j+1),
     &     eestrain(n+2,j+1),eestrain(n,j+1),Y(j+1),mu(j+1) ! plastic Strain
c            print *,'pressure',u(n+1,j),P(n,j+1),bs2,bs3
c             write(33,'(4d26.18)')
c     &   t(n),bs3,r(n,j),bs2
c
      enddo  ! this is j loop
            write(35,'(d26.18)') ke3total
      tskip = tskip + dtskip
      endif
c
cccccccccccccccccccccc Update solution ccccccccccccccccccccccccccccc
c
      do j=0,jj
      if (ibc(j) .ne. 9) then
c
        U(n-1,j)       = U(n+1,j)
        U(n  ,j)       = U(n+2,j)
        phi(n-1,j)     = phi(n+1,j)
        phi(n  ,j)     = phi(n+2,j)
        sigmar(n-1,j)  = sigmar(n+1,j)
        sigmar(n  ,j)  = sigmar(n+2,j)
        sigmao(n-1,j)  = sigmao(n+1,j)
        sigmao(n  ,j)  = sigmao(n+2,j)
        beta(n-1,j)    = beta(n+1,j)
        beta(n  ,j)    = beta(n+2,j)
        V(n-1,j)       = V(n+1,j)
        V(n  ,j)       = V(n+2,j)
        r(n-1,j)       = r(n+1,j)
        r(n  ,j)       = r(n+2,j)
        epsi1(n-1,j)   = epsi1(n+1,j)
        epsi1(n  ,j)   = epsi1(n+2,j)
        epsi2(n-1,j)   = epsi2(n+1,j)
        epsi2(n  ,j)   = epsi2(n+2,j)
        epsip(n-1,1,j) = epsip(n+1,1,j)
        epsip(n-1,2,j) = epsip(n+1,2,j)
        epsip(n-1,3,j) = epsip(n+1,3,j)
        epsip(n  ,1,j) = epsip(n+2,1,j)
        epsip(n  ,2,j) = epsip(n+2,2,j)
        epsip(n  ,3,j) = epsip(n+2,3,j)
        epstrain(n-1,j)= epstrain(n+1,j)
        epstrain(n  ,j)= epstrain(n+2,j)
        eestrain(n-1,j)= eestrain(n+1,j)
        eestrain(n  ,j)= eestrain(n+2,j)
        s1(n-1,j)      = s1(n+1,j)
        s1(n  ,j)      = s1(n+2,j)
        s2(n-1,j)      = s2(n+1,j)
        s2(n  ,j)      = s2(n+2,j)
        s3(n-1,j)      = s3(n+1,j)
        s3(n  ,j)      = s3(n+2,j)
        P(n-1,j)       = P(n+1,j)
        P(n  ,j)       = P(n+2,j)
        q(n-1,j)       = q(n+1,j)
        q(n  ,j)       = q(n+2,j)
        t(n-1)         = t(n+1)
        t(n  )         = t(n+2)
        E(n-1,j)       = E(n+1,j)
        E(n  ,j)       = E(n+2,j)
        K(n-1,j)       = K(n+1,j)
        K(n  ,j)       = K(n+2,j)
        Temp(n-1,j)    = Temp(n+1,j)
        Temp(n  ,j)    = Temp(n+2,j)
        entropy(n-1,j) = entropy(n+1,j)
        entropy(n  ,j) = entropy(n+2,j)
c        icompact(n-1,j)= icompact(n+1,j)
c        icompact(n  ,j)= icompact(n+2,j)

      endif
      enddo
cccccccccccccccccccccc  Main loop closed ccccccccccccccccccccccccccccccccccccc
       if (t(n) .le. tstop) goto 1
c      enddo  ! n or time loop
       close(33)
      print *,'*******************  Finished!!  ******************'
c      read(*,*) foo
       end
