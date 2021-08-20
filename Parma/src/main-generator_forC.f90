! Generate cosmic-ray based on PARMA model
! Forked and modified from 
! https://gitlab.cern.ch/mathusla/CosmicSpectra/-/blob/7b8d8ca669f2e29456f51898e19904694b3fe333/main-generator.f90
! Input order: /nebin, nabin, ip, seed/, /ie, ia/, mass, ehigh, etable, ahigh, atable (dims as below)
subroutine initializegen(iparams, ivars, mass, ehigh, etable, ahigh, atable, flux)
      
      parameter(npart=33) ! number of applicable particle
      implicit real*8 (a-h, o-z)
      dimension IangPart(0:npart)
      ! dimension etable(0:nebin)            ! probability table (0.0 for 0, 1.0 for nebin)
      ! dimension atable(0:nabin,0:nebin)    ! probability table (0.0 for 0, 1.0 for nabin) 
      ! dimension flux(0:nabin,nebin) ! Monte Carlo generated flux
      data IangPart/1,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,5,5,6/ ! Angular particle ID
      character, parameter :: tab = char(9)
      
      integer, intent(in), dimension(0:3) :: iparams
      integer, intent(out), dimension(0:1) :: ivars
      real*8, intent(out) :: mass
      real*8, intent(out), dimension(0:iparams(0)) :: ehigh
      real*8, intent(out), dimension(0:iparams(0)) :: etable
      real*8, intent(out), dimension(0:iparams(1)) :: ahigh
      real*8, intent(out), dimension(0:iparams(1),0:iparams(0)) :: atable
      real*8, intent(out), dimension(0:iparams(1),iparams(0)) :: flux
      dimension emid(iparams(0)) ! higher and middle point of energy bin
      dimension amid(iparams(1)) ! higher and middle point of angular bin

      ! parameter(initial_seed=iparams(3)) ! Seed for RNG
      initial_seed = iparams(3)
      ! parameter(nebin=iparams(0)) ! number of energy mesh (divided by log)
      nebin = iparams(0)
      ! parameter(nabin=iparams(1)) ! number of angle mesh (divided by linear)
      nabin = iparams(1)

      if (initial_seed == 0) then
      write(*,*) 'Please set the random seed'
      stop
      end if

      ! Set number of particle to be generated, and initial random seed
      ! nevent=10000 ! number of particles to be generated
      call sgrnd(initial_seed) ! initial seed of random number (you can change the value)

      ! Set Particle ID
      ip = iparams(2) ! Particle ID (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)
      if(IangPart(ip).eq.0) then
      write(*,*) 'Angular distribution is not available for the particle'
      stop
      endif

      ! Set Conditions (location, time, and local geometry)

      ! Stable P1 data taking
      iyear=2019  ! Year
      imonth=11    ! Month
      iday=1      ! Date

      ! SX1 buffer zone location
      glat=46.23572   ! Latitude (deg), -90 =< glat =< 90
      glong=6.05518 ! Longitude (deg), -180 =< glong =< 180
      Alti=0.439   ! Altitude (km)

      g=100.0      ! Local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin

      s=getHP(iyear,imonth,iday,idummy) ! Solar activity (W index) of the day
      r=getr(glat,glong)                ! Vertical cut-off rigidity (GV)
      d=getd(alti,glat)                 ! Atmospheric depth (g/cm2), set glat = 100 for use US Standard Atmosphere 1976.

      ! Set energy and angle ranges for generation
      if (ip == 0) then
      emin = 1.0e-8 ! Minimum energy of particle (MeV)
      elseif (ip == 33) then
      emin = 10.0 ! Minimum energy of particle (MeV)
      else
      emin = 1.0e3 ! Minimum energy of particle (MeV)
      end if
      if (ip == 29 .or. ip == 30) then
      emax = 1.0e8 ! Maximum energy of particle (MeV)
      else
      emax = 1.0e6 ! Maximum energy of particle (MeV)
      end if
      amin = 0.0   ! Minimum cosine of particle
      amax = 1.0   ! Maximum cosine of particle

      ! Make energy and angle mesh
      elog=log10(emin)
      estep=(log10(emax)-log10(emin))/nebin
      do ie=0,nebin
      ehigh(ie)=10d0**elog
      if(ie.ne.0) emid(ie)=sqrt(ehigh(ie)*ehigh(ie-1))
      elog=elog+estep
      enddo

      astep=(amax-amin)/nabin
      do ia=0,nabin
      ahigh(ia)=amin+astep*ia
      if(ia.ne.0) amid(ia)=(ahigh(ia)+ahigh(ia-1))*0.5
      enddo

      ! Make probability table (absolute value)
      atable(:,:)=0.0d0 ! initialization
      etable(:)=0.0d0
      do ie=1,nebin
      do ia=1,nabin
      atable(ia,ie)=atable(ia-1,ie)+getSpec(ip,s,r,d,emid(ie),g) &
      & *getSpecAngFinal(iangpart(ip),s,r,d,emid(ie),g,amid(ia))*(2.0*acos(-1.0))*(ahigh(ia)-ahigh(ia-1)) ! angular integrated value
      enddo
      enddo
      do ie=1,nebin
      etable(ie)=etable(ie-1)+atable(nabin,ie)*(ehigh(ie)-ehigh(ie-1)) ! energy integrated value
      enddo
      TotalFlux=etable(nebin) ! Total Flux (/cm2/s), used for normalization

      ! Make probability table (normalized to 1)
      do ie=1,nebin
      etable(ie)=etable(ie)/etable(nebin)
      do ia=1,nabin
      atable(ia,ie)=atable(ia,ie)/atable(nabin,ie)
      enddo
      enddo

      write(*,'(a,i2,a,f8.6,a)') '# Particle ', ip, ' total flux = ', TotalFlux * 60.0, ' cm^-2 min^-1'
      write(*,*)
      ! write(*,'(a,i10)') '# Seed = ', initial_seed
      ! write(*,*)
      ! write(*,'(a)') '# Format (by column):'
      ! write(*,*)
      ! write(*,'(a)') '# Particle energy (MeV)'
      ! write(*,'(a)') '# Zenith angle (rad)'
      ! write(*,'(a)') '# Azimuthal angle (rad)'
      ! write(*,*)

      if (ip == 0) then
      mass = 939.565413d0
      else if (ip == 1) then
      mass = 938.272081d0
      else if (ip == 29 .or. ip == 30) then
      mass = 105.6583745d0
      else if (ip == 31 .or. ip == 32) then
      mass = 0.5109989461d0
      else if (ip == 33) then
      mass = 0.0d0
      else
      write(*,*) 'Unknown particle number', ip
      stop
      end if

      ivars(0) = ie
      ivars(1) = ia

end subroutine

! Particle Generation
subroutine getparticle(iparams, ivars, mass, ehigh, etable, ahigh, atable, flux, particle)

      implicit real*8 (a-h, o-z)
      
      integer, intent(in), dimension(0:3) :: iparams
      integer, intent(in), dimension(0:1) :: ivars

      real*8, intent(in) :: mass
      real*8, intent(in), dimension(0:iparams(0)) :: ehigh
      real*8, intent(in), dimension(0:iparams(0)) :: etable
      real*8, intent(in), dimension(0:iparams(1)) :: ahigh
      real*8, intent(in), dimension(0:iparams(1),0:iparams(0)) :: atable
      real*8, intent(out), dimension(0:iparams(1),iparams(0)) :: flux
      real*8, intent(out), dimension(0:3) :: particle

      nebin = iparams(0)
      nabin = iparams(1)
      ip = iparams(2)
      ie = ivars(0)
      ia = ivars(1)

      e=getGeneration(ie,nebin,ehigh(0),etable(0))    ! energy
      phi=2.0*acos(-1.0)*(grnd()-0.5) ! azimuth angle (rad)
      w=getGeneration(ia,nabin,ahigh(0),atable(0,ie)) ! z direction, -1.0:upward, 0.0:horizontal, 1.0:downward
      ! write(*,'(f24.9,a,f8.6,a,f16.6)') e + mass, tab, dacos(w), tab, phi

      particle(0) = ip
      particle(1) = e + mass
      particle(2) = dacos(w)
      particle(3) = phi

      flux(ia,ie)=flux(ia,ie)+TotalFlux/(ehigh(ie)-ehigh(ie-1))/((ahigh(ia)-ahigh(ia-1))*2.0*acos(-1.0)) ! /cm2/s/sr/MeV
      flux(0,ie)=flux(0,ie)+TotalFlux/(ehigh(ie)-ehigh(ie-1)) ! /cm2/s/MeV

end subroutine

!open(2,file='CheckGeneration.out')
!write(2,*) 'Total Flux (/cm2/s) =',TotalFlux
!write(2,'(102i12)') (ia,ia=1,nabin+2)
!write(2,'(''   Emid(MeV)  /cm2/s/MeV'',100f12.4)') (amid(ia),ia=1,nabin)
! do ie=1,nebin
!  write(2,'(102es12.4)') emid(ie),(flux(ia,ie)/nevent,ia=0,nabin)
! enddo

end

function getGeneration(ibin,nbin,high,table)
implicit real*8 (a-h, o-z)
dimension high(0:nbin)
dimension table(0:nbin)

rand=grnd() ! random number
do i=1,nbin-1
 if(rand.le.table(i)) exit
enddo
ibin=i ! bin ID

rand=grnd() ! random number
getGeneration=high(ibin-1)*rand + high(ibin)*(1.0d0-rand)

return

end


! A C-program for MT19937: Real number version
! http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/FORTRAN/fortran.html
!   genrand() generates one pseudorandom real number (double)
! which is uniformly distributed on [0,1]-interval, for each
! call. sgenrand(seed) set initial values to the working area
! of 624 words. Before genrand(), sgenrand(seed) must be
! called once. (seed is any 32-bit integer except for 0).
! Integer generator is obtained by modifying two lines.
!   Coded by Takuji Nishimura, considering the suggestions by
! Topher Cooper and Marc Rieffel in July-Aug. 1997.
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later
! version.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General
! Public License along with this library; if not, write to the
! Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
! 02111-1307  USA
!
! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.
!
!***********************************************************************
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!
!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed
!
! This program uses the following non-standard intrinsics.
!   ishft(i,n): If n>0, shifts bits in i by n positions to left.
!               If n<0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
!
!***********************************************************************
      subroutine sgrnd(seed)
!
      implicit integer(a-z)
!
! Period parameters
      parameter(N     =  624)
!
      dimension mt(0:N-1)
!                     the array for the state vector
      common /block/mti,mt
      save   /block/
!
!      setting initial seeds to mt[N] using
!      the generator Line 25 of Table 1 in
!      [KNUTH 1981, The Art of Computer Programming
!         Vol. 2 (2nd Ed.), pp102]
!
      mt(0)= iand(seed,-1)
      do 1000 mti=1,N-1
        mt(mti) = iand(69069 * mt(mti-1),-1)
 1000 continue
!
      return
      end
!***********************************************************************
      double precision function grnd()
!
      implicit integer(a-z)
!
! Period parameters
      parameter(N     =  624)
      parameter(N1    =  N+1)
      parameter(M     =  397)
      parameter(MATA  = -1727483681)
!                                    constant vector a
      parameter(UMASK = -2147483648)
!                                    most significant w-r bits
      parameter(LMASK =  2147483647)
!                                    least significant r bits
! Tempering parameters
      parameter(TMASKB= -1658038656)
      parameter(TMASKC= -272236544)
!
      dimension mt(0:N-1)
!                     the array for the state vector
      common /block/mti,mt
      save   /block/
      data   mti/N1/
!                     mti==N+1 means mt[N] is not initialized
!
      dimension mag01(0:1)
      data mag01/0, MATA/
      save mag01
!                        mag01(x) = x * MATA for x=0,1
!
      TSHFTU(y)=ishft(y,-11)
      TSHFTS(y)=ishft(y,7)
      TSHFTT(y)=ishft(y,15)
      TSHFTL(y)=ishft(y,-18)
!
      if(mti.ge.N) then
!                       generate N words at one time
        if(mti.eq.N+1) then
!                            if sgrnd() has not been called,
          call sgrnd(4357)
!                              a default initial seed is used
        endif
!
        do 1000 kk=0,N-M-1
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
 1000   continue
        do 1100 kk=N-M,N-2
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
 1100   continue
        y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
        mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti = 0
      endif
!
      y=mt(mti)
      mti=mti+1
      y=ieor(y,TSHFTU(y))
      y=ieor(y,iand(TSHFTS(y),TMASKB))
      y=ieor(y,iand(TSHFTT(y),TMASKC))
      y=ieor(y,TSHFTL(y))
!
      if(y.lt.0) then
        grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
      else
        grnd=dble(y)/(2.0d0**32-1.0d0)
      endif
!
      return
      end
