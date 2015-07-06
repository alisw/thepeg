*CMZ :          11/10/99  20.16.42  by  Tancredi Carli
*CMZ :  1.01/01 28/03/96  15.50.59  by  Tancredi Carli
*-- Author :
      subroutine HZ98051(IFLAG)
****************************************************************************
* Purpose: Produce fraction of events with a forward jet using a cone
* algorithm (R=1).
*
* The event is accepted if it fulfills the requirements:
*   y>0.1, E_positron >10 GeV
*   0.00045 < x < 0.045
*   1 jet with cone algorithm (radius 1, PXCONE) and:
*     Et > 5 GeV
*     0 < eta < 2.6
*     x_jet > 0.036
*     0.5 < E_t**2/Q2 < 2
*     cos(theta_breit) > 0 (jet in fragmentation region of the breit frame)
*
*    Only 1 entry per event (even if there are 2 jets fulfilling the cuts).
*
*    histograms:
*    1    forward jet cross section in x_bjorken bins in the MC.
*   -1    forward jet cross section in x_bjorken bins in the data.
*     The data are corrected to hadron level. The parton level corrections
*     are not included.
*
*     700  eta distribution of the highest xjet jet
*     701  eta distribution of the second highest xjet jet
*     710  Et distribution of the highest xjet jet
*     711  Et distribution of the second highest xjet jet
*     720  Et2/Q2 distribution of the highest xjet jet (no Et2/Q2 cut applied)
*     -720 Et2/Q2 distribution of the highest xjet jet (no Et2/Q2 cut applied)
*          in the data
*     721             "            second      "
*
*     These are parton level results. A direct comparison with the data is
*     bound to hadronization uncertainties.
*
*     300  MEPJET cross sections (scale 0.25*kt**2) for forward jets
*     301  MEPJET cross sections (scale 2*kt**2) for forward jets
*     400  BFKL LO for forward jets
*     401  BFKL born for forward jets
*
*     The program gives the number of events with 1 or 2 forward jets
*     The systematics (in the data) are printed at the end
*
* Event selection:
*
* Running:
*
* Arguments: iflag=1 initialisation
*            iflag=1 filling
*            iflag=3 termination
*
* written by:
****************************************************************************
      IMPLICIT NONE
*KEEP,HZFUNC.
*
* Function declarations for Hztool functions
*
          DOUBLE PRECISION HzPhmang
          DOUBLE PRECISION HzDiskin
          DOUBLE PRECISION HzPhokin
          DOUBLE PRECISION HZETA
          DOUBLE PRECISION HZPHI
          DOUBLE PRECISION HZET
          DOUBLE PRECISION HZPT
          DOUBLE PRECISION HZTHETA
          DOUBLE PRECISION hzeekin
          Integer hzeebeam
          Integer hzeegamn
          Integer HzIpgamn
          Integer HzIdelec
          Integer HzIpgam
          Integer HzIbeam
          Integer HzLchge
          Integer HzLcomp
*
*KEEP,HEPEVTP.
*
* HEP event prime common
* (for explanation see manual)
      Integer NMXHEP
      PARAMETER (NMXHEP=4000)
      Integer NEVHEP,NHEP,ISTHEP,IDHEP
      Integer JMOHEP,JDAHEP
      Double Precision PHEP,VHEP
      COMMON/HEPEVTP/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     & JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
*
*KEEP,HERACMN.
*
* HERA common
*
*     GEN: Name of generator
*     XSEC: total cross section (in pb)
*     IHCHRG: charge of particle/parton times 3
*     NTOT : Number of total events
*     WTX  : event weight
*
      Character*8 Gen
      Double Precision Xsec
      Integer ihchrg
      Real    wtx, Ntot
      Common /HERACMN/ Xsec, Gen, ihchrg(nmxhep), Ntot,wtx
*
*KEND.
*
      Integer iprog
      Real qp_3,qp_2
      Real pthg1,shg1,xpg1,zqg1
      Common/partons/pthg1,shg1,xpg1,zqg1,iprog
c  ntuple jets
      Integer njet,ipro
      Real pthg,shg,xpg,zqg,pth,sh,xp,zq1,zq2,deta
      common/cwnco1/x,q2,y,ipro
      Integer kind(32)
      Integer numjet
      Parameter (numjet=50)
      Double precision pj(numjet,8)
      Double precision pbeam(4),pgam(4),pcm(4),ph(5)
      Double precision pold(4),pnew(4)
*
      Integer i,iflag,ihep,k1
      Character*5 xxxx
*
      Real pi,rd,eb,ee
      Parameter (ee=27.5)
      Integer iel,idum,ibeam,igam,ifj,afj,bfj
      Parameter (pi=3.1415927,rd=180./pi)
      Integer nentry
      Real x,y,q2,enel,thel,efwd,th,lx,lxj
      Real ptj21,pt2q21,etaj1
      Real thj,xj,ptb,costhb
      Integer ibin,nj,ij,naj,naj2,np,np4,np5,nsel
*
      Integer ierr
      Integer maxhi,modjet,maxjet
      Parameter (maxhi=100,modjet=7,maxjet=30)
      Real  selj(maxjet,4),help,empz,empzh,xg,empzl,empzhl,xgl
      Integer ipjet(maxjet),l,k,j
      Real nev(maxhi),lum,nev2(maxhi)
      save nev,nev2
      Double precision rcone
      Parameter (rcone=1.)
      Logical lp,fjet,twojet
      Data lp/.false./
      Real nxg(3)
*

      real q2as,ulalps,a1,a2
      double precision alphas2,alphas
      real xjet1
      real etaj(5),xjet(5),pt2q2(5),ptj2(5)
      real eta0,pt2q20,xjet0,pt0

*********************************************************************
*
      Integer nbin,nbin1
      Parameter (nbin=8)
      Parameter (nbin1=8)
      real xout(nbin),xerr(nbin)
      real ptqout(nbin1),ptqerr(nbin1)
      real fx(nbin+1),ptqbin(nbin1+1)
      Data   fx /0.00045,0.0008,0.0014,0.0025,0.0045,0.008,0.014,
     +            0.025,0.045/
      data ptqbin /0.01,0.03,0.1,0.3,1.,3.,10.,30.,100./
*
*********************************************************************
      Integer nx
      Real n2,nall,nall2
      Real xmin,xmax
      save xxxx,nall
      Data xxxx/'98051'/,NENTRY/0/
*********************************************************************
*
*                      data points from zeus
*
*********************************************************************
      real sup(8),sdown(8)
      real fjxs(8),efjxs(8),sysup(8),sysdown(8)
      real fjetq(8),efjetq(8)
      data fjxs /114.,96.2,77.8,34.4,14.1,6.53,2.65,0.65/
      data efjxs /9.7,6.5,4.7,2.2,1.0,0.54,0.25,0.09/
      data sysup /29.,8.2,5.2,3.8,2.5,0.1,0.3,0.1/
      data sysdown /15.,8.2,6.9,1.9,1.2,0.7,0.3,0.4/
      data fjetq /59.5,164,255,288,190,41.2,2.95,0.120/
      data efjetq /32.3,32.5,21.8,11.6,5.93,1.41,0.19,0.02/
      save sup,sdown
*********************************************************************
*
*                      MEPJET cross sections
*
*********************************************************************
*
* for scale 0.25*kt2
*
      real mep1(8)
      data mep1 /22.90,26.96,20.99,13.46,6.20,2.73,1.24,0.34/
*
* forscale 2*kt2
*
      real mep2(8)
      data mep2 /13.46,21.29,16.80,10.41,5.18,2.43,1.05,0.32/
*********************************************************************
*
*                      BARTELS and Wusthoff cross sections
*
*********************************************************************
*
*     LO BFKL
*
      real bfkl(8)
      data bfkl /147.7,150.9,104.6,46.51,16.26,5.25,1.44,0.23/
*
*     1st term
*
      real born(8)
      data born /28.93,35.50,34.26,18.32,8.38,3.20,1.07,0.19/
c sudakov check
      Integer ihist
      Integer ngl_tot,nsreg
      Real zgen,zminus,qp,qm,ktqt
      common /myzval/zgen(500),zminus(500),qp(500),qm(500),
     +   nsreg(500),ngl_tot,ktqt(500)
      real qplus,qminus
c end
*********************************************************************
*
*                      Initialization
*
*********************************************************************
      NENTRY=NENTRY+1
      wtx=1.0
*
      IF (iflag.eq.1) then
*
         nxg(1)=0.
         nxg(2)=0.
         nxg(3)=0.
         nall=0.
         nall2=0.
         n2=0.
         do i = 1,maxhi
            nev(i) = 0.
            nev2(i) = 0.
         enddo
*
*
* Initialisation: The following MUST always be done
* (i) make subdirectory in PAWC
* - use the name as the xxxxxx in HZxxxxxx subroutine
* (i) make subdirectory in o/p file
*
         Call hcdir('//PAWC',' ')
         call hmdir(xxxx,'S')
         Call hcdir('//HISTO',' ')
         call hmdir(xxxx,'S')
*
* book histos
*
*
         call hbook1(800,'  nglu    ',25,0.,25.,0.)
         call hbook1(900,'  z  ns=0  ',50,0.,1.,0.)
         call hbook1(910,'  z+ all    ',50,0.,1.,0.)
         call hbook1(901,'  z  ns=1  ',50,0.,1.,0.)
         call hbook1(902,'  z  ns=2  ',50,0.,1.,0.)
         call hbook1(903,'  z  ns=3  ',50,0.,1.,0.)
         call hbook1(921,'  z+; kt/qt.ge.1   ',50,0.,1.,0.)
         call hbook1(922,'  z+; z.lt.kt/qt.lt.1   ',50,0.,1.,0.)
         call hbook1(923,'  z+; kt/qt.le.z  ',50,0.,1.,0.)
         call hbook1(930,'  z- all    ',50,0.,1.,0.)
         call hbook1(931,'  z-; kt/qt.ge.1   ',50,0.,1.,0.)
         call hbook1(932,'  z-; z.lt.kt/qt.lt.1   ',50,0.,1.,0.)
         call hbook1(933,'  z-; kt/qt.le.z  ',50,0.,1.,0.)
         call hbook2(940,'  z+ vrs z- all',50,0.,1.,50,0.,1.,0.)
         call hbook2(941,'  z+ vrs z-; kt/qt.ge.1',50,0.,1.,50,0.,1.,0.)
         call hbook2(942,' z+ vrs z-; z.lt.kt/qt.lt.1',50,0.,1.,
     +    50,0.,1.,0.)
         call hbook2(943,'  z+ vrs z-; kt/qt.le.z',50,0.,1.,50,0.,1.,0.)

         call hbook1(951,'  z+; kt**2.gt.z*qt**2   ',50,0.,1.,0.)
         call hbook1(952,'  z+; kt**2.lt.z*qt**2   ',50,0.,1.,0.)
         call hbook1(953,'  z-; kt**2.gt.z*qt**2   ',50,0.,1.,0.)
         call hbook1(954,'  z-; kt**2.lt.z*qt**2   ',50,0.,1.,0.)

         call hbook1(961,'  kt/qt for q+ ordered',50,0.,10.,0.)
         call hbook1(962,'  kt/qt for q+ non-ordered',50,0.,10.,0.)
         call hbook1(963,'  kt/qt for q- ordered',50,0.,10.,0.)
         call hbook1(964,'  kt/qt for q- non-ordered',50,0.,10.,0.)
         call hbook1(971,'  z+ for q+ ordered',50,0.,1.,0.)
         call hbook1(972,'  z+ for q+ non-ordered',50,0.,1.,0.)
         call hbook1(973,'  z+ for q- ordered',50,0.,1.,0.)
         call hbook1(974,'  z+ for q- non-ordered',50,0.,1.,0.)

         call hbook1(31,' nall/cross ',3,0.,3.,0.)
         call hbookb(1 ,' 1 fjets XBj ',nbin,fx,0.)
         call hbookb(-1,' 1 fjets XBj ',nbin,fx,0.)
         call hbook1(700,'ethad 1st fjet',30,0.,30.,0.)
         call hbook1(701,'ethad 2nd fjet',30,0.,30.,0.)
         call hbook1(710,'etahad 1st fjet',15,0.,3.,0.)
         call hbook1(711,'etahad 2nd fjet',15,0.,3.,0.)
         call hbookb(720,'pt2q2 1st fjet',nbin1,ptqbin,0.)
         call hbookb(-720 ,' 1 fjets Et2/Q2 ',nbin1,ptqbin,0.)
         call hbookb(721,'pt2q2 2nd fjet',nbin1,ptqbin,0.)
         call hbookb(300 ,'MEPJET scale 0.25kt2 ',nbin,fx,0.)
         call hbookb(301 ,'MEPJET scale 2kt2',nbin,fx,0.)
         call hbookb(400 ,'BFKL LO ',nbin,fx,0.)
         call hbookb(401 ,'BFKL 1st term',nbin,fx,0.)
*
*      pack data into histograms ; 8 bins
*
*
         Call hpak (-1,fjxs)
         Call hpake(-1,efjxs)
         Do I=1,8
            fjetq(i)=fjetq(i)/1000.
            efjetq(i)=efjetq(i)/1000.
         Enddo
         Call hpak (-720,fjetq)
         Call hpake(-720,efjetq)
*
*     pack theory points
*
         call hpak(300,mep1)
         call hpak(301,mep2)
         call hpak(400,bfkl)
         call hpak(401,born)
*
*      add sytematical and statistical error
*
         Do i=1,8
            sup(i) = sqrt(efjxs(i)**2+sysup(i)**2)
            sdown(i) = sqrt(efjxs(i)**2+sysdown(i)**2)
         enddo
*
*********************************************************************
*
*                      Event Processing
*
*********************************************************************
      Else if(iflag.eq.2) then
         do i = 1,5
            etaj(i) = 0.
            pt2q2(i) = 0.
            xjet(i) = 0.
            ptj2(i) = 0.
         enddo

*
*     Filling: The following MUST always be done
*     (i) move to the correct sub-directory in PAWC
*
         call hcdir('//PAWC/'//xxxx,' ')
*
         nall=nall+wtx
*
         ierr=HZIBEAM(ibeam,idum)
         if (.not.(ierr.eq.1)) then
            write(6,*) 'HZ'//xxxx,' beams not found ! '
            return
         else
            Do i=1,4
               pbeam(i)=PHEP(i,IBEAM)
            enddo
         endif
*
         eb=real(PHEP(4,idum))
         if (abs(abs(eb)-ee).gt.0.1) then
            if (nentry.lt.10) then
               write(6,*) 'HZ'//xxxx,' Electron not at ',ee,' ! ',eb
            elseif (nentry.eq.10) then
               write(6,*) 'HZ'//xxxx,' Electron beam not at ',ee,
     +              ' ...last message ! '
            endif
c            return
         endif
*
* Event selection
*
         q2=real(HZDISKIN(1))
         x =real(HZDISKIN(2))
         y =real(HZDISKIN(3))

         deta = -999.
         njet = 0
         ipro = iprog
*
*
*     get electron variables for kinematic selection
*
         iel=HZIDELEC(idum)
         if (iel.eq.-1) then
            write(*,*) 'Hz'//xxxx,' electron not found '
            return
         endif
         enel=real(PHEP(4,iel))
         thel=real(HZPHMANG(PHEP(3,iel), sqrt(PHEP(1,iel)**2+PHEP(2,
     +   iel)**2)))
*
*     Momentum of the gamma
*
         pgam(1)=-real(PHEP(1,iel))
         pgam(2)=-real(PHEP(2,iel))
         pgam(3)=-27.5 - real(PHEP(3,iel))
         pgam(4)= 27.5 - real(PHEP(4,iel))
*
         if (lp) then
            write(6,*) 'Hz'//xxxx,' enel= ',enel,' thel= ',thel
            write(6,*) 'Hz'//xxxx,' y= ',y,' q2= ',q2
         endif
*
         fjet=.true.
         if (y.lt.0.1) fjet=.false.
         if (enel.lt.10.0) then
            fjet=.false.
         endif
*
* Initialize Breit frame boost
*
         call HZBRTINI(pbeam,pgam,ierr)
         if (ierr.eq.1) then
           write(*,*) 'Hz'//xxxx,' Boost failed'
         endif
*
* find jets
*
         if (fjet) then
*
            call hzjtfind(modjet,rcone,nj,pj)
*
            if (lp) then
               write(6,*) 'forward jet analysis: nj= ',nj
            endif
            if (nj.gt.numjet) then
               write(6,*) 'hz',xxxx,' too many jets found ! '
               return
            endif
*
            naj=0
            ifj=0
            do 100 ij=1,nj
*
               pold(1) = pj(ij,5)
               pold(2) = pj(ij,6)
               pold(3) = pj(ij,7)
               pold(4) = pj(ij,4)
*
*     Boost jet in the breit frame
*
               call HZBRT(pold,pnew,ierr)
               if (ierr.eq.1) then
                 write(*,*) 'Hz'//xxxx,' Boost failed'
               endif
               ptb = sqrt(real(pnew(1))**2+real(pnew(2))**2)
               costhb = real(pnew(3)) / sqrt(ptb**2 + real(pnew(3))**2)
*
               if (lp) write(6,'(a,4f9.3)') ' Jet= ',(real(pj(ij,i)),i=
     +         4,7)
*
               ptj21=real(pj(ij,3))
               pt2q21=ptj21**2/q2
               thj=real(HZPHMANG(pj(ij,7),pj(ij,3)) )
               etaj1 = pj(ij,1)
*
*     Get the forward jets
*
               if ( ptj21.gt.5.and. etaj1.lt.2.6.and. etaj1.gt.0.)
     +         then
                  xjet1=real(pj(ij,4))/820.*cos(thj)
                  if (costhb.gt.0) then
                     ifj = ifj + 1
                     xjet(ifj) = xjet1
                     ptj2(ifj) = ptj21
                     pt2q2(ifj) = pt2q21
                     etaj(ifj) = etaj1
                  endif
               endif
  100       continue
            if (ifj.eq.0) return
*
*     Save an event if there is a f.j. with 0.5<Et2/Q2<2
*     and xjet>0.036
*
            afj = 0
            bfj = 0
            do i = 1,ifj
               if (xjet(i).gt.0.036) then
                  bfj = bfj + 1
               endif
               if (pt2q2(i).lt.2.and.pt2q2(i).gt.0.5) then
                  if (xjet(i).gt.0.036) then
                     afj = afj + 1
                  endif
               endif
            enddo
            if (afj.ge.1) then
               call hf1(1,x,wtx)
c store z values from SMALLX
               if(gen(1:3).eq.'SMA') then
                  qplus = qp(1)
                  qminus = qm(1)
                  call hfill(800,real(ngl_tot+1.3),0.,wtx)
                  do i=1,ngl_tot
c			write(6,*) ' in loop ',i
                     call hfill(910,zgen(i),0.,wtx)
                     call hfill(930,zminus(i),0.,wtx)
                     call hfill(940,zgen(i),zminus(i),wtx)
                     if(nsreg(i).eq.0) then
                        call hfill(900,zgen(i),0.,wtx)
                     elseif(nsreg(i).eq.1) then
                        call hfill(901,zgen(i),0.,wtx)
                     elseif(nsreg(i).eq.2) then
                        call hfill(902,zgen(i),0.,wtx)
                     elseif(nsreg(i).eq.3) then
                        call hfill(903,zgen(i),0.,wtx)
                     endif
                     if(ktqt(i).ge.1.) then
                        call hfill(921,zgen(i),0.,wtx)
                        call hfill(931,zminus(i),0.,wtx)
                        call hfill(941,zgen(i),zminus(i),wtx)
                     elseif(ktqt(i).lt.1..and.ktqt(i).gt.zgen(i)) then
                        call hfill(922,zgen(i),0.,wtx)
                        call hfill(932,zminus(i),0.,wtx)
                        call hfill(942,zgen(i),zminus(i),wtx)
                     elseif(ktqt(i).le.zgen(i)) then
                        call hfill(923,zgen(i),0.,wtx)
                        call hfill(933,zminus(i),0.,wtx)
                        call hfill(943,zgen(i),zminus(i),wtx)
                     endif

c			    write(6,*) ' check ktqt ',ktqt(i),' z ',zgen(i),i
                     if(ktqt(i)**2.ge.zgen(i)) then
c			      write(6,*) ' ktqt > z ',i
                        call hfill(951,zgen(i),0.,wtx)
                        call hfill(953,zminus(i),0.,wtx)
                     else
c			      write(6,*) ' ktqt < z ',i
                        call hfill(952,zgen(i),0.,wtx)
                        call hfill(954,zminus(i),0.,wtx)
                     endif
c        write(6,*) ' hzgluon qplus qminus ',qp(i),qm(i)
                     if(qp(i).le.qplus) then
                        call hfill(961,real(ktqt(i)),0.,wtx)
                        call hfill(971,zgen(i),0.,wtx)
                     else
			     write(6,*) ' non ordering ',qp(i),qplus
                        call hfill(962,real(ktqt(i)),0.,wtx)
                        call hfill(972,zgen(i),0.,wtx)
                     endif
                     qplus=qp(i)
c			   write(6,*) ' check for non ordering ',qm(i),qminus
                     if(qm(i).ge.qminus) then
                        call hfill(963,real(ktqt(i)),0.,wtx)
                        call hfill(973,zgen(i),0.,wtx)
                     else
			     write(6,*) ' non ordering ',qm(i),qminus
                        call hfill(964,real(ktqt(i)),0.,wtx)
                        call hfill(974,zgen(i),0.,wtx)
                     endif
                     qminus=qm(i)
                  enddo
               endif


            endif
*
            if (x.gt.fx(1).and.x.lt.fx(9)) then
               if (afj.eq.1) nev(1)=nev(1)+wtx
               if (afj.eq.2) nev2(1)=nev2(1)+wtx
            endif
*
*     Sort jets by xjet
*
            if (ifj.eq.1) goto 12
            xjet0 = 0.
            do j = 1,ifj
               xjet0 = xjet(j)
               l = j+1
               k1 = 0
               do k = l,ifj
                  if (xjet(k).gt.xjet0) then
                     xjet0 = xjet(k)
                     eta0 = etaj(k)
                     pt0 = ptj2(k)
                     pt2q20 = pt2q2(k)
                     k1 = k
                  endif
               enddo
               if (k1.eq.0) goto 12
               xjet(k1) = xjet(j)
               etaj(k1) = etaj(j)
               ptj2(k1) = ptj2(j)
               pt2q2(k1) = pt2q2(j)
               xjet(j) = xjet0
               etaj(j) = eta0
               ptj2(j) = pt0
               pt2q2(j) = pt2q20
            enddo
   12       continue
*
*     Store the histograms
*     highest xjet
*
            if (xjet(1).gt.0.036) then
****           if (x.gt.0.00025.and.x.lt.0.08) then
               if
     +         (x.gt.0.00025.and.x.lt.0.08.and.q2.gt.10.and.bfj.ge.1)
     +         then
                  call hf1(720,pt2q2(1),wtx)
                  if (pt2q2(1).gt.0.5.and.pt2q2(1).lt.2.) then
                     call hf1(700,ptj2(1),wtx)
                     call hf1(710,etaj(1),wtx)


                  endif
               endif
            endif
*
*     second highest xjet
*
            if (ifj.eq.1) goto 2000
            if (xjet(2).gt.0.036) then
               if (x.gt.fx(1).and.x.lt.fx(9)) then
                  call hf1(721,pt2q2(2),wtx)
                  if (pt2q2(2).gt.0.5.and.pt2q2(2).lt.2.) then
                     call hf1(701,ptj2(2),wtx)
                     call hf1(711,etaj(2),wtx)
                  endif
               endif
            endif
*
 2000       continue
         endif
*
*********************************************************************
*
*     Termination
*
*********************************************************************

      Else if(iflag.eq.3) then

* Termination: The following MUST always be done
* (i) Move to the correct PAW subdirectory
*
         call hcdir('//PAWC/'//xxxx,' ')
*
         lum=999999999.
         if (xsec.ne.0.) then
            lum=real(nall)/real(xsec)
         else
            write(6,*) 'hz',xxxx,' xsec=0 ! '
         endif
*
         write(6,*) 'hz98051: Total Cross= ', real(xsec),' nall= ',
     +   nall
         write(6,*) '  ZEUS: # of 1 fjets= ',nev(1),' lumi= ',lum
         write(6,*) '  ZEUS: # of 2 fjets= ',nev2(1),' lumi= ',lum
         write(6,*) 'Systematics up (+ statistics)', sup
         write(6,*) 'Systematics down (+ statistics)', sdown
*
* convert in nb
*
         lum = lum * 1000.
*
*
*
         call HBOOKB(1111,'Tx',8,fx,0.)
         call HBOOKB(1112,'Tetq',8,ptqbin,0.)
         call HBARX(0)

         do I = 1,8
            xout(i)=fx(i+1)-fx(i)
            xerr(i)=0.
         enddo
         call HPAK(1111,xout)
         call HPAKE(1111,xerr)
*
         do I = 1,8
            ptqout(i)=ptqbin(i+1)-ptqbin(i)
            ptqerr(i)=0.
         enddo
         call HPAK(1112,ptqout)
         call HPAKE(1112,ptqerr)
*
*     Normalise to luminosity
*
         call hopera(1,'+E',1,1,1./lum,0.)
         call hopera(700,'+E',700,700,1./lum,0.)
         call hopera(710,'+E',710,710,1./lum,0.)
         call hopera(720,'+E',720,720,1./lum,0.)
         call hopera(701,'+E',701,701,1./lum,0.)
         call hopera(711,'+E',711,711,1./lum,0.)
         call hopera(721,'+E',721,721,1./lum,0.)
         call hopera(1,'/E',1111,1,1.,1.)
         call hopera(720,'/E',1112,720,1.,1.)
         call hzhinfo(1,int(nev(1)))
*
         do i=0,3
            ihist = 900+i
            call hopera(ihist,'+E',ihist,ihist,1./lum,0.)
         enddo
         do i=1,3
            ihist = 920+i
            call hopera(ihist,'+E',ihist,ihist,1./lum,0.)
         enddo
         call hopera(910,'+E',910,910,1./lum,0.)
         do i=1,3
            ihist = 930+i
            call hopera(ihist,'+E',ihist,ihist,1./lum,0.)
         enddo
         call hopera(930,'+E',930,930,1./lum,0.)
      endif
*
      RETURN
      END
