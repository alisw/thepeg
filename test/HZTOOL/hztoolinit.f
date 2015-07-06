      subroutine hztoolinit(filename)

      PARAMETER(NWPAWC=1000000)
      PARAMETER (NMXHEP=4000)

      COMMON /pawc/h(NWPAWC)
      COMMON /HERACMN/ Xsec, Gen, ichrg(nmxhep), Ntot, wtx, beams
      CHARACTER*8 beams
      DOUBLE PRECISION Xsec
      CHARACTER*8 Gen
      REAL wtx, Ntot
      INTEGER ichrg
      CHARACTER*(*) filename

      COMMON /HEPEV2/nmulti,jmulti(4000)
      COMMON /HEPEV3/nevmulti(16),itrkmulti(16),mltstr(16)
      COMMON /HEPEV4/eventweightlh,alphaqedlh,alphaqcdlh,scalelh(10),
     $     spinlh(3,4000),icolorflowlh(2,4000),idruplh
      COMMON /HEPEV5/eventweightmulti(16),alphaqedmulti(16),
     $     alphaqcdmulti(16),scalemulti(10,16),idrupmulti(16)

      DOUBLE PRECISION eventweightlh,alphaqedlh,alphaqcdlh,scalelh,
     $     spinlh,eventweightmulti,alphaqedmulti,alphaqcdmulti,
     $     scalemulti

      Integer NEVHEP,NHEP,ISTHEP,IDHEP
      Integer JMOHEP,JDAHEP
      Double Precision PHEP,VHEP
      COMMON/HEPEVTP/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     & JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)

      GEN='HRW'
      CALL hlimit(NWPAWC)

      CALL hropen(45,'HISTO',filename,'N', 1024,istat)

      return
      end

      subroutine hztoolredirect(filename)

      CHARACTER*(*) filename

      close(6)
      open(6,file=filename,status='UNKNOWN')

      return
      end

      subroutine hztoolfinish

      CALL hcdir('//PAWC',' ')
      CALL hcdir('//HISTO',' ')
      CALL hrout(0,icycle,'T')
      CALL hrend('HISTO')

      return
      end

      subroutine hztoolsetxsec(x, sumw)

      PARAMETER (NMXHEP=4000)
      COMMON /HERACMN/ Xsec, Gen, ichrg(nmxhep), Ntot, wtx, beams
      DOUBLE PRECISION Xsec
      CHARACTER*8 Gen
      CHARACTER*8 beams
      REAL wtx, Ntot
      INTEGER ichrg
      DOUBLE PRECISION x, sumw

      xsec=x
      ntot=sumw

      return
      end

      subroutine hztoolfixdaughters

      PARAMETER (NMXHEP=4000)

      Integer NEVHEP,NHEP,ISTHEP,IDHEP
      Integer JMOHEP,JDAHEP
      Double Precision PHEP,VHEP
      COMMON/HEPEVTP/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     & JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)

      do 100 i=1,NHEP
        if ( JMOHEP(1,i).le.0 ) goto 100
        j1=JMOHEP(1,i)
        j2=max(j1,JMOHEP(2,i))
        do 110 j=j1,j2
          if ( JDAHEP(1,j).le.0 ) then
            JDAHEP(1,j)=i
          else
            JDAHEP(1,j)=min(JDAHEP(1,j),i)
          endif
          JDAHEP(2,j)=max(JDAHEP(2,j),i)
 110    continue
 100  continue

      return
      end

      subroutine hztoolfixboson(px,py,pz,pe,pm)

      PARAMETER (NMXHEP=4000)

      Integer NEVHEP,NHEP,ISTHEP,IDHEP
      Integer JMOHEP,JDAHEP
      Double Precision PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     & JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)

      double precision px,py,pz,pe,pm

      NHEP=NHEP+1
      ISTHEP(NHEP)=3
      IDHEP(NHEP)=23
      PHEP(1,NHEP)=px
      PHEP(2,NHEP)=py
      PHEP(3,NHEP)=pz
      PHEP(4,NHEP)=pe
      PHEP(5,NHEP)=pm

      return
      end

      subroutine hztoolfixbeams

      PARAMETER (NMXHEP=4000)

      Integer NEVHEP,NHEP,ISTHEP,IDHEP
      Integer JMOHEP,JDAHEP
      Double Precision PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     & JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)

      do 100 i=1,NHEP
        if ( ISTHEP(i).eq.2.and.JMOHEP(1,i).eq.0 ) ISTHEP(i)=101
 100  continue

      return
      end

      subroutine hztoolfixmirdir

      PARAMETER (NMXHEP=4000)

      Integer NEVHEP,NHEP,ISTHEP,IDHEP
      Integer JMOHEP,JDAHEP
      Double Precision PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     & JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)

      do 100 i=1,NHEP
        PHEP(3,i)=-PHEP(3,i)
 100  continue

      return
      end


      integer function printhepevt()
      PARAMETER (NMXHEP=4000)

      Integer NEVHEP,NHEP,ISTHEP,IDHEP
      Integer JMOHEP,JDAHEP
      Double Precision PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     & JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)


      printhepevt=nhep

      do 100 i=1,NHEP
        write(0,600) i,isthep(i),idhep(i),jmohep(1,i),jmohep(2,i),
     $       jdahep(1,i),jdahep(2,i),(phep(j,i),j=1,5)
 100  continue

 600  format(i3,i4,i8,4i4,5f9.2)

      return
      end

      integer function printhepevtp()
      PARAMETER (NMXHEP=4000)

      Integer NEVHEP,NHEP,ISTHEP,IDHEP
      Integer JMOHEP,JDAHEP
      Double Precision PHEP,VHEP
      COMMON/HEPEVTP/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     & JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)


      printhepevtp=nhep

      do 100 i=1,NHEP
        write(0,600) i,isthep(i),idhep(i),jmohep(1,i),jmohep(2,i),
     $       jdahep(1,i),jdahep(2,i),(phep(j,i),j=1,5)
 100  continue

 600  format(i3,i4,i8,4i4,5f9.2)

      return
      end


      integer function printhepevtf()
      PARAMETER (NMXHEP=4000)

      Integer NEVHEP,NHEP,ISTHEP,IDHEP
      Integer JMOHEP,JDAHEP
      Double Precision PHEP,VHEP,psum(4)
      COMMON/HEPEVTP/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     & JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)


      printhepevtf=0
      do 10 j=1,4
        psum(j)=0.0d0
 10   continue

      do 100 i=1,NHEP
        if ( isthep(i).ne.1 ) goto 100
        printhepevtf=printhepevtf+1
        do 110 j=1,4
          psum(j)=psum(j)+phep(j,i)
 110    continue
        write(0,600) i,isthep(i),idhep(i),jmohep(1,i),jmohep(2,i),
     $       jdahep(1,i),jdahep(2,i),(phep(j,i),j=1,5)
 100  continue
      write(0,*)
      write(0,601) (psum(j),j=1,4),
     $     sqrt(max(psum(4)**2-psum(3)**2-psum(2)**2-psum(1)**2,0.0d0))

 600  format(i3,i4,i8,4i4,5f9.2)
 601  format(31x,5f9.2)

      return
      end


