
      program makeinput                !Main program start
	implicit none                  !define implicit variables
	integer::i,j,k,iso,a,ocount=0    !i,j,k,iso are integer
	character(4),dimension(10000)::oyear,omonth*2,oday*2, 
     & otime*6
	real,dimension(10000)::pwc,ozone,soil,water,veg,utc
      real,dimension(10000,4)::opt,ssa,asy
      character(80)::header
      real::snow=0.0
	character(len=2)::months(6)=(/'03','04','05','06','07','08'/)
	

c     ===========Instrument type==============
c     Terms for INPUTFILE(SBDART model)
      integer::idatm  !Mid-Lat Winter for Atmospheric profile Mar-Aug
	real(8)::wlinf=0.20,wlsup=4.0,alat=43.577,alon=104.419
	integer::iaer=5  !aerosol=5,aerosol-free=0
      
c     ========================================
	
c         Read observation data
c         --------------------------
        !do i=1,6
	    !write(*,*) months(i)
	  !end do
	  !write(*,*) months
	  open(0,file="AllDust_Dalanzadgad",status="old",action="read")
        read(0,'(A80)') header
         do i=1,10000
	    read(0,5,iostat=iso) oyear(i),omonth(i),oday(i),otime(i),
     &	 utc(i),(opt(i,j),j=1,4),(ssa(i,j),j=1,4),(asy(i,j),j=1,4),
     &      pwc(i),ozone(i),soil(i),water(i),veg(i)
	
           if (iso.gt.0) then 
	      print *, "error while reading Dust data",oyear(i),omonth(i)
     &		  ,oday(i),otime(i) 
              stop
           else if (iso.lt.0) then
	       close (0)
	       exit
           else 
	        ocount=ocount+1
	     end if
		 
	   end do 



c        Make input file for sbdart aerosol model
         do i=1,ocount

          if (omonth(i)(1:1).eq." ") then 
	       omonth(i)(1:1)="0"
	    end if 
          if (oday(i)(1:1).eq." ") then 
	       oday(i)(1:1)="0"
	    end if
          if (otime(i)(1:1).eq." ") then 
	       otime(i)(1:1)="0"
	    end if


c          Write output results
c          ---------------------
           open (1,file="f"//oyear(i)//omonth(i)//oday(i)//otime(i)(1:2)
     &		//otime(i)(4:5),status="new",action="write")
	       do a=1,6
	           if (omonth(i).eq.months(a)) then 
	               idatm=2
	               exit 
	           else
	                idatm=3
	           end if
	        end do

c     Atmospheric profile, wavelength limits, filter function specification etc.

		   write(1,'(1X,A6)') "&INPUT"
		   write(1,'(1X,A9,2X,I1)')      "idatm   =",idatm  
             write(1,'(1X,A9,2X,F5.3)')    "wlinf   =",wlinf
	       write(1,'(1X,A9,2X,F5.3)')    "wlsup   =",wlsup
	       write(1,'(1X,A12)')           "nf      =  1"
	       write(1,'(1X,A12)')           "isat    =  0"
	       write(1,'(1X,A14)')           "wlinc   =  0.0"

c       Solar geometry
            write(1,'(1X,A9,2X,I3)')      "iday    =",int(utc(i))
            write(1,'(1X,A9,2X,F5.2)')    "time    =",
     &		 (utc(i)-real(int(utc(i))))*24.
	      write(1,'(1X,A9,2X,F6.3)')    "alat    =",alat
	      write(1,'(1X,A9,2X,F7.3)')    "alon    =",alon

c     Surface reflectance properties
           write(1,'(1X,A13)')           "isalb   =  10"
	     write(1,'(1X,A9,2X,3(F4.2,A1),F4.2)') "sc      =",snow,",",
     &      water(i),",",soil(i),",",veg(i) 

c      Precipitable water content, ozone etc. 
	   write(1,'(1X,A9,2X,F6.3)')  "uw      =",pwc(i)
	   write(1,'(1X,A9,2X,F6.3)')  "uo3     =",ozone(i)/1000.

c      Stratospheric aerosols
        write(1,'(1X,A14)')        "jaer    =  5*1"    !aerosol=1,aerosol-free=0


c     Boundary layer aerosols
        write(1,'(1X,A9,2X,I1)')    "iaer    =",iaer

	  write(1,'(A37)') "wlbaer  =  0.440, 0.675, 0.870, 1.020"
        write(1,'(1X,A9,1X,4(1X,F6.4,A1))')  "qbaer   =",
     &opt(i,1),",",opt(i,2),",",opt(i,3),",",opt(i,4)
        write(1,'(1X,A9,1X,4(1X,F6.4,A1))') "wbaer   =",ssa(i,1)
     &,",",ssa(i,2),",",ssa(i,3),",",ssa(i,4)
        write(1,'(1X,A9,1X,4(1X,F6.4,A1))') "gbaer   =",asy(i,1)
     &,",",asy(i,2),",",asy(i,3),",",asy(i,4)
	  write(1,'(1X,A13)') "iout    =  10"
        write(1,'(1X,A1)') "/"
	  close(1)

	  end do

    5 FORMAT(A4,T8,A2,T13,A2,1X,A6,1X,F9.5,T32,12(2X,F7.5),T142,F6.4,
     &T149,F5.1,T154,3(1X,F5.3))      

      end 
