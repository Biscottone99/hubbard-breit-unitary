program basis
  implicit none
  integer:: nsiti,i,max,min, check1, check2, nso, count, config, nf, a,b, n, spin, ne, dimm, dimp, dimz, dimm2, dim2
 ! real*8::spin
  character:: array(8)
  logical::bool
  nsiti=4 !siti di disposizione elettroni
  nso=nsiti*2
  ne=4
  
  open(1,file='basis.dat')
  open(2,file='configurations.dat')
  open(3,file='dim2.dat')

  max=0
  do i=nso-1, nso-ne,-1
     max=max+2**i
  enddo
  min=0
  do i=0,ne-1
     min=min+2**i
  enddo


!!!! Starting whit writing of configurations
  nf=0
  do n=min,max
     count=0
     config=0
     a=0
     b=0
     spin=0
     do i=0,nso-1
        bool=btest(n,i)
        if(bool)then
           array(i+1)='1'
           count=count+1
           if(i/2*2.eq.i)then
              a=a+1
           else
              b=b+1
           endif
        else
           array(i+1)='0'          
        endif
     enddo
     spin=(a-b)*0.5d0
     if((count.eq.ne))then
        config=0
        do i=0,nso-1
           if(array(i+1).eq.'1')then
              config=config+2**i
           endif
        enddo
       
        write(1,*) config
        write(2,*) n, (array(i),i=1, nso), spin
        nf=nf+1
     endif
  enddo
     
  
  write(*,*)'numero di funzioni di base uguale', nf
  write(3,*) nf

  

endprogram basis
