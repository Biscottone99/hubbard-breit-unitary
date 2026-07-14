program vb
  use newmodule
  implicit none
  Integer :: i, n, j, nsiti, nso, dim, lwork, lrwork, liwork, info, npoints, nstate, dimm2, dimm, dimz, dimp, dimp2, dimred, state1, state2, multiply
  integer, allocatable :: basis(:), iwork(:),  basism2(:), basism(:), basisz(:), basisp(:), basisp2(:), old(:), nz(:)
  complex*16, allocatable :: hamiltonian(:,:), coup_mono(:,:), soc(:,:), sso(:,:), soo(:,:), work(:), mu(:,:,:), sq(:,:),sqrot(:,:), mono_rot(:,:), sso_rot(:,:), soo_rot(:,:), soc_rot(:,:)
  complex*16 :: imag, pf, unit, temp2
  real*8, allocatable :: coord(:,:), u(:), esite(:), r(:,:), hop_so(:,:), pot(:), hop(:,:), eigenvalue(:), rwork(:), hop_use(:,:), carica(:,:), charges(:,:), spectral(:,:), utility(:)
  real*8 :: length, t, me, gs, e, e0, pi, cl, temp, norm
  logical :: PPPflag, hubbardflag, SOCflag, bool, gslabel, SOCtot, SOCzflag, SOCmono
  complex*16,allocatable::hamm2(:,:), hamm(:,:), hamz(:,:), hamp(:,:), hamp2(:,:), coupx(:,:), coupy(:,:), coupz(:,:), soox(:,:), sooy(:,:), sooz(:,:), ssox(:,:), ssoy(:,:), ssoz(:,:)
  complex*16, allocatable:: coupx_r(:,:), coupy_r(:,:), coupz_r(:,:), soox_r(:,:), sooy_r(:,:), sooz_r(:,:), ssox_r(:,:), ssoy_r(:,:), ssoz_r(:,:), socx_r(:,:), socy_r(:,:), socz_r(:,:)
  character(1) :: jobz, uplo
  character(3) :: hopflag
  character(1), allocatable :: state(:)
  character(10),allocatable:: LABEL(:)
  real*8,allocatable:: ss22(:), sorted(:), vec_out(:),  molt(:), moltmag(:), ordine(:),  hopping_so(:,:), hopping(:,:),dipole(:,:), hop_donor(:,:), hop_acceptor(:,:), hopdr(:,:), hopar(:,:)
  complex*16,allocatable::mat_out(:,:),v1(:), v2(:), sz(:,:), szrot(:,:), num(:,:,:), num_rot(:,:,:), EIG(:), psi0(:), sdr(:,:), muz(:,:), muzrot(:,:),tempo(:)
  complex*16,allocatable:: hop_rot(:,:), hop_cplx(:,:), muy(:,:), muyrot(:,:),mux(:,:), muxrot(:,:), spinpol(:,:), temp_inp(:,:), temp_out(:,:), denmat(:,:)
  !=========================CONSTANTS===========================================================================================================================================================================
  imag = cmplx(0.0, 1.0)
  unit = cmplx(1.0, 0.0)
  me = 9.1093837015d-31
  gs = 2.00231930436256
  e = 1.602176634d-19
  e0 = 8.8541878128d-12
  pi = dacos(-1.0d0)
  cl = 299792458
  pf = 2.25d-2!((gs * e**2) / (8 * pi * e0 * me * cl**2)) * 10.0d10
  write(*,*) pf

  !======================READING INPUT==========================================================================================================================================================================
  open(10,file='dim2.dat')
  read(10,*) dim
  close(10)
  open(1, file='input.inp')
  read(1,*) nsiti
  read(1,*) length ! Armstrong
  read(1,*) t ! eV
  read(1,*) PPPflag
  read(1,*) hubbardflag
  read(1,*) hopflag
  read(1,*) SOCflag
  read(1,*) SOCtot
  read(1,*) socmono
  read(1,*) multiply
  read(1,*) gslabel
  nso = 2 * nsiti
  allocate(coord(nsiti, 3), basis(dim), u(nsiti), esite(nsiti), r(nsiti, nsiti), nz(nsiti))
  open(11,file='output.out')
  write(11,*) 'MODEL PARAMETERS (Site, U, site energi, Z)'
  do i = 1, nsiti
     read(1,*) u(i), esite(i), nz(i)
     write(11,'(I, 2x, 2(f10.5, 2x), I)') i, u(i), esite(i), nz(i)
  end do
  close(1)
  write(11,*) '====================================================================='
  write(11,*) ' GEOMETRY'
  open(2,file='geom.dat')
  write(11,'(A,2x,F5.2,1x,A)') 'Bond length=', length, 'Angstrom'
  write(11,*)
  WRITE(11,*) 'Coordinates (x,y,z)'
  do i = 1, nsiti
     read(2,*) coord(i, 1), coord(i, 2), coord(i, 3)
     write(11,'(I, 2x, 3(f10.5, 2x))') i, coord(i, 1), coord(i, 2), coord(i, 3)
  end do
  close(2)
  ! Computing the distance matrix
  do i = 1, nsiti
     do j = 1, nsiti
        r(i, j) = dsqrt((coord(i, 1) - coord(j, 1))**2 + &
             (coord(i, 2) - coord(j, 2))**2 + &
             (coord(i, 3) - coord(j, 3))**2)
     end do
  end do
  write(11,*)
  write(11,*) 'Distances matrix'
  do i = 1, nsiti
     write(11,'(4(f10.5,2x))') (r(i,j), j = 1, nsiti)
  enddo
  write(11,*) '====================================================================='
  write(11,*)
  write(11,*) 'INPUT PART'
  write(11,*) 'Number of basis functions=', dim
  write(11,*) 'nsiti=', nsiti
  write(11,'(1x,A,2X,F5.2,1x,A)') 't=', t, 'eV'
  write(11,*) 'PPP = ', PPPflag
  write(11,*) 'Hubbard = ', hubbardflag
  write(11,*) 'SOC =', SOCflag
  write(11,*) 'SOC total =', soctot
  write(11,*) 'SOC only z = ', soczflag
  write(11,*) 'SOC only mono =', socmono
  write(11,*) 'Multply factor for SOC, x10^', multiply
  write(11,*) '====================================================================='

  !======================LETTURA BASE===========================================================================================================================================================================
  open(1,file='basis.dat')
  do i=1,dim
     read(1,*) basis(i)
  enddo
  close(1)
!!$  !Ho appena letto la base
!!$  allocate(utility(dim))
!!$
!!$  call sz_on_conf(dim, nso, basis, utility)
!!$  call sort_by_real(basis,utility,dim)
!!$  write(11,*)
!!$  write(11,*) 'Basis functions and their S_z'
!!$  do i = 1, dim
!!$     write(11,'(I2, 2x, I10, 2x, f10.5)') i, basis(i), utility(i)
!!$  enddo
!!$  deallocate(utility)
  !ho appena organizzato la base in funzione di sz crescente
  !=====================================================================================================================
  allocate(num(dim,dim,nso), spinpol(dim,dim))
  num=0
  do n=1,dim
     do i=0,nso-1
        if(btest(basis(n),i))num(n,n,i+1)=num(n,n,i+1)+1
     enddo
  enddo

  spinpol=0d0
  do i = 1, dim
     spinpol(i,i) = (num(i,i,nso-1)-num(i,i,nso))-(num(i,i,1)-num(i,i,2))
  enddo

  allocate(charges(dim,nsiti))
  call charge(charges, basis, nz, dim, nso)
  allocate(sq(dim,dim), sz(dim,dim))
  call s2_realspace(dim, nso, basis, sz, sq)
  !Sulla base ho calcolato l'operatore numero, la spin pol, le cariche, sz e s2
  !====================START MAKING THE OPERATOR FOR DYNAMIC==================================================================================================================================================
!!$  allocate(hopping(nsiti,nsiti))
!!$  hopping=0d0
!!$  hopping(1,2)=t*0.1d0
!!$  hopping(2,1)=t*0.1d0
!!$
!!$  allocate(hopping_so(nso,nso),hop_donor(dim,dim))
!!$  hopping_so = 0d0
!!$  hop_donor=0d0
!!$  call op_siti_2_so(hopping_so, hopping, nso, nsiti)
!!$  call sq_oe_op_real(nso, dim, hopping_so, hop_donor, basis)
!!$  !Ho Definito la matrice di hopping tra il sito 1 e 2 e tra 3 e 4 per la dinamica e l'ho calcolato sulla base RS
!!$  !FINISHED HOPPING FROM DONOR TO B1
!!$  hopping=0d0
!!$  hopping(nsiti-1,nsiti)=t*0.1d0
!!$  hopping(nsiti,nsiti-1)=t*0.1d0
!!$  allocate(hop_acceptor(dim,dim))
!!$  hopping_so = 0d0
!!$  hop_acceptor=0d0
!!$  call op_siti_2_so(hopping_so, hopping, nso, nsiti)
!!$  call sq_oe_op_real(nso, dim, hopping_so, hop_acceptor, basis)
!!$  deallocate(hopping,hopping_so)
!!$  call  check_hermitian(HOP_acceptor*unit, dim, bool)
!!$  if(.not.bool)write(*,*) 'Problem hop acceptor'
!!$  call  check_hermitian(hop_donor*unit, dim, bool)
!!$  if(.not.bool)write(*,*) 'Problem hop donor'
!!$  !FINISHED HOPPING B2 TO DONOR

  !=========================START WRITING THE HAMILTONIAN=======================================================================================================================================================
  allocate( eigenvalue(dim),hamiltonian(dim,dim),pot(dim),hop(dim,dim),dipole(dim,3), muz(dim,dim), mux(dim,dim), muy(dim,dim))
  if(hubbardflag)call site_energy_u(nso, dim, esite, u, basis, pot)
  if(PPPflag) call ppp_diag(dim, nsiti, coord, esite, basis, u, nz, pot)
  !Calcolo la diagonale dell'hamiltoniano
  !====================WRITING OFF-DIAGONAL PART================================================================================================================================================================
  allocate(hop_use(nsiti, nsiti))
  do i = 1, nsiti
     do j = 1, nsiti
        if(hopflag.eq.'NCE')then
           if((i.eq.1).or.(i.eq.nsiti).or.(j.eq.1).or.(j.eq.nsiti))then
              if(i.ne.j)hop_use(i, j) = 0.1*t * dexp(length - r(i,j))
           else
              if(i.ne.j)hop_use(i, j) = t * dexp(length - r(i,j))
           endif
        endif
        if(hopflag.eq.'E')then
           if(i.ne.j)hop_use(i, j) = t * dexp(length - r(i,j))
        endif
        if(hopflag.eq.'D')then
           if(i.ne.j)hop_use(i, j) = t * 10d0**(length - r(i,j))
        endif
        if(hopflag.eq.'N')then
           if( (dabs(r(i,j)-length).le.1d-4).and.(i.ne.j))hop_use(i,j)=t
        endif
        if(hopflag.eq.'NCN')then
           if((i.eq.1).or.(i.eq.nsiti).or.(j.eq.1).or.(j.eq.nsiti))then
             if( (dabs(r(i,j)-length).le.1d-4).and.(i.ne.j)) hop_use(i,j)= 0.1*t
           else
              if( (dabs(r(i,j)-length).le.1d-4).and.(i.ne.j)) hop_use(i,j)= t 
           endif
        endif
     end do
  end do
  write(11,*) 'Hopping'
  do i = 1, nsiti
     write(11,'(<nsiti>(f10.5, 2x))')(hop_use(i,j) , j =1, nsiti)
  enddo
  allocate(hop_so(nso, nso))
  hop_so = 0.0d0
  call op_siti_2_so(hop_so, hop_use, nso, nsiti)
  call sq_oe_op_real(nso, dim, hop_so, hop, basis)
  !Calcolo la parte off-diagonal
  !==============================================================================================================================================================================================================
  call dipole_moment(dipole, charges, coord, dim, nsiti)
  muz=0d0
  mux=0d0
  muy=0d0
  do i = 1, dim
     muz(i,i)= dipole(i,3) * unit
     mux(i,i) = dipole(i,1)*unit
     muy(i,i) = dipole(i,2)*unit     
  enddo
  !=========================COMPUTE SOC========================================================================================================================================================================
  allocate(soc(dim, dim))
  soc = 0.0d0
  if (SOCflag) then

     allocate(coup_mono(dim, dim), sso(dim, dim), soo(dim, dim))
     if(soctot.or.socmono)call compute_soc_mono(nsiti, dim, nz, coord, hop_use, basis, pf, coup_mono)
     if(SOCtot)call compute_sso(nsiti, dim, coord, hop_use, basis, pf, sso)
     if(soctot)call compute_soo(nsiti, dim, coord, hop_use, basis, pf, soo)
     do i = 1, dim
        do j = 1, dim
           if(soctot)soc(i, j) =(coup_mono(i, j)- sso(i, j) - soo(i, j))
           if(socmono)soc(i, j) =coup_mono(i, j)
        end do
     end do

  end if
  if(SOCflag)then
     open(12,file='soc.dat')
     open(13,file='coupling.dat')
     do i = 1, dim
        do j = 1, dim
           write(12, '(I,2x, I, 2x, 10(f10.5, 2x))') i, j, dreal(coup_mono(i,j))*8065.54,dimag(coup_mono(i,j))*8065.54, dreal(soo(i,j))*8065.54,&
               dimag(soo(i,j))*8065.54!, dreal(sso(i,j))*8065.54,dimag(sso(i,j))*8065.54, dreal(soc(i,j))*8065.54,dimag(soc(i,j))*8065.54
           if((zabs(coup_mono(i,j)*8065.54).ge.1d-7).and.(zabs(soc(i,j)*8065.54).le.1d-7))write(13,'(I,2x, I, 2x, 10(f10.5, 2x))') basis(i),basis(j), dreal(coup_mono(i,j))*8065.54,&
                dimag(coup_mono(i,j))*8065.54, dreal(soo(i,j))*8065.54,dimag(soo(i,j))*8065.54, dreal(sso(i,j))*8065.54,dimag(sso(i,j))*8065.54, &
                dreal(soc(i,j))*8065.54,dimag(soc(i,j))*8065.54
        enddo
     end do
     close(12)
     close(13)
  endif

  !=========================TOTAL HAMILTONIAN==================================================================================================================================================================
  hamiltonian=0d0
  do i = 1, dim
     hamiltonian(i,i) = hamiltonian(i,i) + pot(i)
     do j =1, dim
        hamiltonian(i,j) = hamiltonian(i,j) - hop(i,j)+ soc(i,j)*10**multiply
     end do
  end do
  !========================DIAGONALIZATION===============================================================================================================================================================
  lrwork = (1 + 5*dim + 2*dim**2)
  liwork = (3 + 5*dim)
  lwork = (2*dim + dim**2)
  allocate(work(max(1, lwork)), rwork(lrwork), iwork(max(1, liwork)),state(dim))
  call  check_hermitian(hamiltonian, dim, bool)
  if(bool)write(*,*) 'HAMILTONIAN HERMITIAN'
  call zheevd('V', 'U', dim, hamiltonian, dim, eigenvalue, work, lwork, rwork, lrwork, iwork, liwork, info)
  temp = eigenvalue(1)
  eigenvalue = eigenvalue -temp
  !=========================PROPERTIES CALCULATION================================================================================================================================================
  call eigenvalues(dim,1d-8,eigenvalue,state)
  allocate(carica(dim,nsiti))
  call rotate_real(dim, carica, charges, nsiti, hamiltonian)

  allocate(sqrot(dim,dim), szrot(dim,dim), num_rot(dim,dim, nso), muzrot(dim,dim),muxrot(dim,dim), muyrot(dim,dim),sdr(dim,dim))
  sqrot=0d0
  szrot=0d0
  muzrot=0d0
  muxrot=0d0
  muyrot=0d0
  sdr = 0d0
  call rotate_cplx_2x2(dim, sqrot, sq, hamiltonian)
  call rotate_cplx_2x2(dim, szrot, sz, hamiltonian)
  call rotate_cplx_2x2(dim, muzrot, muz, hamiltonian)
  
  call rotate_cplx_2x2(dim, muyrot, muy, hamiltonian)
  call rotate_cplx_2x2(dim, muxrot, mux, hamiltonian)
  allocate(temp_inp(dim,dim), temp_out(dim,dim))
  do i = 1, nso
     temp_inp = num(:,:,i)
     temp_out=0d0
     call rotate_cplx_2x2(dim, temp_out, temp_inp, hamiltonian)
     num_rot(:,:,i) = temp_out(:,:)
  enddo
  deallocate(temp_inp, temp_out)
  !call rotate_cplx_2x2(dim, sdr, spinpol, hamiltonian)
  sdr = 0d0
  sdr = (num_rot(:,:,nso-1)-num_rot(:,:,nso))-(num_rot(:,:,1)-num_rot(:,:,2))

!!$  allocate(hopdr(dim,dim), hopar(dim,dim))
!!$  hopdr=0d0
!!$  hopar=0d0
!!$  call rotate_real_2x2(dim, hopar, hop_acceptor, hamiltonian)
!!$  call rotate_real_2x2(dim, hopdr, hop_donor, hamiltonian)
  open(345, file='number.dat')
  do i = 1, dim
     write(345,'(I2, 2x, <nso>(f17.13,2x))') i, (dreal(num_rot(i,i,j)), j = 1, nso)
  enddo
  close(345)

  open(345, file='soc_check.dat')
  do i = 1, dim
     write(345,*) basis(i)
  enddo
  do i = 1, dim
     do j = 1, dim
        write(345,*) i, j, dreal(soc(i,j)), dimag(soc(i,j))
     enddo
  enddo

  
  allocate(hop_cplx(dim,dim), hop_rot(dim,dim))
  hop_cplx=0d0
  hop_cplx = hop*unit
  call rotate_cplx_2x2(dim, hop_rot, hop_cplx, hamiltonian)
  deallocate(hop, hop_cplx)

  deallocate(muy, mux, sz, sq,num)!, hop_acceptor, hop_donor)
  write(11,*)
  write(11,*) 'EIGENVALUES'
  do i = 1, dim
     write(11,'(I2, 2x, <3>(f17.10, 2x))') i, eigenvalue(i), dreal(sqrot(i,i)), dreal(szrot(i,i))
  enddo

  write(11,*)
  write(11,*) 'CHARGES'
  do i = 1, dim
     write(11,'(I2, 2x, <10>(f10.5, 2x), A)') i, eigenvalue(i), dreal(sqrot(i,i)), (carica(i,j), j = 1, nsiti)
  enddo

  write(11,*)
  write(11,*) 'NUMBERS'
  do i = 1, dim
     write(11,'(I2, 4x, <nso>(f10.5, 2x))') i, (dreal(num_rot(i,i,j)), j = 1, nso)
  enddo
  write(11,*)
  write(11,*) 'TRANSITION DIPOLE MOMENT'
  write(11,*) 'In order: x, y, z'
  do i = 2, dim
     if((zabs(muxrot(1,i)).ge.1d-5).or.(zabs(muyrot(1,i)).ge.1d-5).or.(zabs(muzrot(1,i)).ge.1d-5)) write(11,'(I2, 4x, <nso>(f10.5, 2x))') i, zabs(muxrot(1,i))**2, zabs(muyrot(1,i))**2, zabs(muzrot(1,i))**2
  enddo

  write(11,*)
  write(11,*) 'SPIN-DENSITY'

  do i = 1, dim
     write(11,'(I2, 2x, <3>(f10.5, 2x),E10.2)') i, eigenvalue(i), dreal(sqrot(i,i)), dreal(szrot(i,i)),(dreal(sdr(i,i)))
  enddo
    
  !=========================REDFIELD===========================================================================================================================================================================
  allocate(psi0(dim))
  do i = 1, dim
     psi0(i)=muzrot(1,i) !sbagliato
    ! psi0(i)=muzrot(i,1)
  enddo
  norm=0.d0
  do i=1,dim
     norm=norm+dconjg(psi0(i))*psi0(i)
  enddo
  psi0=psi0/dsqrt(norm)
  temp=0d0

 
  temp = 0d0
  do i = 1,dim
     temp = temp + zabs(psi0(i))**2 * eigenvalue(i)
  enddo
  write(*,*) 'Energy (eV)=', temp
  temp = 0d0
  do i = 1,dim
     temp = temp + zabs(psi0(i))*sdr(i,i)
  enddo
  write(*,*) 'Spin density=', temp
  allocate(denmat(dim,dim))
  !normalization
  temp = 0.d0
  do i = 1,dim
     temp = temp + dconjg(psi0(i))*psi0(i)
  end do
  temp = dsqrt(temp)
  psi0 = psi0/temp

  !writing the density matrix at t=0
  denmat = (0.d0,0.d0)
  do i = 1,dim
     do j = 1,dim
        denmat(i,j) = dconjg(psi0(i))*psi0(j)
     end do
  end do
  open(66,file='denmat.dat')
  do i = 1, dim
     do j = 1, dim
        write(66,*) i,j,denmat(i,j)
     enddo
  enddo
  if(.not.gslabel)then
     open(66,file='unitary/rho.bin',form="unformatted", access='stream')
     write(66) denmat(2:dim,2:dim)
     close(66)
     open(66,file='unitary/spin-density.bin',form="unformatted", access='stream')
     write(66) sdr(2:dim,2:dim)
     close(66)
     open(66,file='unitary/eigen.bin',form="unformatted", access='stream')
     write(66) eigenvalue(2:dim)
     close(66)
!!$     open(88,file='unitary/psi0.bin',form="unformatted")
!!$     write(88) psi0(2:dim)
!!$     close(88)
!!$     open(99,file='unitary/num.bin',form="unformatted")
!!$     write(99) num_rot(2:dim,2:dim,:)
!!$     close(99)
!!$     open(99,file='unitary/op1.bin',form="unformatted")
!!$     write(99) hopdr(2:dim,2:dim)*unit
!!$     close(99)
!!$     open(99,file='unitary/op2.bin',form="unformatted")
!!$     write(99) hopar(2:dim,2:dim)*unit
!!$     close(99)

!!$     open(55,file='unitary/system_input.dat')
!!$     write(55,*) dim-1
!!$     write(55,*) nsiti
  endif


 
  if(gslabel)then
     open(66,file='unitary/rho.bin',form="unformatted", access='stream')
     write(66) denmat
     close(66)
     open(66,file='unitary/spin-density.bin',form="unformatted", access='stream')
     write(66) sdr(1:dim,1:dim)
     close(66)
     open(66,file='unitary/eigen.bin',form="unformatted", access='stream')
     write(66) eigenvalue(1:dim)
     close(66)
!!$     open(66)
!!$     write(66) sqrot(1:dim,1:dim)
!!$     close(66)
!!$     open(88,file='unitary/psi0.bin',form="unformatted")
!!$     write(88) psi0(1:dim)
!!$     close(88)
!!$     open(99,file='unitary/num.bin',form="unformatted")
!!$     write(99) num_rot(1:dim,1:dim,:)
!!$     close(99)
!!$     open(99,file='unitary/op1.bin',form="unformatted")
!!$     write(99) hopdr(1:dim,1:dim)*unit
!!$     close(99)
!!$     open(99,file='unitary/op2.bin',form="unformatted")
!!$     write(99) hopar(1:dim,1:dim)*unit
     close(99)
!!$
!!$     open(55,file='unitary/system_input.dat')
!!$     write(55,*) dim
!!$     write(55,*) nsiti
  endif

  write(11,*)
  write(11,*) 'DYNAMIC INPUT'
  write(11,*) 'GS in Dynamic=', gslabel
  write(11,*) 'States number in dynamic=', DIMRED
end program vb
