module newmodule
  implicit none
contains
subroutine dipole_moment(dipole, carica, nuclei, dim2, nsiti)
    implicit none
    integer, intent(in)::dim2, nsiti
    real*8, intent(in):: carica(dim2,nsiti), nuclei(nsiti,3)
    real*8, intent(out):: dipole(dim2,3)
    integer::n,k,j

    dipole = 0.0d0

    do n = 1, dim2
       do k = 1, 3
          do j = 1, nsiti
             dipole(n, k) = dipole(n, k) + carica(n, j) * nuclei(j, k)
          enddo
       enddo
    enddo
  end subroutine dipole_moment

  subroutine check_hermitian(matrix, n, is_hermitian)
    implicit none
    integer, intent(in) :: n
    complex*16, intent(in) :: matrix(n, n)
    logical, intent(out) :: is_hermitian

    integer :: i, j
    complex*16, allocatable:: conj_transpose_matrix(:, :)
    allocate(conj_transpose_matrix(n, n))
    conj_transpose_matrix=0

    ! Calcola la trasposta coniugata della matrice
    do i = 1, n
       do j = 1, n
          conj_transpose_matrix(i, j) = dconjg(matrix(j, i))
       end do
    end do

    ! Verifica se la matrice è Hermitiana
    is_hermitian = all(zabs(matrix - conj_transpose_matrix).le.1d-10)

  end subroutine check_hermitian

  subroutine charge(carica, vecconfig, nz, dim2, nso)
    implicit none
    integer, intent(in):: dim2, nso, vecconfig(dim2), nz(nso/2)
    real*8, intent(out):: carica(dim2,nso/2)
    integer:: n, i, sito,a , b
    logical::bool, bool1

    carica=0
    do n=1,dim2
       do i=0,nso-2,2
          sito=(i+2)/2
          bool=btest(vecconfig(n), i)
          bool1=btest(vecconfig(n), i+1)
          if(bool)then
             a=1
          else
             a=0
          endif

          if(bool1)then
             b=1
          else
             b=0
          endif
          carica(n,sito)=nz(sito)-(a+b)
       enddo
    enddo
  end subroutine charge

  subroutine s2_realspace(dim, nso, basis, szo, sq)
    integer, intent(in) :: dim, nso, basis(dim)
    complex*16, intent(out) :: sq(dim, dim), szo(dim,dim)

    integer :: i, j
    complex*16 :: cplx, spin(2,2,3)
    complex*16, allocatable :: sx(:,:), sy(:,:), sz(:,:), sxrs(:,:), syrs(:,:), szrs(:,:)
    logical :: is_hermitian

    ! Definisco cplx come unità immaginaria
    cplx = cmplx(0.d0, 1.d0)

    ! Inizializzo la matrice degli operatori di spin
    spin = 0
    spin(1,2,1) = 1.d0
    spin(2,1,1) = 1.d0

    spin(1,2,2) = -cplx
    spin(2,1,2) = cplx

    spin(1,1,3) = 1.d0
    spin(2,2,3) = -1.d0

    !========================= Operatori di spin =========================
    allocate(sx(nso, nso), sy(nso, nso), sz(nso, nso))
    sx = 0
    sy = 0
    sz = 0

    ! Costruisco l'operatore sx
    do i = 1, nso - 1, 2
       sx(i,i) = spin(1,1,1)
       sx(i+1,i+1) = spin(2,2,1)
       sx(i,i+1) = spin(1,2,1)
       sx(i+1,i) = spin(2,1,1)
    end do

    ! Costruisco l'operatore sy
    do i = 1, nso - 1, 2
       sy(i,i) = spin(1,1,2)
       sy(i+1,i+1) = spin(2,2,2)
       sy(i,i+1) = spin(1,2,2)
       sy(i+1,i) = spin(2,1,2)
    end do

    ! Costruisco l'operatore sz
    do i = 1, nso - 1, 2
       sz(i,i) = spin(1,1,3)
       sz(i+1,i+1) = spin(2,2,3)
       sz(i,i+1) = spin(1,2,3)
       sz(i+1,i) = spin(2,1,3)
    end do

    !========================= Trasformazione in real space =========================
    allocate(sxrs(dim, dim), syrs(dim, dim), szrs(dim,dim))

    ! Trasformazione di sx
    call sq_oe_op_compl(nso, dim, sx, sxrs, basis)
    call check_hermitian(sxrs, dim, is_hermitian)
    if (.not. is_hermitian) write(*,*) 'Problem sxrs'

    ! Trasformazione di sy
    call sq_oe_op_compl(nso, dim, sy, syrs, basis)
    call check_hermitian(syrs, dim, is_hermitian)
    if (.not. is_hermitian) write(*,*) 'Problem syrs'

    ! Trasformazione di sz
    call sq_oe_op_compl(nso, dim, sz, szrs, basis)
    call check_hermitian(szrs, dim, is_hermitian)
    if (.not. is_hermitian) write(*,*) 'Problem szrs'

    !========================= Calcolo degli operatori al quadrato =========================
    sxrs = 0.5 * sxrs
    call square_complex_matrix(dim, sxrs)

    syrs = 0.5 * syrs
    call square_complex_matrix(dim, syrs)

    szrs = 0.5 * szrs
    do i = 1, dim
       szo(i,i) = szrs(i,i)
    enddo
    call square_complex_matrix(dim, szrs)

    ! Controllo se le matrici al quadrato sono hermitiane
    call check_hermitian(sxrs, dim, is_hermitian)
    if (.not. is_hermitian) write(*,*) 'Not hermitian sx^2'

    call check_hermitian(syrs, dim, is_hermitian)
    if (.not. is_hermitian) write(*,*) 'Not hermitian sy^2'

    call check_hermitian(szrs, dim, is_hermitian)
    if (.not. is_hermitian) write(*,*) 'Not hermitian sz^2'

    !========================= Calcolo di S^2 =========================
    sq = 0
    do i = 1, dim
       do j = 1, dim
          sq(i,j) = sq(i,j) + sxrs(i,j) + syrs(i,j) + szrs(i,j)
       end do
    end do

  end subroutine s2_realspace


  subroutine op_siti_2_so(hop_so, hop_use, nso, nsiti)
    implicit none
    integer, intent(in) :: nsiti,nso
    real*8, intent(in) :: hop_use(nsiti, nsiti)
    real*8, intent(out) :: hop_so(nso, nso)
    integer :: i, j, isito, jsito

    hop_so = 0.0d0
    do i = 1, nso
       do j = 1, nso
          isito = (i + 1) / 2
          jsito = (j + 1) / 2
          if ((hop_use(isito, jsito) /= 0.0d0) .and. (mod(i, 2) == mod(j, 2))) then
             hop_so(i, j) = hop_use(isito, jsito)
          end if
       end do
    end do
  end subroutine op_siti_2_so

  subroutine site_energy_u(nso, dim2, esite, u, vecconfig, energy)
    implicit none
    integer :: n, i, sito
    real*8, intent(in) :: esite(nso/2), u(nso/2)
    integer, intent(in) :: dim2, vecconfig(dim2), nso
    real*8, intent(out) :: energy(dim2)

    ! Initialize energy array to zero
    energy = 0.0d0

    do n = 1, dim2
       ! Add single site energies
       do i = 0, nso-1
          sito = (i+2)/2
          if (btest(vecconfig(n), i)) energy(n) = energy(n) + esite(sito)
       end do

       ! Add interaction energies
       do i = 0, nso-2, 2
          sito = (i+2)/2
          if (btest(vecconfig(n), i) .and. btest(vecconfig(n), i+1)) energy(n) = energy(n) + u(sito)
       end do
    end do
  end subroutine site_energy_u


  subroutine ppp_diag(dim2, nsiti, nuclei, esite, vecconfig, u, nz, pot)
    implicit none
    integer, intent(in) :: dim2, nsiti
    integer, intent(in) :: vecconfig(dim2), nz(nsiti)
    real*8, intent(in) :: nuclei(nsiti, 3), u(nsiti), esite(nsiti)
    real*8, intent(out) :: pot(dim2)
    integer :: n, i, j, p, sito, occupazioni(nsiti), a, b
    real*8 :: PPP, r(nsiti, nsiti), dx, dy, dz

    ! Initialize pot array to zero
    pot = 0.0d0
    ! Initialize r array to zero
    r = 0.0d0

    do i = 1, nsiti
       do j = i + 1, nsiti
          dx = nuclei(i, 1) - nuclei(j, 1)
          dy = nuclei(i, 2) - nuclei(j, 2)
          dz = nuclei(i, 3) - nuclei(j, 3)

          r(i, j) = dsqrt(dx**2 + dy**2 + dz**2)
          r(j, i) = r(i, j)
       end do
    end do

    do n = 1, dim2
       sito = 0
       do i = 0, 2 * nsiti - 2, 2
          sito = (i + 2) / 2
          a = 0
          b = 0
          if (btest(vecconfig(n), i)) a = 1
          if (btest(vecconfig(n), i + 1)) b = 1
          occupazioni(sito) = a + b
       end do

       do i = 1, nsiti
          pot(n) = pot(n) + esite(i) * occupazioni(i)
       end do

       PPP = 0.0d0
       do i = 1, nsiti
          do p = 1, nsiti
             if (i /= p) PPP = PPP + (14.397d0 / dsqrt(r(i, p)**2 + (28.794d0 / (u(i) + u(p)))**2)) * (nz(i) - occupazioni(i)) * (nz(p) - occupazioni(p))
             if ((i == p) .and. (occupazioni(i) == 2)) PPP = PPP + 2.d0 * u(i)
          end do
       end do
       pot(n) = pot(n) + 0.5d0 * PPP
    end do
  end subroutine ppp_diag


  subroutine sq_oe_op_real(nso,dim,op_tb,op,basis)
    implicit none
    integer :: iso,jso,conta,i,j,istate,jstate,step,a
    integer, intent(in) :: dim,nso,basis(dim)
    double precision, intent(in) :: op_tb(nso,nso)
    double precision, intent(out) :: op(dim,dim)
    double precision :: phase

    op = 0.d0
    !$omp parallel do default(none) private(iso,jso,conta,i,j,istate,jstate,step,a,phase) shared(dim,nso,basis,op_tb,op)
    do j = 1,dim !col index
       do iso = 0,nso-1 ! creation op index
          do jso = 0,nso-1 ! annihilation op index
             jstate = basis(j)
             if (btest(jstate,jso)) then
                istate = ibclr(jstate,jso)
                if (.not.btest(istate,iso)) then

                   istate = ibset(istate,iso)

                   i = binary_search(basis, istate, 1, dim) !row index
                   if (i/=0) then
                      !determine the phase
                      !get direction from iso to jso
                      if (jso>iso) step = -1
                      if (iso>jso) step = 1

                      if (iso==jso) then
                         phase = 1.d0
                      else
                         conta = 0
                         do a = jso+step, iso-step, step
                            if (btest(istate,a)) conta = conta + 1
                         end do

                         if (mod(conta,2)==0) then
                            phase = 1.d0
                         else
                            phase = -1.d0
                         end if
                      end if

                      op(i,j) = op(i,j) + phase * op_tb(iso+1,jso+1)
                   end if
                end if
             end if

          end do
       end do
    end do
    !$omp end parallel do
  end subroutine sq_oe_op_real

  !generates a DOUBLE COMPLEX one-electron operator in second quantization
  !given a real space basis and one-electron operator in first quantization
  !nso    number of spin-orbitals
  !dim    real-space dimension
  !op_tb  operator in the first quantization (tight binding basis)
  !op     operator in the real space basis
  !basis  1D array containing the integers that describe the real space basis in bit representaion
  !NB: the same spin-orbital ordering used to describe the real space configurations
  !    must be used for the first quantization basis

  subroutine sq_oe_op_compl(nso,dim,op_tb,op,basis)
    implicit none
    integer iso,jso,conta,i,j,istate,jstate,step,a
    integer, intent (in) :: dim,nso,basis(dim)
    double complex, intent (in) :: op_tb(nso,nso)
    double complex, intent (out) :: op(dim,dim)
    double precision phase

    op = 0.d0
    !$omp parallel do default(none) private(iso,jso,conta,i,j,istate,jstate,step,a,phase) shared(dim,nso,basis,op_tb,op)
    do j = 1,dim !col index
       do iso = 0,nso-1 ! creation op index
          do jso = 0,nso-1 ! annihilation op index
             jstate = basis(j)
             if (btest(jstate,jso)) then
                istate = ibclr(jstate,jso)
                if (.not.btest(istate,iso)) then

                   istate = ibset(istate,iso)

                   i = binary_search(basis, istate, 1, dim) !row index
                   if (i/=0) then
                      !determine the phase
                      !get direction from iso to jso
                      if (jso>iso) step = -1
                      if (iso>jso) step = 1

                      if (iso==jso) then
                         phase = 1.d0
                         goto 1000
                      end if

                      conta = 0
                      do a = jso+step, iso-step, step
                         if (btest(istate,a)) conta = conta + 1
                      end do

                      if (conta/2*2==conta) then
                         phase = 1.d0
                      else
                         phase = -1.d0
                      end if

1000                  continue

                      !write(*,*) phase,i,j,(btest(istate,a),a=0,nso-1), 0,0, (btest(jstate,a),a=0,nso-1)
                      op(i,j) = op(i,j) + phase * op_tb(iso+1,jso+1)
                   end if
                end if
             end if

          end do
       end do
       !write(*,*) ''

    end do
    !$omp end parallel do
  end subroutine sq_oe_op_compl



  subroutine eigenvalues(dim2,tollerance,w,state)


    integer, intent(in) :: dim2
    real*8, intent(in) :: w(dim2)
    real*8, intent(in) :: tollerance
    character(len=1), intent(out) :: state(dim2)

    integer :: i, count_result, j

    do i = 1, dim2
       ! Inizializza il contatore
       count_result = 0

       ! Conta gli elementi che soddisfano la condizione
       do j = 1, dim2
          if (abs(w(j) - w(i)) < tollerance) then
             count_result = count_result + 1
          endif
       end do

       ! Determina lo stato in base al risultato del conteggio
       if (count_result == 3) then
          state(i) = 'T'
       elseif (count_result == 5) then
          state(i) = 'Q'
       elseif (count_result == 4) then
          state(i) = '4'
       elseif (count_result == 2) then
          state(i) = 'D'
       else
          state(i) = 'S'
       endif
    end do
  end subroutine eigenvalues




  subroutine rotate_real(dim2, out, in, col, rot)
    implicit none
    integer :: i, j, k
    integer, intent(in) :: dim2, col
    complex*16, intent(in) :: rot(dim2, dim2)
    real*8, intent(in) :: in(dim2, col)
    real*8, intent(out) :: out(dim2, col)

    !$omp parallel do default(none) private(i, j, k) shared(dim2, col, rot, in, out)
    do i = 1, dim2
       do j = 1, dim2
          do k = 1, col
             out(i, k) = out(i, k) + dconjg(rot(j, i)) * rot(j, i) * in(j, k)
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine rotate_real

  subroutine rotate_cplx(dim2, out, in, col, rot)
    implicit none
    integer :: i, j, k
    integer, intent(in) :: dim2, col
    complex*16, intent(in) :: rot(dim2, dim2)
    complex*16, intent(in) :: in(dim2, col)
    complex*16, intent(out) :: out(dim2, col)

    !$omp parallel do default(none) private(i, j, k) shared(dim2, col, rot, in, out)
    do i = 1, dim2
       do j = 1, dim2
          do k = 1, col
             out(i, k) = out(i, k) + dconjg(rot(j, i)) * rot(j, i) * in(j, k)
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine rotate_cplx



  function binary_search(array, target, usefull, n) result(pos)
    implicit none
    integer, intent(in) :: array(:), target, usefull, n
    integer :: pos
    integer :: i

    ! Inizializza la posizione come -1 (non trovato)
    pos = 0

    ! Ciclo attraverso il vettore per trovare il target
    do i = 1, n
       if (array(i) == target) then
          pos = i
          return
       end if
    end do
  end function binary_search


  subroutine rotate_cplx_2x2(dim2, coupling, coup, ham)
    implicit none
    integer :: i, l, j, k
    integer, intent(in) :: dim2
    complex*16, intent(out) :: coupling(dim2, dim2)
    complex*16, intent(in) :: coup(dim2, dim2)
    complex*16, intent(in) :: ham(dim2, dim2)

    !$omp parallel do default(none) private(i, l, j, k) shared(coupling, coup, ham, dim2)
    do i = 1, dim2
       do l = 1, dim2
          do j = 1, dim2
             do k = 1, dim2
                coupling(i, l) = coupling(i, l) + dconjg(ham(j, i)) * ham(k, l) * coup(j, k)
             end do
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine rotate_cplx_2x2


  subroutine rotate_real_2x2(dim2, coupling, coup, ham)
    implicit none
    integer :: i, l, j, k
    integer, intent(in) :: dim2
    real*8, intent(out) :: coupling(dim2, dim2)
    real*8, intent(in) :: coup(dim2, dim2)
    complex*16, intent(in) :: ham(dim2, dim2)

    !$omp parallel do default(none) private(i, l, j, k) shared(coupling, coup, ham, dim2)
    do i = 1, dim2
       do l = 1, dim2
          do j = 1, dim2
             do k = 1, dim2
                coupling(i, l) = coupling(i, l) + dconjg(ham(j, i)) * ham(k, l) * coup(j, k)
             end do
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine rotate_real_2x2


  subroutine compute_soc_mono(nsiti, dim2,nz, nuclei, hop, vecconfig, pf, coup)
    implicit none
    integer, intent(in) :: nsiti, dim2, nz(nsiti)
    real*8, intent(in) :: nuclei(nsiti, 3), hop(nsiti, nsiti)
    integer, intent(in) :: vecconfig(dim2)
    complex*16, intent(in) :: pf
    complex*16, intent(out) :: coup(dim2, dim2)
    complex*16 :: coupx(dim2,dim2), coupy(dim2,dim2), coupz(dim2,dim2)
    integer :: i, j, k, isito, sitoi, sitoj, si, sj, nso
    complex*16 :: mom(nsiti, nsiti, 3)
    complex*16,allocatable:: soc_a(:,:,:), soc_b(:,:,:), soc_mono(:,:,:), hamsoc(:,:), tbcoupx(:,:), tbcoupy(:,:), tbcoupz(:,:)
    real*8 :: dr,  radius

    logical :: bool, is_hermitian
    complex*16 :: spin(2, 2, 3), cplx,vec1(3), vec2(3), cp(3)


    ! Define number of spin orbitals
    nso = 2 * nsiti
    allocate( soc_a(3, nso, nso), soc_b(3, nso, nso), soc_mono(3, nso, nso), hamsoc(nso, nso))
    allocate(tbcoupx(nso,nso), tbcoupy(nso,nso), tbcoupz(nso,nso))
    ! Define the complex unit
    cplx = (0.0d0, 1.0d0)

    ! Define the spin matrices
    spin = 0
    spin(1,2,1) = 1d0
    spin(2,1,1) = 1d0

    spin(1,2,2) = -cplx
    spin(2,1,2) =  cplx

    spin(1,1,3) = 1d0
    spin(2,2,3) = -1d0



    ! Initialize the mom array
    mom = (0.0d0, 0.0d0)
    do i = 1, nsiti
       do j = 1, nsiti
          do k = 1, 3
             dr = nuclei(i, k) - nuclei(j, k)
             mom(i, j, k) = cplx * dr * hop(i, j)
          end do
       end do
    end do

    ! Check if mom is Hermitian
    do k = 1, 3
       call check_hermitian(mom(:, :, k), nsiti, bool)
       if (.not. bool) write(*, *) 'PROBLEMI MOM K=', k
    end do

    ! Initialize soc_a, soc_b, soc_mono arrays to zero
    soc_a = (0.0d0, 0.0d0)
    soc_b = (0.0d0, 0.0d0)
    soc_mono = (0.0d0, 0.0d0)

    ! Calculate soc_a
    do i = 1, nso !elettroni
       do j = 1, nso !elettroni
          do isito = 1, nsiti !nuclei
             if (isito /= (i + 1) / 2) then    !isito è il nucleo
                vec1 = 0.0d0
                vec2 = 0.0d0
                cp = 0.0d0
                sitoi = (i + 1) / 2 !elettroni
                sitoj = (j + 1) / 2 !elettroni
                do k = 1, 3
                   vec1(k) = nuclei(sitoi, k) - nuclei(isito, k)
                   vec2(k) = mom(sitoi, sitoj, k)
                end do
                cp = cross_product(vec1, vec2)

                radius = 0.0d0
                do k = 1, 3
                   radius = radius + dreal(vec1(k))**2
                end do
                radius = (dsqrt(radius))**3
                cp = cp*nz(isito) / radius

                if (mod(i, 2) == 0) then !check
                   si = 2
                else
                   si = 1
                end if

                if (mod(j, 2) == 0) then
                   sj = 2
                else
                   sj = 1
                end if

                do k = 1, 3
                   soc_a(k, i, j) = soc_a(k, i, j) + cp(k) * spin(si, sj, k)
                end do
             end if
          end do
       end do
    end do

    ! Calculate soc_b
    do i = 1, nso
       do j = 1, nso
          do isito = 1, nsiti
             if (isito /= (j + 1) / 2) then
                vec1 = 0.0d0
                vec2 = 0.0d0
                cp = 0.0d0
                sitoi = (i + 1) / 2
                sitoj = (j + 1) / 2
                do k = 1, 3
                   vec1(k) = nuclei(sitoj, k) - nuclei(isito, k)
                   vec2(k) = mom(sitoi, sitoj, k)
                end do
                cp  = cross_product(vec1, vec2)

                radius = 0.0d0
                do k = 1, 3
                   radius = radius + dreal(vec1(k))**2
                end do
                radius = (dsqrt(radius))**3
                cp = cp * nz(isito) / radius

                if (mod(i, 2) == 0) then
                   si = 2
                else
                   si = 1
                end if

                if (mod(j, 2) == 0) then
                   sj = 2
                else
                   sj = 1
                end if

                do k = 1, 3
                   soc_b(k, i, j) = soc_b(k, i, j) + cp(k) * spin(si, sj, k)
                end do
             end if
          end do
       end do
    end do

    ! Calculate soc_mono
    do i = 1, nso
       do j = 1, nso
          do k = 1, 3
             soc_mono(k, i, j) = 0.5d0 * (soc_a(k, i, j) + soc_b(k, i, j))
          end do
       end do
    end do

    ! Initialize and calculate hamsoc
    hamsoc = (0.0d0, 0.0d0)
    do i = 1, nso
       do j = 1, nso
          do k = 1, 3
             hamsoc(i, j) = hamsoc(i, j) + soc_mono(k, i, j)
          end do
       end do
    end do
    open(1111,file='soc_mono_so.dat')
    do i = 1, nso
       do j = 1, nso
          write(1111,'(I5,I5,2ES20.12)') i, j, dreal(hamsoc(i,j)), dimag(hamsoc(i,j))
       end do
    end do
    coupx=0
    coupy=0
    coupz=0
    tbcoupx = soc_mono(1, :, :)
    tbcoupy = soc_mono(2, :, :)
    tbcoupz = soc_mono(3, :, :)



    ! Calculate coup
    coup = 0.0d0
    call sq_oe_op_compl(nso, dim2, hamsoc, coup, vecconfig)
    call check_hermitian(coup, dim2, is_hermitian)
    if (is_hermitian) write(*, *) 'COUP HERMITIANA'

    call sq_oe_op_compl(nso, dim2, tbcoupx, coupx, vecconfig)
    call check_hermitian(coupx, dim2, is_hermitian)
    if (.not.is_hermitian) write(*, *) 'COUPX Problem'

    call sq_oe_op_compl(nso, dim2, tbcoupy, coupy, vecconfig)
    call check_hermitian(coupy, dim2, is_hermitian)
    if (.not.is_hermitian) write(*, *) 'COUPY Problem'

    call sq_oe_op_compl(nso, dim2, tbcoupz, coupz, vecconfig)
    call check_hermitian(coupz, dim2, is_hermitian)
    if (.not.is_hermitian) write(*, *) 'COUPZ Problem'

    coup = pf * coup
    coupx = pf * coupx
    coupy = pf * coupy
    coupz = pf * coupz

  end subroutine compute_soc_mono


  subroutine compute_sso(nsiti, dim2, nuclei, hop, vecconfig, pf, sso)
    implicit none
    integer, intent(in) :: nsiti, dim2
    real*8, intent(in) :: nuclei(nsiti, 3), hop(nsiti, nsiti)
    integer, intent(in) :: vecconfig(dim2)
    complex*16, intent(in) :: pf
    complex*16, intent(out) :: sso(dim2, dim2)
    complex*16 :: ssox(dim2, dim2), ssoy(dim2, dim2), ssoz(dim2, dim2)

    complex*16,allocatable:: hssotb(:,:,:,:,:), ssotb(:,:,:,:), ppso(:,:,:),ssotbx(:,:,:,:),ssotby(:,:,:,:),ssotbz(:,:,:,:), hssotb2(:,:,:,:,:)
    complex*16 ::  spin(2, 2, 3), cplx, vec1(3), vec2(3), cp(3), cp2(3)
    real*8 :: dist(nsiti, nsiti, 3), radius
    integer :: i, j, k, a, b, c, d, si, sj, nso, asito, bsito, csito, dsito
    logical :: bool

    nso = 2 * nsiti
    allocate(hssotb(3, nso, nso, nso, nso), ssotb(nso, nso, nso, nso), ppso(3, nsiti, nsiti),hssotb2(3, nso, nso, nso, nso))
    allocate(ssotbx(nso, nso, nso, nso),ssotby(nso, nso, nso, nso),ssotbz(nso, nso, nso, nso))

    ! Define the complex unit
    cplx = (0.0d0, 1.0d0)

    ! Define the spin matrices
    spin = 0
    spin(1, 2, 1) = 1d0
    spin(2, 1, 1) = 1d0

    spin(1, 2, 2) = -cplx
    spin(2, 1, 2) = cplx

    spin(1, 1, 3) = 1d0
    spin(2, 2, 3) = -1d0

    ! Allocate and initialize the dist array
    ppso=0d0
    do k=1,3
       do i=1,nsiti
          do j=1,nsiti
             ppso(k,i,j)=cplx*hop(i,j)*(nuclei(i,k)-nuclei(j,k))
          enddo
       enddo
    enddo

    ! Initialize hssotb array
    hssotb = (0.0d0, 0.0d0)

    ! Calculate hssotb
    do a = 1, nso
       do b = 1, nso
          do c = 1, nso
             asito = (a+1)/2
             bsito = (b+1)/2
             csito = (c+1)/2

             if (asito.ne.bsito)then

                vec1 = 0.0d0
                vec2 = 0.0d0
                cp = 0.0d0
                do k = 1, 3
                   vec1(k) = nuclei(asito, k) - nuclei(bsito, k) !r_ab
                   vec2(k) = ppso(k, asito, csito) !p_ac
                end do
                cp  = cross_product(vec1, vec2)

                radius = 0.0d0
                do k = 1, 3
                   radius = radius + dreal(vec1(k))**2
                end do
                radius = (dsqrt(radius))**3
                cp = cp / radius
             else
                cp =0d0
             endif

             if(bsito.ne.csito)then
                vec1 = 0.0d0
                vec2 = 0.0d0
                cp2 = 0.0d0
                do k = 1, 3
                   vec2(k) = nuclei(csito, k) - nuclei(bsito, k) !r_cb
                   vec1(k) = ppso(k, asito, csito) !p_ac
                end do
                cp2  = cross_product(vec1, vec2)

                radius = 0.0d0
                do k = 1, 3
                   radius = radius + dreal(vec2(k))**2
                end do
                radius = (dsqrt(radius))**3
                cp2 = cp2 / radius
             else
                cp2=0d0
             endif

             if (mod(a, 2) == 0) then  !spin_a
                si = 2
             else
                si = 1
             end if

             if (mod(c, 2) == 0) then !spin_b
                sj = 2
             else
                sj = 1
             end if

             do k = 1, 3
                hssotb(k, a, b, c, b) = hssotb(k, a, b, c, b) + (cp(k)-cp2(k)) * spin(si, sj, k)
             end do

          end do
       end do
    end do
    ! Initialize ssotb array
    ssotb = (0.0d0, 0.0d0)

    ! Calculate ssotb

    do a = 1, nso
       do b = 1, nso
          do c = 1, nso
             do d = 1, nso
                do k = 1, 3
                   ssotb(a, b, c, d) = ssotb(a, b, c, d) + 0.5d0 * (hssotb(k, a, b, c, d))
                end do
                !  if(zabs(ssotb(a,b,c,d)).ge.1d-10) write(1,*) a, b, c, d, ssotb(a, b, c, d)
             end do
          end do
       end do
    end do
    do a = 1, nso
       do b = 1, nso
          do c = 1, nso
             do d = 1, nso
                ssotbx(a, b, c, d) = ssotbx(a, b, c, d) + 0.5d0 * (hssotb(1, a, b, c, d))
                ssotby(a, b, c, d) = ssotby(a, b, c, d) + 0.5d0 * (hssotb(2, a, b, c, d) )
                ssotbz(a, b, c, d) = ssotbz(a, b, c, d) + 0.5d0 * (hssotb(3, a, b, c, d) )
                ! if(zabs(ssotbz(a,b,c,d)).ge.1d-10) write(1,*) a, b, c, d, ssotbz(a, b, c, d)
             end do
          end do
       end do
    end do

    ! Initialize and calculate sso
    sso = 0.0d0
    call bielectron(dim2, nso, pf, vecconfig, ssotb, sso)
    call bielectron(dim2, nso, pf, vecconfig, ssotbx, ssox)
    call bielectron(dim2, nso, pf, vecconfig, ssotby, ssoy)
    call bielectron(dim2, nso, pf, vecconfig, ssotbz, ssoz)

    call check_hermitian(sso, dim2, bool)
    if ( bool) write(*, *) 'SSO Hermitian'

    call check_hermitian(ssoy, dim2, bool)
    if ( .not.bool) write(*, *) 'SSOY Problem'

    call check_hermitian(ssox, dim2, bool)
    if ( .not.bool) write(*, *) 'SSOX Problem'

    call check_hermitian(ssoz, dim2, bool)
    if ( .not.bool) write(*, *) 'SSOZ Problem'
  end subroutine compute_sso



  subroutine compute_soo(nsiti, dim2, nuclei, hop, vecconfig, pf, soo)
    implicit none
    integer, intent(in) :: nsiti, dim2
    real*8, intent(in) :: nuclei(nsiti, 3), hop(nsiti,nsiti)
    complex*16, intent(in) :: pf
    integer, intent(in) :: vecconfig(dim2)
    complex*16 :: soo(dim2, dim2), soox(dim2,dim2), sooy(dim2,dim2), sooz(dim2,dim2)

    complex*16, allocatable:: hsootb(:,:,:,:,:), sootb(:,:,:,:), ppso(:,:,:), check(:,:),  conj_transpose_matrix(:, :, :, :), hsootb2(:,:,:,:,:)
    complex*16,allocatable::sootbx(:,:,:,:),sootby(:,:,:,:),sootbz(:,:,:,:)
    complex*16 :: vec1(3), vec2(3), cp(3), cp2(3)
    complex*16 :: spin(2, 2, 3), cplx
    real*8 :: dist(nsiti, nsiti, 3), radius
    integer :: i, j, k, a, b, c, d, si, sj, nso, asito, bsito, csito, dsito, sa, sb
    logical :: bool


    nso=nsiti*2
    allocate(hsootb(3, nso, nso, nso, nso), sootb(nso, nso, nso, nso), ppso(3, nsiti, nsiti),check(nso,nso))
    allocate(sootbx(nso, nso, nso, nso),sootby(nso, nso, nso, nso),sootbz(nso, nso, nso, nso))
    ! Define the complex unit
    cplx = (0.0d0, 1.0d0)

    ! Define the spin matrices
    spin = 0
    spin(1, 2, 1) = 1d0
    spin(2, 1, 1) = 1d0

    spin(1, 2, 2) = -cplx
    spin(2, 1, 2) = cplx

    spin(1, 1, 3) = 1d0
    spin(2, 2, 3) = -1d0

    ! Allocate arrays
    ppso=0d0
    do i =1, nsiti
       do j = 1, nsiti
          do k = 1, 3
             ppso(k,i,j)= cplx * hop(i,j) * (nuclei(i,k)-nuclei(j,k))
          enddo
       enddo
    enddo


    ! Initialize hsootb array
    hsootb = (0.0d0, 0.0d0)

    ! Calculate hsootb
    do a = 1, nso
       do b = 1, nso
          do c = 1, nso
             do  d = 1, nso
                asito = (a+1)/2
                bsito = (b+1)/2
                csito = (c+1)/2
                dsito = (d+1)/2
                if(mod(a,2).eq.0)then
                   sa = 2
                else
                   sa = 1
                endif

                if(mod(c,2).eq.0)then
                   sb = 2
                else
                   sb = 1
                endif

                if((sa.eq.sb).and.(bsito.eq.dsito))then
                   if(asito.ne.bsito)then
                      vec1 = 0.0d0
                      vec2 = 0.0d0
                      cp = 0.0d0
                      do k = 1, 3
                         vec1(k) = nuclei(asito, k) - nuclei(bsito, k)
                         vec2(k) = ppso(k, asito,csito)
                      end do
                      cp = cross_product(vec1, vec2)

                      radius = 0.0d0
                      do k = 1, 3
                         radius = radius + dreal(vec1(k))**2
                      end do
                      radius = (dsqrt(radius))**3
                      cp = cp / radius
                   else
                      cp=0d0
                   endif

                   if(csito.ne.bsito)then
                      vec1 = 0.0d0
                      vec2 = 0.0d0
                      cp2 = 0.0d0
                      do k = 1, 3
                         vec1(k) = ppso(k, asito, csito)
                         vec2(k) =nuclei( csito, k) - nuclei( bsito, k)
                      end do
                      cp2 = cross_product(vec1, vec2)

                      radius = 0.0d0
                      do k = 1, 3
                         radius = radius + dreal(vec2(k))**2
                      end do
                      radius = (dsqrt(radius))**3
                      cp2 = cp2 / radius
                   else
                      cp2=0d0
                   endif

                   if (mod(b, 2) == 0) then
                      si = 2
                   else
                      si = 1
                   end if

                   if (mod(d, 2) == 0) then
                      sj = 2
                   else
                      sj = 1
                   end if
                   do k = 1, 3
                      hsootb(k, a, b, c, d) = hsootb(k, a, b, c, d) + (cp(k) - cp2(k)) * spin(si, sj, k)
                   end do
                endif
             enddo
          end do
       enddo
    enddo



    ! Initialize sootb array
    sootb = (0.0d0, 0.0d0)
    sootbx = (0.0d0, 0.0d0)
    sootby = (0.0d0, 0.0d0)
    sootbz = (0.0d0, 0.0d0)

    ! Calculate sootb

    do a = 1, nso
       do b = 1, nso
          do c = 1, nso
             do d = 1, nso
                do k = 1, 3
                   ! sootb(a, b, c, d) = sootb(a, b, c, d) + 0.5d0 * (hsootb(k, a, b, c, d) + dconjg(hsootb(k, c, d, a, b)))
                   sootb(a, b, c, d) = sootb(a, b, c, d) + 0.5*hsootb(k, a, b, c, d)
                end do
                ! if(zabs(sootb(a,b,c,d)).ge.1d-4) write(1,'(<4>(I3, 2x), 2x, <2>(f10.5, 2x))') a, b, c, d, dreal(sootb(a,b,c,d)), dimag(sootb(a,b,c,d))
             end do
          end do
       end do
    end do


    do a = 1, nso
       do b = 1, nso
          do c = 1, nso
             do d = 1, nso
                sootbx(a, b, c, d) = sootbx(a, b, c, d) + 0.5d0 * (hsootb(1, a, b, c, d))
                sootby(a, b, c, d) = sootby(a, b, c, d) + 0.5d0 * (hsootb(2, a, b, c, d))
                sootbz(a, b, c, d) = sootbz(a, b, c, d) + 0.5d0 * (hsootb(3, a, b, c, d))
             end do
          end do
       end do
    end do

    allocate(conj_transpose_matrix(nso, nso, nso, nso))
    do a = 1, nso
       do b = 1, nso
          do c = 1, nso
             do d = 1, nso
                conj_transpose_matrix(a, b, c, d) = dconjg(sootb(c, d, a, b))
             end do
          end do
       enddo
    enddo
    bool = all(zabs(sootb - conj_transpose_matrix).le.1d-8)
    IF(bool)write(*,*) 'SOOTB hermitian'
    ! Initialize and calculate soo
    soo = 0.0d0
    soox = 0d0
    sooy = 0d0
    sooz = 0d0
    call bielectron(dim2, nso, pf, vecconfig, sootb, soo)
    call bielectron(dim2, nso, pf, vecconfig, sootbx, soox)
    call bielectron(dim2, nso, pf, vecconfig, sootby, sooy)
    call bielectron(dim2, nso, pf, vecconfig, sootbz, sooz)

    call check_hermitian(soo, dim2, bool)
    if (bool) write(*, *) 'SOO hermitian'

    call check_hermitian(sooy, dim2, bool)
    if ( .not.bool) write(*, *) 'SOOY Problem'

    call check_hermitian(soox, dim2, bool)
    if ( .not.bool) write(*, *) 'SOOX Problem'

    call check_hermitian(sooz, dim2, bool)
    if ( .not.bool) write(*, *) 'SOOZ Problem'
    soo = 2*soo
    soox = 2*soox
    sooy = 2*sooy
    sooz = 2*sooz
  end subroutine compute_soo



  subroutine bielectron(dim2, nso, pf, vecconfig, sootb, soo)
    implicit none
    integer, intent(in) :: dim2, nso
    integer, intent(in) :: vecconfig(dim2)
    complex*16, intent(out) :: soo(dim2, dim2)
    complex*16, intent(in) :: sootb(nso, nso, nso, nso)
    complex*16 :: phase
    complex*16, intent(in) :: pf
    integer :: n, a, b, c, d, m, temp, perm, i
    logical:: bool

    soo = 0.0d0

    !$omp parallel do default(none) private(n, a, b, c, d, m, temp, perm, phase, i) shared(dim2, nso, pf, vecconfig, sootb, soo)
    do n = 1, dim2
       do d = 0, nso - 1
          do c = 0, nso - 1
             do b = 0, nso - 1
                do a = 0, nso - 1
                   perm = 0
                   if (btest(vecconfig(n), d)) then
                      if (d /= nso - 1) then
                         do i = nso - 1, d + 1, -1
                            if (btest(vecconfig(n), i)) perm = perm + 1
                         enddo
                      endif
                      temp = ibclr(vecconfig(n), d)
                      if (btest(temp, c)) then
                         if (c /= nso - 1) then
                            do i = nso - 1, c + 1, -1
                               if (btest(temp, i)) perm = perm + 1
                            enddo
                         endif
                         temp = ibclr(temp, c)
                         if (.not. btest(temp, b)) then
                            if (b /= nso - 1) then
                               do i = nso - 1, b + 1, -1
                                  if (btest(temp, i)) perm = perm + 1
                               enddo
                            endif
                            temp = ibset(temp, b)
                            if (.not. btest(temp, a)) then
                               if (a /= nso - 1) then
                                  do i = nso - 1, a + 1, -1
                                     if (btest(temp, i)) perm = perm + 1
                                  enddo
                               endif
                               temp = ibset(temp, a)
                               m = binary_search(vecconfig, temp, 1, dim2)
                               if (m /= 0) then
                                  if (mod(perm, 2) == 0) then
                                     phase = +1
                                  else
                                     phase = -1
                                  endif
                                  soo(n, m) = soo(n, m) + pf * phase * sootb(a + 1, b + 1, c + 1, d + 1)
                               endif
                            endif
                         endif
                      endif
                   endif
                enddo
             enddo
          enddo
       enddo
    end do
    !$omp end parallel do
    ! call check_hermitian(soo, dim2, bool)
    ! if(.not.bool) write(*,*) 'Not bool'
  end subroutine bielectron

  function cross_product(v1, v2) result(result_vector)
    implicit none
    complex*16, dimension(3), intent(in) :: v1, v2
    complex*16, dimension(3) :: result_vector

    result_vector(1) = v1(2) * v2(3) - v1(3) * v2(2)
    result_vector(2) = v1(3) * v2(1) - v1(1) * v2(3)
    result_vector(3) = v1(1) * v2(2) - v1(2) * v2(1)
  end function cross_product
  subroutine assign_labels(dim2, nsiti, charge, label)
    implicit none
    integer, intent(in) :: dim2, nsiti
    real(8), intent(in) :: charge(dim2, nsiti)
    character(len=10), intent(out) :: label(dim2)
    integer :: i
    real(8), dimension(nsiti) :: state

    do i = 1, dim2
       state = charge(i, :)
       if (all(dabs(state) < 3.0e-1)) then
          label(i) = 'N'

       else if (all(dabs(state - (/1.0d0, 0.0d0, 0.0d0, -1.0d0/)) < 3.0e-1)) then
          label(i) = 'LRCT'

       else if (all(dabs(state - (/1.0d0, -1.0d0, 0.0d0, 0.0d0/)) < 3.0e-1)) then
          label(i) = 'DB1CT'

       else if (all(dabs(state - (/1.0d0, 0.0d0, -1.0d0, 0.0d0/)) < 3.0e-1)) then
          label(i) = 'DB2CT'
       else
          label(i) = 'ALTRO'
       end if
    end do
  end subroutine assign_labels

subroutine square_complex_matrix(n, A)
    implicit none
    integer, intent(in) :: n
    complex*16, intent(inout) :: A(n,n)
    complex*16 :: temp(n,n)
    integer :: i, j, k

    ! Calcola il prodotto di A per se stessa
    do i = 1, n
       do j = 1, n
          temp(i,j) = (0.0d0, 0.0d0)
          do k = 1, n
             temp(i,j) = temp(i,j) + A(i,k) * A(k,j)
          end do
       end do
    end do

    ! Copia il risultato nel parametro di input 
    do i = 1, n
       do j = 1, n
          A(i,j) = temp(i,j)
       end do
    end do

  end subroutine square_complex_matrix
 subroutine sz_on_conf(dim, nso, basis, sz)
    integer, intent(in)::dim, nso, basis(dim)
    real*8, intent(out):: sz(dim)
    integer:: i, j
    sz=0d0
    do i = 1, dim
       do j = 0, nso-1
          if(btest(basis(i),j))then
             if(mod(j,2).eq.0)then
                sz(i)=sz(i)+0.5
             else
                sz(i)=sz(i)-0.5
             end if
          end if
       end do
    enddo
  end subroutine sz_on_conf
subroutine sort_by_real(arr_i, arr_r, n)
  implicit none
  integer, intent(inout) :: arr_i(:)
  real*8,    intent(inout) :: arr_r(:)
  integer, intent(in)    :: n
  integer :: i, j, tmp_i
  real*8    :: tmp_r

  do i = 1, n-1
     do j = i+1, n
        if (arr_r(j) < arr_r(i)) then
           ! scambia real
           tmp_r    = arr_r(i)
           arr_r(i) = arr_r(j)
           arr_r(j) = tmp_r
           ! scambia integer corrispondente
           tmp_i    = arr_i(i)
           arr_i(i) = arr_i(j)
           arr_i(j) = tmp_i
        end if
     end do
  end do
end subroutine sort_by_real
end module newmodule
