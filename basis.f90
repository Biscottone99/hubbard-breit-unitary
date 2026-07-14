program basis
  implicit none
  integer :: nsiti, i, max, min, count, config, nf, a, b, n, ne, nso
  logical :: bool
  
  ! Variabili per la gestione in memoria prima della scrittura
  integer, allocatable :: basis_array(:)
  integer, allocatable :: bit_array(:)
  integer, allocatable :: sz2_array(:)
  character(len=8), allocatable :: char_arrays(:)
  character(len=8) :: temp_char
  integer :: dim_max

  ! Variabili per l'ordinamento
  integer :: j, temp_basis, temp_bit, temp_sz2

  ! Parametri del sistema
  nsiti = 4     ! siti di disposizione elettroni
  nso = nsiti * 2
  ne = 4
  
  ! Calcolo dei limiti per la rappresentazione decimale (bit representation)
  max = 0
  do i = nso-1, nso-ne, -1
     max = max + 2**i
  enddo

  min = 0
  do i = 0, ne-1
     min = min + 2**i
  enddo

  ! Dimensione massima teorica (nso su ne) per allocare gli array
  ! Per nso=8, ne=4, dim_max = 70. 
  ! Usiamo un'allocazione abbondante basata su max-min
  dim_max = max - min + 1
  allocate(basis_array(dim_max))
  allocate(bit_array(dim_max))
  allocate(sz2_array(dim_max))
  allocate(char_arrays(dim_max))

  !!!! 1. Generazione delle configurazioni in memoria
  nf = 0
  do n = min, max
     count = 0
     config = 0
     a = 0
     b = 0
     temp_char = '00000000' ! Inizializza la stringa di 8 caratteri
     
     do i = 0, nso-1
        bool = btest(n, i)
        if (bool) then
           temp_char(i+1:i+1) = '1'
           count = count + 1
           ! Convenzione: spin SU sui siti pari, spin GIU sui siti dispari
           if (mod(i, 2) == 0) then
              a = a + 1
           else
              b = b + 1
           endif
        endif
     enddo
     
     ! Se abbiamo il numero corretto di elettroni, salviamo in memoria
     if (count == ne) then
        nf = nf + 1
        
        ! Calcola la configurazione intera
        config = 0
        do i = 0, nso-1
           if (temp_char(i+1:i+1) == '1') then
              config = config + 2**i
           endif
        enddo
        
        ! Salva i dati nei vettori
        basis_array(nf) = config
        bit_array(nf) = n
        char_arrays(nf) = temp_char
        sz2_array(nf) = a - b ! 2*Sz
     endif
  enddo

  !!!! 2. Ordinamento (Selection Sort) basato sullo spin (Sz)
  do i = 1, nf - 1
     do j = i + 1, nf
        if (sz2_array(j) < sz2_array(i)) then
           ! Scambia 2*Sz
           temp_sz2 = sz2_array(i)
           sz2_array(i) = sz2_array(j)
           sz2_array(j) = temp_sz2
           
           ! Scambia configurazione intera (basis)
           temp_basis = basis_array(i)
           basis_array(i) = basis_array(j)
           basis_array(j) = temp_basis
           
           ! Scambia numero intero di partenza (n)
           temp_bit = bit_array(i)
           bit_array(i) = bit_array(j)
           bit_array(j) = temp_bit
           
           ! Scambia array di caratteri
           temp_char = char_arrays(i)
           char_arrays(i) = char_arrays(j)
           char_arrays(j) = temp_char
        endif
     enddo
  enddo

  !!!! 3. Scrittura definitiva sui file
  open(unit=1, file='basis.dat', status='replace')
  open(unit=2, file='configurations.dat', status='replace')
  open(unit=3, file='dim2.dat', status='replace')

  write(3,*) nf

  do i = 1, nf
     ! Scrive la base ordinata
     write(1,*) basis_array(i)
     
     ! Scrive configurations.dat come lo volevi tu: n, stringa separata da spazi, spin reale
     write(2,*) bit_array(i), (char_arrays(i)(j:j), j=1, nso), dble(sz2_array(i)) * 0.5d0
  enddo

  close(1)
  close(2)
  close(3)

  write(*,*) 'Numero di funzioni di base totale:', nf
  write(*,*) 'I file basis.dat e configurations.dat sono stati generati e ordinati per Sz.'

  ! Pulizia memoria
  deallocate(basis_array)
  deallocate(bit_array)
  deallocate(sz2_array)
  deallocate(char_arrays)

end program basis
