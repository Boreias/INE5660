program TCCMain
  use ModuleGlobalData
  use TCCModule
  implicit none

  ! Variáveis de controle
  integer :: argc
  character(len=256) :: arg1, arg2
  integer :: ios, case
  character(len=100) :: input_filename, intermediary_data_path, full_path

  ! Variáveis do arquivo intermediário
  integer :: number_of_nodes, number_of_triangles, number_of_arcs
  integer :: number_of_boundary_nodes, lnew, swap_node
  logical :: must_swap
  real(4) :: min_x, max_x, min_y, max_y

  ! Vetores e matrizes
  real(4), allocatable :: node_x(:), node_y(:), node_z(:)
  integer, allocatable :: list(:), lptr(:), lend(:), near(:), next_tri(:)
  integer, allocatable :: bnodes(:)
  real(4), allocatable :: dist(:)
  integer, allocatable :: ltri(:,:)

  ! Caminhos de dados
  ! intermediary_data_path = "../../../Data/Intermediary/" ! Sem executar pelo arquivo .sh
  intermediary_data_path = "./Data/Intermediary/" ! Executando pelo arquivo .sh

  argc = command_argument_count()

  if (argc < 2) then
      print *, "Uso: ./meu_programa <case> <arquivo_dados>"
      stop 1
  end if

  call get_command_argument(1, arg1)
  read(arg1, *, iostat=ios) case
  if (ios /= 0) then
      print *, "Erro: o primeiro argumento deve ser um inteiro."
      stop 1
  end if

  call get_command_argument(2, arg2)

  input_filename = trim(arg2)

  full_path = trim(intermediary_data_path) // trim(input_filename)


  call ReadIntermediaryFile(full_path, &
      number_of_nodes, number_of_triangles, number_of_arcs, &
      number_of_boundary_nodes, lnew, must_swap, swap_node, &
      min_x, max_x, min_y, max_y, &
      node_x, node_y, node_z, &
      list, lptr, lend, near, next_tri, ltri, bnodes, dist)

  call ConstructTriangulationXYZ(case, number_of_nodes, number_of_triangles, number_of_arcs, &
       number_of_boundary_nodes, lnew, must_swap, swap_node, &
       min_x, max_x, min_y, max_y, &
       node_x, node_y, node_z, &
       list, lptr, lend, near, next_tri, ltri, bnodes, dist, 0.0)

contains

  subroutine ReadIntermediaryFile(filename, &
    NumberOfNodes, NumberOfTriangles, NumberOfArcs, NumberOfBoundaryNodes, &
    Lnew, MustSwap, SwapNode, MinX, MaxX, MinY, MaxY, &
    NodeX, NodeY, NodeZ, List, Lptr, Lend, Near, NextTri, Ltri, Bnodes, Dist)
  implicit none

  ! ---------- interface ----------
  character(len=*), intent(in) :: filename

  integer, intent(out) :: NumberOfNodes, NumberOfTriangles, NumberOfArcs, NumberOfBoundaryNodes
  integer, intent(out) :: Lnew, SwapNode
  logical, intent(out) :: MustSwap
  real(4), intent(out) :: MinX, MaxX, MinY, MaxY

  real(4), allocatable, intent(out) :: NodeX(:), NodeY(:), NodeZ(:), Dist(:)
  integer, allocatable, intent(out) :: List(:), Lptr(:), Lend(:), Near(:), NextTri(:), Bnodes(:)
  integer, allocatable, intent(out) :: Ltri(:,:)  ! matrix 9 x NumberOfTriangles

  ! ---------- locais ----------
  character(len=512) :: line
  character(len=512) :: s
  integer :: ios, unit
  integer :: p, iidx, jidx, val
  real(8) :: rv8
  real(4) :: rv4
  integer :: pos1, pos2, start1, start2
  integer :: r, c

  ! inicializações seguras
  NumberOfNodes = 0
  NumberOfTriangles = 0
  NumberOfArcs = 0
  NumberOfBoundaryNodes = 0
  Lnew = 0
  SwapNode = 0
  MustSwap = .false.
  MinX = 0.0
  MaxX = 0.0
  MinY = 0.0
  MaxY = 0.0

  if (allocated(NodeX)) deallocate(NodeX)
  if (allocated(NodeY)) deallocate(NodeY)
  if (allocated(NodeZ)) deallocate(NodeZ)
  if (allocated(List)) deallocate(List)
  if (allocated(Lptr)) deallocate(Lptr)
  if (allocated(Lend)) deallocate(Lend)
  if (allocated(Near)) deallocate(Near)
  if (allocated(NextTri)) deallocate(NextTri)
  if (allocated(Bnodes)) deallocate(Bnodes)
  if (allocated(Dist)) deallocate(Dist)
  if (allocated(Ltri)) deallocate(Ltri)

  ! -----------------------------
  ! 1ª passagem: ler escalares (para alocação)
  ! -----------------------------
  open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
  if (ios /= 0) then
    print *, 'ReadIntermediaryFile: erro ao abrir ', trim(filename)
    stop
  end if

  do
    read(unit, '(A)', iostat=ios) line
    if (ios /= 0) exit
    s = adjustl(trim(line))
    if (len_trim(s) == 0) cycle

    if (index(s, 'number_of_nodes') > 0) then
      call read_int_after_label(s, 'number_of_nodes', NumberOfNodes)
    end if
    if (index(s, 'number_of_triangles') > 0) then
      call read_int_after_label(s, 'number_of_triangles', NumberOfTriangles)
    end if
    if (index(s, 'number_of_arcs') > 0) then
      call read_int_after_label(s, 'number_of_arcs', NumberOfArcs)
    end if
    if (index(s, 'number_of_boundary_nodes') > 0) then
      call read_int_after_label(s, 'number_of_boundary_nodes', NumberOfBoundaryNodes)
    end if
    if (index(s, 'Lnew:') > 0 .or. index(s, 'Lnew') > 0) then
      call read_int_after_label(s, 'Lnew', Lnew)
    end if
    if (index(s, 'swap_node:') > 0 .or. index(s, 'swap_node') > 0) then
      call read_int_after_label(s, 'swap_node', SwapNode)
    end if
    if (index(s, 'must_swap') > 0) then
      call read_logical_after_label(s, 'must_swap', MustSwap)
    end if
    if (index(s,'min_x') > 0 .and. index(s,'max_x') > 0) then
      call read_two_reals_after_labels(s, 'min_x', 'max_x', MinX, MaxX)
    else
      if (index(s,'min_x') > 0) call read_real_after_label(s, 'min_x', MinX)
      if (index(s,'max_x') > 0) call read_real_after_label(s, 'max_x', MaxX)
    end if
    if (index(s,'min_y') > 0 .and. index(s,'max_y') > 0) then
      call read_two_reals_after_labels(s, 'min_y', 'max_y', MinY, MaxY)
    else
      if (index(s,'min_y') > 0) call read_real_after_label(s, 'min_y', MinY)
      if (index(s,'max_y') > 0) call read_real_after_label(s, 'max_y', MaxY)
    end if
  end do

  close(unit)

  ! Agora aloca os vetores principais sabendo os tamanhos
  if (NumberOfNodes > 0) then
    allocate(NodeX(NumberOfNodes))
    allocate(NodeY(NumberOfNodes))
    allocate(NodeZ(NumberOfNodes))
    allocate(Lend(NumberOfNodes))
    allocate(Near(NumberOfNodes))
    allocate(NextTri(NumberOfNodes))
    allocate(Bnodes(NumberOfNodes))
    allocate(List(6*NumberOfNodes-12))
    allocate(Lptr(6*NumberOfNodes-12))
    allocate(Ltri(9, 2*NumberOfNodes))
    NodeX = 0.0; NodeY = 0.0; NodeZ = 0.0; Lend = 0; Near = 0; NextTri = 0
    Bnodes = 0; List = 0; Lptr = 0; Ltri = 0
  end if

  ! -----------------------------
  ! 2ª passagem: preencher vetores / matriz
  ! -----------------------------
  open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
  if (ios /= 0) then
    print *, 'ReadIntermediaryFile: erro ao reabrir ', trim(filename)
    stop
  end if

  do
    read(unit, '(A)', iostat=ios) line
    if (ios /= 0) exit
    s = adjustl(trim(line))
    if (len_trim(s) == 0) cycle

    ! --- Node coordinates XT(i) YT(i) ZT(i) ---
    if (index(s, 'XT(') > 0) then
      call parse_indexed_real_single(s, iidx, rv4)
      if (iidx >= 1 .and. iidx <= NumberOfNodes) NodeX(iidx) = rv4
      cycle
    end if
    if (index(s, 'YT(') > 0) then
      call parse_indexed_real_single(s, iidx, rv4)
      if (iidx >= 1 .and. iidx <= NumberOfNodes) NodeY(iidx) = rv4
      cycle
    end if
    if (index(s, 'ZT(') > 0) then
      call parse_indexed_real_single(s, iidx, rv4)
      if (iidx >= 1 .and. iidx <= NumberOfNodes) NodeZ(iidx) = rv4
      cycle
    end if

    ! --- List, Lptr, Lend, Near, NextTri, BNodes (inteiros indexados) ---
    if (index(s, 'List(') > 0) then
      call parse_indexed_int(s, iidx, val)
      if (iidx >= 1 .and. iidx <= size(List)) List(iidx) = val
      cycle
    end if
    if (index(s, 'Lptr(') > 0) then
      call parse_indexed_int(s, iidx, val)
      if (iidx >= 1 .and. iidx <= size(Lptr)) Lptr(iidx) = val
      cycle
    end if
    if (index(s, 'Lend(') > 0) then
      call parse_indexed_int(s, iidx, val)
      if (iidx >= 1 .and. iidx <= size(Lend)) Lend(iidx) = val
      cycle
    end if
    if (index(s, 'Near(') > 0) then
      call parse_indexed_int(s, iidx, val)
      if (iidx >= 1 .and. iidx <= size(Near)) Near(iidx) = val
      cycle
    end if
    if (index(s, 'NextTri(') > 0) then
      call parse_indexed_int(s, iidx, val)
      if (iidx >= 1 .and. iidx <= size(NextTri)) NextTri(iidx) = val
      cycle
    end if
    if (index(s, 'BNodes(') > 0) then
      call parse_indexed_int(s, iidx, val)
      if (iidx >= 1 .and. iidx <= size(Bnodes)) Bnodes(iidx) = val
      cycle
    end if

    ! --- Dist(i) (real) ---
    if (index(s, 'Dist(') > 0) then
      call parse_indexed_real_single(s, iidx, rv4)
      if (.not. allocated(Dist)) then
        ! aloca dinamicamente caso necessário (tamanho aproximado NumberOfArcs ou NumberOfNodes)
        allocate(Dist(max(1, NumberOfArcs, NumberOfNodes)))
        Dist = 0.0
      end if
      if (iidx >= 1 .and. iidx <= size(Dist)) Dist(iidx) = rv4
      cycle
    end if

    ! --- Ltri(i,j) = val  (Matriz 9 x NumberOfTriangles) ---
    if (index(s, 'Ltri(') > 0) then
      call parse_ltri_entry(s, r, c, val)
      if (r >= 1 .and. c >= 1 .and. allocated(Ltri)) then
        if (r <= 9 .and. c <= NumberOfTriangles) then
          Ltri(r, c) = val
        end if
      end if
      cycle
    end if

  end do

  close(unit)
end subroutine ReadIntermediaryFile

  subroutine read_int_after_label(line_in, label, outint)
    character(len=*), intent(in) :: line_in, label
    integer, intent(out) :: outint
    integer :: p, iosl, colon
    character(len=256) :: rest, tmp

    outint = 0
    p = index(line_in, label)
    if (p > 0) then
      rest = adjustl(line_in(p + len_trim(label):))
      ! Remove o ":" se houver
      colon = index(rest, ":")
      if (colon > 0) then
        tmp = adjustl(rest(colon + 1:))
      else
        tmp = rest
      end if
      read(tmp, *, iostat=iosl) outint
      if (iosl /= 0) outint = 0
    end if
  end subroutine

  subroutine read_real_after_label(line_in, label, outreal)
    character(len=*), intent(in) :: line_in, label
    real(4), intent(out) :: outreal
    integer :: p, iosl, colon
    character(len=256) :: rest, tmp
    outreal = 0.0
    p = index(line_in, label)
    if (p > 0) then
      rest = adjustl(line_in(p + len_trim(label):))
      colon = index(rest, ":")
      if (colon > 0) then
        tmp = adjustl(rest(colon + 1:))
      else
        tmp = rest
      end if
      read(tmp, *, iostat=iosl) outreal
      if (iosl /= 0) outreal = 0.0
    end if
  end subroutine

  subroutine read_two_reals_after_labels(line_in, label1, label2, out1, out2)
    implicit none
    character(len=*), intent(in) :: line_in, label1, label2
    real(4), intent(out) :: out1, out2
    integer :: p1, p2
    character(len=512) :: part1, part2

    out1 = 0.0
    out2 = 0.0

    p1 = index(line_in, label1)
    p2 = index(line_in, label2)

    if (p1 > 0 .and. p2 > p1) then
      ! extrai substring entre o fim do label1 e o inicio do label2
      part1 = adjustl(line_in(p1 + len_trim(label1) : p2 - 1))
      call normalize_after_label(part1, part1)
      call tokenize_first_real(part1, out1)

      ! extrai substring após label2 até o fim da linha
      part2 = adjustl(line_in(p2 + len_trim(label2) : ))
      call normalize_after_label(part2, part2)
      call tokenize_first_real(part2, out2)
    end if
  end subroutine read_two_reals_after_labels

  subroutine normalize_after_label(str_in, str_out)
    implicit none
    character(len=*), intent(in) :: str_in
    character(len=*), intent(out) :: str_out
    integer :: pos_colon, pos_eq, pos_first

    pos_colon = index(str_in, ':')
    pos_eq    = index(str_in, '=')

    ! se houver ':' ou '=', e estiverem no começo (ou antes de tokens), pular até depois do símbolo
    if (pos_colon > 0 .and. (pos_eq == 0 .or. pos_colon < pos_eq)) then
      ! pega substring após o primeiro ':'
      if (pos_colon < len_trim(str_in)) then
        str_out = adjustl(str_in(pos_colon + 1 : len_trim(str_in)))
      else
        str_out = ''
      end if
    else if (pos_eq > 0) then
      if (pos_eq < len_trim(str_in)) then
        str_out = adjustl(str_in(pos_eq + 1 : len_trim(str_in)))
      else
        str_out = ''
      end if
    else
      ! nada a pular, só ajustar
      str_out = adjustl(str_in)
    end if
  end subroutine normalize_after_label


  subroutine tokenize_first_real(str_in, outreal)
    implicit none
    character(len=*), intent(in) :: str_in
    real(4), intent(out) :: outreal
    character(len=512) :: tmp
    integer :: i, L, iosl, start_pos, end_pos
    character(len=1) :: ch

    tmp = adjustl(str_in)
    L = len_trim(tmp)
    outreal = 0.0

    if (L == 0) return

    ! Encontrar início do token numérico (pule espaços e caracteres não numéricos até achar +,-,digit ou '.')
    start_pos = 1
    do while (start_pos <= L)
      ch = tmp(start_pos:start_pos)
      if (index('0123456789+-.' , ch) > 0) exit
      start_pos = start_pos + 1
    end do
    if (start_pos > L) return

    ! Encontrar fim do token numérico: primeiro caractere que não pertence a [0-9+-.Ee]
    end_pos = start_pos
    do while (end_pos <= L)
      ch = tmp(end_pos:end_pos)
      if (index('0123456789+-.Ee' , ch) == 0 .and. ch /= ' ') exit
      end_pos = end_pos + 1
    end do
    end_pos = end_pos - 1
    if (end_pos < start_pos) return

    ! Substring contendo o token candidato
    tmp = tmp(start_pos:end_pos)

    ! Fazer leitura interna
    read(tmp,*,iostat=iosl) outreal
    if (iosl /= 0) outreal = 0.0
  end subroutine tokenize_first_real

  subroutine read_logical_after_label(line_in, label, outlog)
    character(len=*), intent(in) :: line_in, label
    logical, intent(out) :: outlog
    integer :: p, colon
    character(len=256) :: rest, tmp
    outlog = .false.
    p = index(line_in, label)
    if (p > 0) then
      rest = adjustl(line_in(p + len_trim(label):))
      colon = index(rest, ":")
      if (colon > 0) then
        tmp = adjustl(rest(colon + 1:))
      else
        tmp = rest
      end if
      tmp = adjustl(tmp)
      if (index(tmp, 't') == 1 .or. index(tmp, '.true.') > 0) outlog = .true.
    end if
  end subroutine

  pure function to_lower(str) result(res)
    implicit none
    character(len=*), intent(in) :: str
    character(len=len(str)) :: res
    integer :: i, c

    res = str
    do i = 1, len(str)
      c = iachar(str(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) then
        res(i:i) = achar(c + 32)
      end if
    end do
  end function to_lower

  subroutine parse_indexed_real_single(line_in, idx_out, val_out)
    character(len=*), intent(in) :: line_in
    integer, intent(out) :: idx_out
    real(4), intent(out) :: val_out
    integer :: p1, p2, eq, iosl
    character(len=256) :: sidx, sval
    idx_out = 0; val_out = 0.0
    p1 = index(line_in,'(')
    p2 = index(line_in,')')
    eq = index(line_in,'=')
    if (p1>0 .and. p2>p1 .and. eq>p2) then
      sidx = adjustl(line_in(p1+1:p2-1))
      sval = adjustl(line_in(eq+1:))
      read(sidx,*,iostat=iosl) idx_out
      if (iosl == 0) then
        call tokenize_first_real(sval, val_out)
      end if
    end if
  end subroutine

  subroutine parse_indexed_int(line_in, idx_out, val_out)
    implicit none
    character(len=*), intent(in) :: line_in
    integer, intent(out) :: idx_out, val_out
    integer :: p1, p2, eq, iosl
    character(len=256) :: sidx, sval

    idx_out = 0
    val_out = 0

    p1 = index(line_in, '(')
    p2 = index(line_in, ')')
    eq = index(line_in, '=')

    if (p1 > 0 .and. p2 > p1 .and. eq > p2) then
      sidx = adjustl(line_in(p1+1 : p2-1))
      sval = adjustl(line_in(eq+1 :))
      ! Lê o índice
      read(sidx, *, iostat=iosl) idx_out
      if (iosl == 0) then
        ! Lê o valor numérico (em vez de usar int() em string)
        read(sval, *, iostat=iosl) val_out
        if (iosl /= 0) val_out = 0
      end if
    end if
  end subroutine parse_indexed_int

  subroutine parse_ltri_entry(line_in, i_out, j_out, v_out)
    implicit none
    character(len=*), intent(in) :: line_in
    integer, intent(out) :: i_out, j_out, v_out
    integer :: p1, pc, p2, eq, iosl
    character(len=256) :: s1, s2, s3

    ! Inicializa saídas
    i_out = 0
    j_out = 0
    v_out = 0

    ! Localiza posições dos parênteses, vírgula e sinal de igual
    p1 = index(line_in, '(')
    pc = index(line_in, ',')
    p2 = index(line_in, ')')
    eq = index(line_in, '=')

    ! Verifica se todos os delimitadores foram encontrados na ordem esperada
    if (p1 > 0 .and. pc > p1 .and. p2 > pc .and. eq > p2) then
      ! Extrai substrings correspondentes aos números
      s1 = adjustl(line_in(p1+1 : pc-1))
      s2 = adjustl(line_in(pc+1 : p2-1))
      s3 = adjustl(line_in(eq+1 : len_trim(line_in)))

      ! Lê os três números
      read(s1, *, iostat=iosl) i_out
      if (iosl /= 0) i_out = 0

      read(s2, *, iostat=iosl) j_out
      if (iosl /= 0) j_out = 0

      read(s3, *, iostat=iosl) v_out
      if (iosl /= 0) v_out = 0
    end if
  end subroutine parse_ltri_entry

end program TCCMain
