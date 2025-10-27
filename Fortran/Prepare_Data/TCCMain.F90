program TCCMain
  use ModuleGlobalData
  use TCCModule
  implicit none

  integer :: argc
  character(len=256) :: arg1
  integer :: i, ios, NumNodes
  real(4), allocatable :: NodeX(:), NodeY(:), NodeZ(:)
  character(len=100) :: input_filename, data_path, full_path
  integer :: unit = 10

  ! data_path = "../../Data/Input/" ! Sem executar pelo arquivo .sh
  data_path = "./Data/Input/" ! Executando pelo arquivo .sh

  argc = command_argument_count()

  if (argc < 1) then
      print *, "Uso: ./meu_programa <arquivo_dados>"
      stop 1
  end if

  call get_command_argument(1, arg1)
  input_filename = trim(arg1)

  full_path = trim(data_path) // trim(input_filename)

  ! Abrir arquivo de entrada
  open(unit=unit, file=trim(full_path), status="old", action="read", iostat=ios)
  if (ios /= 0) then
    print *, "Erro ao abrir o arquivo de entrada."
    stop
  end if

  ! Contar número de linhas
  NumNodes = 0
  do
    read(unit, *, iostat=ios)
    if (ios /= 0) exit
    NumNodes = NumNodes + 1
  end do

  rewind(unit)

  print *, 'NumNodes:', NumNodes
  if (NumNodes <= 0) stop 'Número de nós inválido.'


  ! Alocar arrays
  allocate(NodeX(NumNodes))
  allocate(NodeY(NumNodes))
  allocate(NodeZ(NumNodes))

  ! Ler os dados
  do i = 1, NumNodes
    read(unit, *, iostat=ios) NodeX(i), NodeY(i), NodeZ(i)
    if (ios /= 0) then
      print *, "Erro na leitura da linha ", i
      stop
    end if
  end do
  close(unit)

  call PrepareDataConstructTriangulationXYZ(0, NumNodes, real(NodeX, kind=4), real(NodeY, kind=4), real(NodeZ, kind=4), 0.0, input_filename)

  ! Liberar memória
  deallocate(NodeX, NodeY, NodeZ)

end program TCCMain
