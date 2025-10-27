!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Triangulation
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module to perform Delaunay Triangulation 
!
!------------------------------------------------------------------------------
!
!This program is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License 
!version 2, as published by the Free Software Foundation.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software
!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!
!------------------------------------------------------------------------------

Module TCCModule

    use ModuleGlobalData

    implicit none 

    private

    !Constructor
    public  :: ConstructTriangulationXYZ
    private ::      AllocateInstance
    private ::      ConstructDirectedEdges
    private ::          AlreadyInList
    private ::          AddEdge
    private ::      ConstructNodes
    private ::      ConstructTriangles

    !Auxilary
    private ::      circum
    private ::      nbcnt


    !Parameter                              
    integer, parameter                              :: pnRow        = 9  !See Robert Renka Code

    integer, parameter                              :: Pit          = 1
    integer, parameter                              :: Unflooded    = 2
    integer, parameter                              :: Flooded      = 3
    integer, parameter                              :: CurrentLake  = 4

    !Types---------------------------------------------------------------------

    type T_DirectedEdge
        integer                                     :: EdgeID
        integer                                     :: StartNode
        integer                                     :: EndNode
        real                                        :: VoronoiX
        real                                        :: VoronoiY
        real                                        :: Length
        real                                        :: Slope
        logical                                     :: Boundary
        type (T_DirectedEdge), pointer              :: CounterClockEdge
    end type T_DirectedEdge

    type T_Node
        real                                        :: X
        real                                        :: Y
        real                                        :: Z
        integer                                     :: nNeighbor = 0
        logical                                     :: Boundary
        integer                                     :: State
        real                                        :: VoronoiArea
        type (T_DirectedEdge), pointer              :: FirstEdge
        type (T_DirectedEdge), pointer              :: FlowExit
    end type T_Node

    type T_Triangles
        real                                        :: CenterX, CenterY
        real                                        :: Radius, Area
        real                                        :: AspectRatio
    end type T_Triangles

    type T_Reach
        type (T_DirectedEdge), pointer              :: Edge
        integer                                     :: nStrahler            = null_int
        real                                        :: DrainageArea         = null_real
        type (T_Reach), pointer                     :: Next
    end type T_Reach

    type       T_Triangulation
        !Instance ID
        integer                                     :: InstanceID
        real,    dimension(:), pointer              :: XT, YT, ZT
        integer, dimension(:), pointer              :: List, Lptr
        integer, dimension(:), pointer              :: Lend
        integer                                     :: Lnew
        integer, dimension(:), pointer              :: Near, NextTri
        integer, dimension(:, :), pointer           :: Ltri
        integer, dimension(:), pointer              :: BNodes
        real, dimension(:), pointer                 :: Dist
        integer                                     :: NumberOfNodes
        integer                                     :: NumberOfTriangles
        integer                                     :: NumberOfArcs
        integer                                     :: NumberOfBoundaryNodes
        logical                                     :: MustSwap = .false.
        integer                                     :: SwapNode
        real                                        :: MinX  = -null_real
        real                                        :: MaxX  =  null_real
        real                                        :: MinY  = -null_real
        real                                        :: MaxY  =  null_real
        type (T_Node), dimension(:), pointer        :: Nodes
        type (T_DirectedEdge), dimension(:), pointer :: DirectedEdges
        type (T_Triangles), dimension(:), pointer   :: Triangles
        logical                                     :: HaveHeightValues     = .false.
        type (T_Reach), pointer                     :: FirstReach
        integer                                     :: nReaches             = 0
        type (T_Triangulation), pointer             :: Next
    end type T_Triangulation

    !Global Variables
    type (T_Triangulation), pointer                 :: FirstTriangulation
    type (T_Triangulation), pointer                 :: Me

    !--------------------------------------------------------------------------

    contains

    subroutine ConstructTriangulationXYZ(TriangulationID,                       &
                                         NumberOfNodes, NumberOfTriangles, NumberOfArcs, &
                                          NumberOfBoundaryNodes, Lnew, MustSwap, SwapNode, &
                                          MinX, MaxX, MinY, MaxY, &
                                          NodeX, NodeY, NodeZ, &
                                          List, Lptr, Lend, Near, NextTri, Ltri, BNodes, Dist,    &
                                         Tolerance, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TriangulationID
        real                                        :: Tolerance
        integer, optional, intent(OUT)              :: STAT     

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: STAT_

        integer :: NumberOfNodes, NumberOfTriangles, NumberOfArcs
        integer :: NumberOfBoundaryNodes, Lnew, SwapNode
        logical :: MustSwap
        real(4) :: MinX, MaxX, MinY, MaxY

        ! Vetores e matrizes
        real(4), allocatable :: NodeX(:), NodeY(:), NodeZ(:)
        integer, allocatable :: List(:), Lptr(:), Lend(:), Near(:), NextTri(:)
        integer, allocatable :: BNodes(:)
        real(4), allocatable :: Dist(:)
        integer, allocatable :: Ltri(:,:)

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call AllocateInstance

        Me%NumberOfNodes = NumberOfNodes

        !Allocates variables of ObjTriangulation
        allocate(Me%XT(1:NumberOfNodes), STAT = STAT_CALL)
        allocate(Me%YT(1:NumberOfNodes), STAT = STAT_CALL)
        allocate(Me%ZT(1:NumberOfNodes), STAT = STAT_CALL)
        allocate(Me%List(1:6*NumberOfNodes-12), STAT = STAT_CALL)
        allocate(Me%Lptr(1:6*NumberOfNodes-12), STAT = STAT_CALL)
        allocate(Me%Lend(1:NumberOfNodes), STAT = STAT_CALL)
        allocate(Me%Near(1:NumberOfNodes), STAT = STAT_CALL)
        allocate(Me%NextTri(1:NumberOfNodes), STAT = STAT_CALL)
        allocate(Me%Dist(1:NumberOfNodes), STAT = STAT_CALL)
        allocate(Me%BNodes(1:NumberOfNodes), STAT = STAT_CALL)
        allocate(Me%Ltri(pnRow, 2*NumberOfNodes), STAT = STAT_CALL)

        allocate(NodeX(1:NumberOfNodes), STAT = STAT_CALL)
        allocate(NodeY(1:NumberOfNodes), STAT = STAT_CALL)
        allocate(NodeZ(1:NumberOfNodes), STAT = STAT_CALL)
        allocate(List(1:6*NumberOfNodes-12), STAT = STAT_CALL)
        allocate(Lptr(1:6*NumberOfNodes-12), STAT = STAT_CALL)
        allocate(Lend(1:NumberOfNodes), STAT = STAT_CALL)
        allocate(Near(1:NumberOfNodes), STAT = STAT_CALL)
        allocate(NextTri(1:NumberOfNodes), STAT = STAT_CALL)
        allocate(Dist(1:NumberOfNodes), STAT = STAT_CALL)
        allocate(BNodes(1:NumberOfNodes), STAT = STAT_CALL)
        allocate(Ltri(pnRow, 2*NumberOfNodes), STAT = STAT_CALL)

        Me%InstanceID = TriangulationID
        Me%XT = NodeX
        Me%YT = NodeY
        Me%ZT = NodeZ
        Me%List = List
        Me%Lptr = Lptr
        Me%Lend = Lend
        Me%Near = Near
        Me%NextTri = NextTri
        Me%Dist = Dist
        Me%BNodes = BNodes
        Me%Ltri = Ltri
        Me%NumberOfTriangles = NumberOfTriangles
        Me%NumberOfArcs = NumberOfArcs
        Me%NumberOfBoundaryNodes = NumberOfBoundaryNodes
        Me%Lnew = Lnew
        Me%MustSwap = MustSwap
        Me%SwapNode = SwapNode
        Me%MinX = MinX
        Me%MaxX = MaxX
        Me%MinY = MinY
        Me%MaxY = MaxY

        call ConstructNodes         (NodeX, NodeY, NodeZ)

        call ConstructTriangles

        call ConstructDirectedEdges

        STAT_ = SUCCESS_

        if (present(STAT)) STAT = STAT_

    end subroutine ConstructTriangulationXYZ

!     !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Triangulation), pointer                :: NewTriangulation
        type (T_Triangulation), pointer                :: PreviousTriangulation


        !Allocates new instance
        allocate (NewTriangulation)
        nullify  (NewTriangulation%Next)
        nullify  (NewTriangulation%FirstReach)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstTriangulation)) then
            FirstTriangulation    => NewTriangulation
            Me                    => NewTriangulation
        else
            PreviousTriangulation => FirstTriangulation
            Me                    => FirstTriangulation%Next
            do while (associated(Me))
                PreviousTriangulation  => Me
                Me                     => Me%Next
            enddo
            Me                         => NewTriangulation
            PreviousTriangulation%Next => NewTriangulation
        endif

        Me%InstanceID = RegisterNewInstance (mTRIANGULATION_)

    end subroutine AllocateInstance
    
    !--------------------------------------------------------------------------

    subroutine ConstructDirectedEdges
    
        !Arguments-------------------------------------------------------------


        !Local-----------------------------------------------------------------
        integer                                     :: v1, v2, v3
        integer                                     :: n1, n2, n3
        integer                                     :: e1, e2, e3
        integer                                     :: EdgeSoFar, iT, iE, iE2, iN
        integer                                     :: lpl, lp, nd, k, i
        integer, dimension(1000)                    :: nabor
        integer                                     :: EndNodeCounterClock
        integer                                     :: StartNode, EndNode

        !Allocates Variable                        
        allocate(Me%DirectedEdges(2*Me%NumberOfArcs))

        !Stores Directed Edge Information
        EdgeSoFar = 0
        do iT = 1, Me%NumberOfTriangles

            !Vertex Pointers
            v1 = Me%LTRI(1, iT)
            v2 = Me%LTRI(2, iT)
            v3 = Me%LTRI(3, iT)
            
            !Neighbor Pointers
            n1 = Me%LTRI(4, iT)
            n2 = Me%LTRI(5, iT)
            n3 = Me%LTRI(6, iT)

            !Edge Pointers    
            e1 = Me%LTRI(7, iT)
            e2 = Me%LTRI(8, iT)
            e3 = Me%LTRI(9, iT)

            !Test all three edges of the triangle if they are already in the Directed 
            !Edge list. If not add them

            !Vertex 1
            if (.not. AlreadyInList(EdgeSoFar, v2, v3)) then
                call AddEdge (EdgeSoFar, e1, v2, v3, n1, iT)
            endif

            !Vertex 2
            if (.not. AlreadyInList(EdgeSoFar, v3, v1)) then
                call AddEdge (EdgeSoFar, e2, v3, v1, n2, iT)
            endif

            !Vertex 3
            if (.not. AlreadyInList(EdgeSoFar, v1, v2)) then
                call AddEdge (EdgeSoFar, e3, v1, v2, n3, iT)
            endif

        enddo

        !Updates information about counter clockwise edge
        do iE = 1, 2 * Me%NumberOfArcs

            !Find the first CounterClockwise Neighbor of the Starting Node
            StartNode = Me%DirectedEdges(iE)%StartNode
            EndNode   = Me%DirectedEdges(iE)%EndNode

            !Robert Renka Stuff
            lpl       = Me%Lend(StartNode)
            lp        = lpl
            k         = 0

            do
                k        = k + 1
                lp       = Me%lptr(lp)
                nd       = Me%list(lp)
                nabor(k) = nd
                if (lp == lpl ) then
                    exit
                end if
            end do
            
            !NODE is a boundary node.
            if (nd <= 0) then
                nabor(k) = -nd
            end if
            
            !Find the node index of counterclock edge the list of Boundary Nodes
            do i = 1, k
                if (nabor(i) == EndNode) then
                    if (i == k) then
                        EndNodeCounterClock = nabor(1)
                    else
                        EndNodeCounterClock = nabor(i+1)
                    endif
                    exit
                endif
            enddo

            !Points to the CounterClock Edge
            do iE2 = 1, 2 * Me%NumberOfArcs
                if (Me%DirectedEdges(iE2)%StartNode   == StartNode) then
                    if (Me%DirectedEdges(iE2)%EndNode == EndNodeCounterClock) then
                        Me%DirectedEdges(iE)%CounterClockEdge =>           &
                            Me%DirectedEdges(iE2)
                        exit
                    endif
                endif
            enddo

        enddo

        !Updates PointerMe%s to first edges
        do iN = 1, Me%NumberOfNodes
doEdge:     do iE = 1, 2 * Me%NumberOfArcs
                if (iN == Me%DirectedEdges(iE)%StartNode) then
                    Me%Nodes(iN)%FirstEdge => Me%DirectedEdges(iE)
                    exit doEdge
                endif
            enddo doEdge
        enddo


    end subroutine ConstructDirectedEdges

    ! --------------------------------------------------------------------------

    logical function AlreadyInList (EdgeSoFar, StartNode, EndNode)

        !Arguments-------------------------------------------------------------
        integer                                     :: EdgeSoFar, StartNode, EndNode

        !Local-----------------------------------------------------------------
        integer                                     :: iE

        AlreadyInList = .false.
        do iE = EdgeSoFar, 1, -1 
            if (Me%DirectedEdges(iE)%StartNode   == StartNode) then
                if (Me%DirectedEdges(iE)%EndNode == EndNode  ) then
                    AlreadyInList = .true.
                    return
                endif
            endif
        enddo

    end function AlreadyInList

    !--------------------------------------------------------------------------

    subroutine AddEdge (EdgeSoFar, EdgeID, StartNode, EndNode, Nb, iT)

        !Arguments-------------------------------------------------------------
        integer                                     :: EdgeSoFar, EdgeID
        integer                                     :: StartNode, EndNode
        integer                                     :: Nb, iT

        !Local-----------------------------------------------------------------
        real                                        :: dx, dy, dz, Length, Slope
        logical                                     :: BoundaryEdge

        !Length / Slope between nodes
        dx = Me%Nodes(StartNode)%X - Me%Nodes(EndNode)%X
        dy = Me%Nodes(StartNode)%Y - Me%Nodes(EndNode)%Y
        dz = Me%Nodes(StartNode)%Z - Me%Nodes(EndNode)%Z

        Length = sqrt(dx**2. + dy**2.)
        Slope  = dz / Length

        !Verifies if Edge is a Boundary edge
        if (Me%Nodes(StartNode)%Boundary .and.                             &
            Me%Nodes(EndNode)%Boundary .and. Nb == 0) then
            BoundaryEdge = .true.
        else
            BoundaryEdge = .false.
        endif


        !Increase Number of Edges and stores the Edge
        EdgeSoFar = EdgeSoFar + 1

        Me%DirectedEdges(EdgeSoFar)%EdgeID    = EdgeID
        Me%DirectedEdges(EdgeSoFar)%StartNode = StartNode
        Me%DirectedEdges(EdgeSoFar)%EndNode   = EndNode
        Me%DirectedEdges(EdgeSoFar)%Length    = Length
        Me%DirectedEdges(EdgeSoFar)%Slope     = Slope
        Me%DirectedEdges(EdgeSoFar)%Boundary  = BoundaryEdge
        !Voronoi Center of the right hand side 
        if (Nb /= 0) then
            Me%DirectedEdges(EdgeSoFar)%VoronoiX =                         &
                                                    Me%Triangles(Nb)%CenterX
            Me%DirectedEdges(EdgeSoFar)%VoronoiY =                         &
                                                    Me%Triangles(Nb)%CenterY
        else
            Me%DirectedEdges(EdgeSoFar)%VoronoiX  = null_real
            Me%DirectedEdges(EdgeSoFar)%VoronoiY  = null_real
        endif

        !Stores the complementary Edge
        EdgeSoFar = EdgeSoFar + 1
        Me%DirectedEdges(EdgeSoFar)%EdgeID    = EdgeID
        Me%DirectedEdges(EdgeSoFar)%StartNode = EndNode
        Me%DirectedEdges(EdgeSoFar)%EndNode   = StartNode
        Me%DirectedEdges(EdgeSoFar)%Length    = Length
        Me%DirectedEdges(EdgeSoFar)%Slope     = -1. * Slope
        Me%DirectedEdges(EdgeSoFar)%Boundary  = BoundaryEdge

        !Voronoi Center of the right hand side 
        Me%DirectedEdges(EdgeSoFar)%VoronoiX =                         &
                                                Me%Triangles(iT)%CenterX
        Me%DirectedEdges(EdgeSoFar)%VoronoiY =                         &
                                                Me%Triangles(iT)%CenterY

    
    end subroutine AddEdge

    !--------------------------------------------------------------------------

    subroutine ConstructNodes (NodeX, NodeY, NodeZ)

        !Arguments-------------------------------------------------------------
        real, dimension(Me%NumberOfNodes)               :: NodeX, NodeY, NodeZ
            
        !Local-----------------------------------------------------------------
        integer                                         :: iN, iB

        !Allocates Nodes
        allocate(Me%Nodes(Me%NumberOfNodes))

        do iN = 1, Me%NumberOfNodes
            nullify (Me%Nodes(iN)%FlowExit )
            nullify (Me%Nodes(iN)%FirstEdge)
        enddo
        

        Me%HaveHeightValues   = .true.
        Me%Nodes(:)%nNeighbor = 0
        
        do iN = 1, Me%NumberOfNodes
            Me%Nodes(iN)%X         = NodeX(iN)
            Me%Nodes(iN)%Y         = NodeY(iN)
            Me%Nodes(iN)%Z         = NodeZ(iN)
            Me%Nodes(iN)%nNeighbor = nbcnt(Me%lend(iN), Me%lptr)
            Me%Nodes(iN)%Boundary  = .false.
            do iB = 1, Me%NumberOfBoundaryNodes
                if (Me%BNodes(iB) == iN) then
                   Me%Nodes(iN)%Boundary = .true.
                   exit
                endif
            enddo
            
            !Checks Min and Maximum
            if (Me%Nodes(iN)%X < Me%MinX) Me%MinX = Me%Nodes(iN)%X
            if (Me%Nodes(iN)%X > Me%MaxX) Me%MaxX = Me%Nodes(iN)%X
            if (Me%Nodes(iN)%Y < Me%MinY) Me%MinY = Me%Nodes(iN)%Y
            if (Me%Nodes(iN)%Y > Me%MaxY) Me%MaxY = Me%Nodes(iN)%Y

        enddo

    end subroutine ConstructNodes
    
    !--------------------------------------------------------------------------

    subroutine ConstructTriangles

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                         :: iT

        !Allocates Triangles
        allocate(Me%Triangles(Me%NumberOfTriangles))

        !Calculates the Area, Aspect Ratio, CircumCenter of each Triangle
        do iT = 1, Me%NumberOfTriangles
            call Circum (x1 = Me%Nodes(Me%Ltri(1, iT))%X,     &
                         y1 = Me%Nodes(Me%Ltri(1, iT))%Y,     &
                         x2 = Me%Nodes(Me%Ltri(2, iT))%X,     &
                         y2 = Me%Nodes(Me%Ltri(2, iT))%Y,     &
                         x3 = Me%Nodes(Me%Ltri(3, iT))%X,     &
                         y3 = Me%Nodes(Me%Ltri(3, iT))%Y,     &
                         ratio = .true.,                                                  &
                         xc = Me%Triangles(iT)%CenterX,                     &
                         yc = Me%Triangles(iT)%CenterY,                     &
                         cr = Me%Triangles(iT)%Radius,                      &
                         sa = Me%Triangles(iT)%Area,                        &
                         ar = Me%Triangles(iT)%AspectRatio)
        enddo
        


    end subroutine ConstructTriangles

    !--------------------------------------------------------------------------

    subroutine circum ( x1, y1, x2, y2, x3, y3, ratio, xc, yc, cr, sa, ar )
    !
    !*******************************************************************************
    !
    !! CIRCUM determines the circumcenter (and more) of a triangle.
    !
    !
    !  Discussion:
    !
    !    Given three vertices defining a triangle, this subrou-
    !    tine returns the circumcenter, circumradius, signed
    !    triangle area, and, optionally, the aspect ratio of the
    !    triangle.
    !
    !  Author:
    !
    !    Robert Renka,
    !    Department of Computer Science,
    !    University of North Texas,
    !    renka@cs.unt.edu
    !
    !  Parameters:
    !
    !    Input, real X1, Y1, X2, Y2, X3, Y3, the coordinates of the vertices.
    !
    !    Input, logical RATIO, is TRUE if and only if the aspect ratio is 
    !    to be computed.
    !
    !    Output, real XC, YC, coordinates of the circumcenter (center of the 
    !    circle defined by the three points) unless SA = 0, in which XC and YC
    !    are not altered.
    !
    !    Output, real CR, the circumradius (radius of the circle defined by
    !    the three points) unless SA = 0 (infinite radius), in which case 
    !    CR is not altered.
    !
    !    Output, real SA, the signed triangle area with positive value if
    !    and only if the vertices are specified in counterclockwise order:  
    !    (X3,Y3) is strictly to the left of the directed line from (X1,Y1)
    !    toward (X2,Y2).
    !
    !    Output, real AR, the aspect ratio r/CR, where r is the radius of the
    !    inscribed circle, unless RATIO = FALSE, in which case AR is not 
    !    altered.  AR is in the range 0 to 0.5, with value 0 iff SA = 0 and
    !    value 0.5 iff the vertices define an equilateral triangle.
    !
      implicit none
    !
      real ar
      real cr
      real ds(3)
      real fx
      real fy
      logical ratio
      real sa
      real u(3)
      real v(3)
      real x1
      real x2
      real x3
      real xc
      real y1
      real y2
      real y3
      real yc
    !
    !  Set U(K) and V(K) to the x and y components, respectively,
    !  of the directed edge opposite vertex K.
    !
      u(1) = x3 - x2
      u(2) = x1 - x3
      u(3) = x2 - x1
      v(1) = y3 - y2
      v(2) = y1 - y3
      v(3) = y2 - y1
    !
    !  Set SA to the signed triangle area.
    !
      sa = ( u(1) * v(2) - u(2) * v(1) ) / 2.0E+00

      if ( sa == 0.0E+00 ) then
        if ( ratio ) then
          ar = 0.0E+00
        end if
        return
      end if
    !
    !  Set DS(K) to the squared distance from the origin to vertex K.
    !
      ds(1) = x1 * x1 + y1 * y1
      ds(2) = x2 * x2 + y2 * y2
      ds(3) = x3 * x3 + y3 * y3
    !
    !  Compute factors of XC and YC.
    !
      fx = - dot_product ( ds(1:3), v(1:3) )
      fy =   dot_product ( ds(1:3), u(1:3) )

      xc = fx / ( 4.0E+00 * sa )
      yc = fy / ( 4.0E+00 * sa )
      cr = sqrt ( ( xc - x1 )**2 + ( yc - y1 )**2 )

      if ( .not. ratio ) then
        return
      end if
    !
    !  Compute the squared edge lengths and aspect ratio.
    !
      ds(1:3) = u(1:3)**2 + v(1:3)**2

      ar = 2.0E+00 * abs ( sa ) / &
          ( ( sqrt ( ds(1) ) + sqrt ( ds(2) ) + sqrt ( ds(3) ) ) * cr )

      return
    end subroutine circum

    integer function nbcnt ( lpl, lptr )
    !
    !*******************************************************************************
    !
    !! NBCNT returns the number of neighbors of a node.
    !
    !
    !  Discussion:
    !
    !    This function returns the number of neighbors of a node
    !    N0 in a triangulation created by Subroutine TRMESH (or
    !    TRMSHR).
    !
    !  Author:
    !
    !    Robert Renka,
    !    Department of Computer Science,
    !    University of North Texas,
    !    renka@cs.unt.edu
    !
    !  Parameters:
    !
    !    Input, integer LPL, the LIST pointer to the last neighbor of N0.
    !    LPL = LEND(N0).
    !
    !    Input, integer LPTR(*), pointers associated with LIST.
    !
    !    Output, integer NBCNT, the  number of neighbors of N0.
    !
      implicit none
    !
      integer k
      integer lp
      integer lpl
      integer lptr(*)
    !
      lp = lpl
      k  = 1

      do
        lp = lptr(lp)
        if ( lp == lpl ) then
          exit
        end if
        k = k + 1

      end do

      nbcnt = k

      return
    end function nbcnt

end Module TCCModule

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior T�cnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
