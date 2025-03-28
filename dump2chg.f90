Module inp
  Implicit None
  Character (Len=3) :: itype
  Character (Len=3), Allocatable :: elent(:)
  Real :: slat(3, 3)
  Real :: scell(3, 3)

Contains
  Subroutine readinp
    Implicit None
    Integer :: ios, ntype
    Character (Len=10) :: line

! Read inp.dat
    Open (1001, File='inp.dat')
    Do While (.True.)
      Read (1001, *, Iostat=ios) line
      If (ios/=0) Then
        Write (*, '(a)') 'No ntype iterm!'
        Exit
      End If
      If (index(line,'ntype')/=0) Then
        Read (1001, *) ntype
        Write (*, '(a15, i10)') 'ntype:', ntype
        Exit
      End If
    End Do
    Allocate (elent(ntype))

    Rewind (1001)
    Do While (.True.)
      Read (1001, *, Iostat=ios) line
      If (ios/=0) Then
        Write (*, '(a)') 'No element iterm!'
        Exit
      End If
      If (index(line,'element')/=0) Then
        Read (1001, *) elent(:)
        Write (*, '(a15, 10a)') 'element: ', elent
        Exit
      End If
    End Do

    Rewind (1001)
    Do While (.True.)
      Read (1001, *, Iostat=ios) line
      If (ios/=0) Then
        Write (*, '(a)') 'No itype iterm!'
        Exit
      End If
      If (index(line,'itype')/=0) Then
        Read (1001, *) itype
        Write (*, '(a15,a)') 'itype: ', itype
        Exit
      End If
    End Do

    Rewind (1001)
    Do While (.True.)
      Read (1001, *, Iostat=ios) line
      If (ios/=0) Then
        Write (*, '(a)') 'No slat iterm!'
        Exit
      End If
      If (index(line,'slat')/=0) Then
        Read (1001, *) slat(1, :)
        Read (1001, *) slat(2, :)
        Read (1001, *) slat(3, :)
        Write (*, '(a15)') 'S lattice:'
        Write (*, '(3f10.5)') slat
        Exit
      End If
    End Do

    Rewind (1001)
    Do While (.True.)
      Read (1001, *, Iostat=ios) line

      If (ios/=0) Then
        Write (*, '(a)') 'No scell iterm!'
        Exit
      End If
      If (index(line,'scell')/=0) Then
        Read (1001, *) scell(1, :)
        Read (1001, *) scell(2, :)
        Read (1001, *) scell(3, :)
        Write (*, '(a15)') 'Scell:'
        Write (*, '(3f10.5)') scell
        Exit
      End If
    End Do
    Close (1001)
  End Subroutine readinp

  Subroutine inv_matrix(a, b)
    Implicit None
    Real, Intent (In) :: a(3, 3)
    Real, Intent (Out) :: b(3, 3)
    Real :: c(3, 6)
    Real :: identity(3, 3)
    Integer :: i, j, k
    Real :: temp

! Initialize the identity matrix
    identity = 0.0
    Do i = 1, 3
      identity(i, i) = 1.0
    End Do

! Merge A and the identity matrix into an augmented matrix
    Do i = 1, 3
      Do j = 1, 3
        c(i, j) = a(i, j)
      End Do
    End Do

    Do i = 1, 3
      Do j = 1, 3
        c(i, j+3) = identity(i, j)
      End Do
    End Do

! Gauss-jordan elimination method
    Do k = 1, 3
! Normalizes the current row
      temp = c(k, k)
      Do j = 1, 6
        c(k, j) = c(k, j)/temp
      End Do

! Eliminates the current column of the other rows
      Do i = 1, 3
        If (i/=k) Then
          temp = c(i, k)
          Do j = 1, 6
            c(i, j) = c(i, j) - temp*c(k, j)
          End Do
        End If
      End Do
    End Do

! Extract inverse matrix
    Do i = 1, 3
      Do j = 1, 3
        b(i, j) = c(i, j+3)
      End Do
    End Do
  End Subroutine inv_matrix
End Module inp

Program dump2chg
  Use inp
  Implicit None

  Real :: inv_slat(3, 3), lat(3, 3), inv_scell(3, 3), ixyz(3)
  Real :: eps = 1.0e-3
  Integer :: ii, jj, kk, mm, natms
  Integer :: nconfig, ios
  Character (Len=10) :: line
  Character (Len=3) :: chatom

!
  Call readinp
  Call inv_matrix(slat, inv_slat)
  Call inv_matrix(scell, inv_scell)
  lat = matmul(slat, inv_scell)

! Read dump.xyz
  nconfig = 0
  Open (1001, File='dump.xyz')
  Read (1001, *) natms
  Do While (.True.)
    Read (1001, *, Iostat=ios) line
    If (ios/=0) Then
      Exit
    End If
    If (index(line,'Time')/=0 .or. index(line,'pbc')/=0 .or. index(line,'Lattice')/=0) nconfig = nconfig + 1
  End Do
  Write (*, '(a15, i10)') 'nconfig:', nconfig

  Rewind (1001)
  
  Open (1002, File='pos.xyz')
  mm = 0
  Do ii = 1, nconfig
    Read (1001, *)
    Read (1001, *)
    Do jj = 1, natms
	  Read (1001, *) chatom, ixyz(:)
	  If(trim(adjustl(itype)) == trim(adjustl(chatom))) Then
        ixyz = matmul(ixyz, inv_slat)
		ixyz = matmul(ixyz, scell)
        Do kk = 1, 3
		  ixyz(kk) = Mod(ixyz(kk)+1000.0, 1.0)
		End Do
		Write(1002, '(a, 3f20.15)')chatom, ixyz(:)
		mm = mm + 1
	  End If
    End Do
  End Do
  Close(1001)
  Close(1002)
  
  Open (1002, File='unit.xyz')
  Write(1002, *) mm
  Write(1002, '(a,9f20.15,a)') 'pbc="T T T" Lattice="', lat, '" Properties=species:S:1:pos:R:3'
  Call system("cat pos.xyz>>unit.xyz")
  Close(1002)
End Program dump2chg