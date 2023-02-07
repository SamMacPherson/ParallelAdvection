program AdvectionParallel
IMPLICIT NONE
include 'mpif.h'

! Parallel variables
integer	:: myid, ierror,numprocs,status(MPI_STATUS_SIZE)
integer :: tagforward,tagbackward,length, sendRequest, recieveRequest
! Grid variables and time stepping variables
integer 		:: Ngrid,np,si,ei,NTimeSteps
real(kind=8)	:: x0,xmax, dx, dt,TotalTime, BC0,BCEND
! Timing variables
real(kind=8)	:: t1,t2,t1Comm,t2Comm,t1SingleStep,t2SingleStep,sumTime
real(kind=8)	:: sumCommTime,sumSingleStepTime, currentTime,lastTimeStep
! Loop integers
integer::	t,k,n,i,p
! Solution arrays and solution calculation
real(kind=8)::				fstar,fstarminusone
real(kind=8),allocatable::	f(:),fplusone(:),X(:),finalResult(:),XfinalResult(:)
! Strings
character(len=10)::	Func,solver,str1,str2,str3
! Misc.
integer:: 				repeat
real(kind=8) :: 		u, CFL
integer, allocatable ::	idArray(:)


call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, NUMPROCS, ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, MYID, ierror)

tagforward=2001
tagbackward=2002


open (unit=1,file="config",status="old")
  read (1,*) TotalTime
  read (1,*) Ngrid
  read (1,*) x0
  read (1,*) xmax
  read (1,*) u
  read (1,*) solver
  read (1,*) Func
  read (1,*) BC0
  read (1,*) BCEND
  read (1,*) CFL
  read (1,*) repeat
 close(1)



! Set timings initially to zero
sumTime = 0
sumCommTime = 0
sumSingleStepTime = 0

! Can repeat the entire calculation to average timings if wished
! this is a small piece of code so is possible
do p = 1,repeat

	! Grid and time calculations
	dx = (xmax-x0)/(Ngrid-1)
	dt = dx/u*CFL !for non-diffusive advection

	! Number of points per processor
	np = Ngrid/numprocs

	NTimeSteps = TotalTime/dt

	! 0 is master
	! Everything else is worker



	! Calculate the processors individual domain
	if(myid.ne.0) then ! Worker processer
		! Find indices
		si = myid*np !start index
		ei = si+np+1 !end index
		if(myid.eq.(numprocs-1)) then
			ei = Ngrid
		end if
	
		! Allocate arrays
		n= ei-si+1 !length of array
		allocate(f(n),fplusone(n),X(n))
		
		! Generate X grid
		do i = 1,n
			X(i) = x0-dx +i*dx + (si-1)*dx
		end do

		! Enforce boundary condition
		if(myid.eq.(numprocs-1)) then
			f(n) = BCEND
			fplusone(n) = BCEND
		end if
	
	else ! Master processor
		! Array to store Id's later
		allocate(idArray(numprocs))

		! For timing
		t1 = MPI_WTIME()

		! Find indices
		si = 1
		ei = 1+np
		n = ei-si+1
		allocate(f(n),fplusone(n),X(n))

		do i = si,ei
			X(i) = x0-dx +i*dx
		end do

		f(1)=BC0
		fplusone(1) = BC0
	end if
		
	! Enforce initial condition 
	select case(Func)
		case('EXP2')
			call exp2(X,f,n)
		case ('STEP')
			call stepfunction(X,f,n)
	end select

	
	! Main time loop
	currentTime = 0
	do t=1,NTimeSteps+1
		if(t.eq.(NTimeSteps+1)) then
			dt = TotalTime - currentTime
		end if

		if(myid.eq.0) then
			t1SingleStep = MPI_WTIME()
		end if

		currentTime = currentTime+dt

		! Select flux based on solver
		select case(solver)
			case('UPWIND')
				do k = 2,n-1
					fplusone(k) = f(k) - u*dt/dx*(f(k) -f(k-1))
				end do
			case('FTCS')
				do k = 2,n-1
					fplusone(k) = f(k) - u*dt/2/dx*(f(k+1)-f(k-1))
				end do
			case('LAX-F')
				do k = 2,n-1
					fplusone(k) = 0.5*(f(k+1)+f(k-1)) - u*dt/2/dx*(f(k+1)-f(k-1))
				end do
			case('LAX-W')
				do k = 2,n-1
					fplusone(k) = f(k) - u*dt/2/dx*(f(k+1)-f(k-1)) + 0.5*u*u*dt*dt/dx/dx*(f(k+1)-2*f(k)+f(k-1)) 
				end do
			case('MAC')
				do k = 2,n-1
					fstar = f(k) - u*dt/dx*(f(k+1)-f(k))
					fstarminusone = f(k-1)- u*dt/dx*(f(k)-f(k-1))
					fplusone(k) = 0.5*( (f(k)+fstar)-u*dt/dx*(fstar-fstarminusone)     )
				end do
		end select
		
		! Update solution
		f = fplusone	
		! Enforce boundary condition
		if(myid.eq.(numprocs-1)) then
			f(n) = BCEND
			fplusone(n) = BCEND
		end if
		

		! Sending and recieving of internal boundaries
		
		if(myid.ne.0) then ! Worker processor

			! Always send to previous processor
			call MPI_SENDRECV(f(2),1,MPI_DOUBLE_PRECISION,myid-1,tagbackward, &
								f(1),1,MPI_DOUBLE_PRECISION,myid-1,tagforward, &
								MPI_COMM_WORLD,status,ierror)
			
			if(myid.ne.(numprocs-1)) then
				! All but last processor sends forward
				call MPI_SENDRECV(f(n-1),1,MPI_DOUBLE_PRECISION,myid+1,tagforward, &
									f(n),1,MPI_DOUBLE_PRECISION,myid+1,tagbackward, &
									MPI_COMM_WORLD,status,ierror)
			end if
		
		else ! Master processor

			! For timing
			t1Comm = MPI_WTIME()

			call MPI_SENDRECV(f(n-1),1,MPI_DOUBLE_PRECISION,myid+1,tagforward, &
							  f(n),1,MPI_DOUBLE_PRECISION,myid+1,tagbackward, &
							  MPI_COMM_WORLD,status,ierror)
			t2Comm = MPI_WTIME()
			t2SingleStep = MPI_WTIME()
		end if
	! Time loop end
	end do
	
	! Final clean up step
	! If worker, send data to master
	! If master, collect and output data
	if(myid.ne.0) then ! Worker processor

		!send array to master
		if(myid.eq.(numprocs-1)) then
			call MPI_SEND(f(2:n),np,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierror)
		else
			call MPI_SEND(f(2:n-1),np,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierror)	
		end if
		! Send array id to master
		call MPI_SEND(myid,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierror)	

	else ! Master processor

		allocate(finalResult(Ngrid),XfinalResult(Ngrid))
		!add master results to finalResult
		finalResult = 0
		! Add master processor to result
		finalResult(1:np) = f(1:np)
		
					
		! Add worker results to result
		! Collect id of worker processors
		do i = 1,numprocs-1
			call MPI_RECV(finalResult(i*np+1),np,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,status,ierror)	
			call MPI_RECV(idArray(i+1),1,MPI_INTEGER,i,i,MPI_COMM_WORLD,status,ierror)
		end do
		
		! Timing
		t2 = MPI_WTIME()
		! Calculate timings
		sumTime = sumTime + t2 - t1
		sumCommTime = sumCommTime + t2Comm - t1Comm
		sumSingleStepTime = sumSingleStepTime + t2SingleStep - t1SingleStep

	end if

	deallocate(f,fplusone,X)
! Averaging loop
end do	

if (myid.eq.0) then
	call WriteOut(numprocs,Ngrid,TotalTime,solver,func,dx,x0,finalResult,idarray,repeat,sumCommTime,sumSingleStepTime,sumTime)
	deallocate(finalResult,XfinalResult,idarray)
end if


call MPI_FINALIZE(ierror)
end program AdvectionParallel





subroutine stepfunction (X,f,n)
	! Applies the step function initial condition 
	implicit none
	integer, intent(in) 					:: n
	real(kind=8),dimension(n),intent(inout)	:: X,f

	integer :: i

	do i = 1,n
		if (X(i)<0)then
			f(i) = 0
		else
			f(i) = 1
		end if
	end do

end subroutine stepfunction

subroutine exp2 (X,f,n)
	! Applies the exponential initial condition
	implicit none
	integer, intent(in) 					:: n
	real(kind=8),dimension(n),intent(inout)	:: X,f

	integer :: i

	do i=1,n
		f(i) = 0.5*EXP(-X(i)*X(i))
	end do

end subroutine exp2

subroutine WriteOut(numprocs,Ngrid,TotalTime,solver,func,dx,x0,finalResult,idarray,repeat,sumCommTime,sumSingleStepTime,sumTime)
	! Handles output of processors used, and final result array, and timings
	implicit none
	integer, intent(in) 						:: numprocs, Ngrid, repeat
	character(len=10), intent(in) 				:: solver, func
	real(kind=8), intent(in) 					:: dx, TotalTime, x0,sumTime,sumCommTime,sumSingleStepTime
	real(kind=8), dimension(Ngrid), intent(in)	:: finalResult
	integer, dimension(numprocs),intent(in) 	:: idarray

	real(kind=8),dimension(Ngrid)::XfinalResult
	integer :: i
	character(len=10):: str1, str2, str3
	

	select case(numprocs)
	case(1:9)
		write(str3,'(I1)') numprocs
	case(10:99)
		write(str3,'(I2)') numprocs
	case(100:999)
		write(str3,'(I3)') numprocs
	end select

	select case(Ngrid)
	case(1:9)
		write(str1,'(I1)') Ngrid
	case(10:99)
		write(str1,'(I2)') Ngrid
	case(100:999)
		write(str1,'(I3)') Ngrid
	case(1000:9999)
		write(str1,'(I4)') Ngrid
	case(10000:99999)
		write(str1,'(I5)') Ngrid
	case(100000:999999)
		write(str1,'(I6)') Ngrid

	end select
	select case(int(TotalTime))
	case(1:9)
		write(str2, '(F3.1)') TotalTime
	case(10:99)
		write(str2, '(F4.1)') TotalTime

	end select
	open(1,file = 'PDEOutput_'//trim(Func)//'_'//trim(solver)//'_'//trim(str2)//'_'//trim(str1)//'_'//trim(str3)//'procs'//'.dat')

	do i=1,Ngrid
	XfinalResult(i) = x0-dx + i*dx
	write (1,*) XfinalResult(i),finalResult(i)

	end do
	close (1)

	open (2,file = 'Processors_Used'//'_'//trim(str3)//'procs'//'.dat')

	write(2,*) (idArray(i),i=1,numprocs)

	close(2)

	print*,'Simulation Time: ',TotalTime
	print*, 'Grid points:', Ngrid
	print*, 'Solver: ', solver
	print*,'Total run time = ', sumTime/repeat
	print*,'Communication time = ', sumCommTime/repeat
	print*, 'Single step time = ', sumSingleStepTime/repeat
	
	open(1,file = 'Times_'//trim(Func)//'_'//trim(solver)//'_'//trim(str2)//'_'//trim(str1)//'_'//trim(str3)//'procs'//'.dat')
	
	write (1,*) SumTime/real(repeat),sumCommTime/real(repeat),sumSingleStepTime/real(repeat)
		
	
	close (1)

end subroutine WriteOut





