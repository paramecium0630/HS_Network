program HS_Network
    implicit none
    integer :: i, j, k, ii, jj
    integer :: N, nstep, npast, nedge
    integer :: LWORK, LWORK2, LDA, LDVL, LDVR
    integer :: nreal, iseed, idir, ntau, na, nb, NN
    real(8) :: tt, tpast, protocol
    real(8) :: dum, time1, time2
    real(8) :: DD, dt, prob, W_mean, W_sig, r0, w0, period, tau, bb
    real(8), parameter :: pi = 4d0*atan(1d0)

    ! network variables
    integer, allocatable :: G0(:,:), indeg(:), outdeg(:)
    real(8), allocatable :: a(:), r(:), noise(:,:), G(:,:), Q(:,:)
    real(8), allocatable :: inweightdeg(:), outweightdeg(:)

    ! dynamic variables
    real(8), allocatable :: f1(:), x1(:), gau(:)!, xt(:,:), ft(:,:)
    real(8), allocatable :: f_avg(:,:), x_avg(:,:), K_avg(:,:,:)
    real(8), allocatable :: f_real(:,:), x_real(:,:)
    real(8), allocatable :: dxi(:), dqi(:), dsi(:), dwi(:)
    ! real(8), allocatable :: protocol(:)

    ! energetics
    real(8), allocatable :: heat_sim(:), entropy_sim(:), work_sim(:), energy_sim(:)

    call read_variables_from_file()
	N = NN
	LWORK = 67*N
	LWORK2 = N
	LDA = N; LDVL = N; LDVR = N

    ISEED = -ISEED
    dum = RAN3(ISEED) ! initialize random number
    tpast = 5d0
    NSTEP = INT(period/dt)
    NPAST = INT(tpast/dt) 

    if (N <= 0) stop "N must be > 0."
    allocate(G0(N,N), indeg(N), outdeg(N))
    allocate(a(N), r(N), noise(N, nreal), G(N,N), Q(N,N))
    allocate(inweightdeg(N), outweightdeg(N))

    allocate(f1(N), x1(N), gau(N))
    allocate(dxi(N), dqi(N), dsi(N), dwi(N))
    allocate(f_avg(N, -NPAST:NSTEP), x_avg(N, -NPAST:NSTEP), K_avg(N, N, -NPAST:NSTEP))
    allocate(f_real(N, -NPAST:NSTEP), x_real(N, -NPAST:NSTEP))

    OPEN(UNIT=11, FILE='config/G0_connect.csv', STATUS='REPLACE')
	OPEN(UNIT=12, FILE='config/G_connect.csv', STATUS='REPLACE')	     
	OPEN(UNIT=13, FILE='config/noise.csv', STATUS='REPLACE')
	OPEN(UNIT=14, FILE='config/parameter.csv', STATUS='REPLACE')

    OPEN(UNIT=15, FILE='output/x.dat', STATUS='REPLACE')
    OPEN(UNIT=16, FILE='output/f.dat', STATUS='REPLACE')
    OPEN(UNIT=17, FILE='output/protocol.dat', STATUS='REPLACE')
    OPEN(UNIT=18, FILE='output/K.dat', STATUS='REPLACE')

    G = 0d0; G0 = 0
	!random graph
	do i = 1, N-1; do j = i+1, N
	   G0(i,j)=INT(1.d0+prob-ran3(iseed))
	   if(idir.eq.0) G0(j,i)=G0(i,j)
	   if(idir.eq.1) G0(j,i)=INT(1.d0+prob-ran3(iseed))	      
	enddo; enddo	
	DO i=1,N-1; DO j=i+1,N	   
	   G(i,j)=G0(i,j)*GAUDEV(w_mean, w_sig, iseed)
	   if(idir.eq.0) G(j,i)=G(i,j)
	   if(idir.eq.1) G(j,i)=G0(j,i)*GAUDEV(w_mean, w_sig, iseed)
	ENDDO; ENDDO

!dec$if(.true.)
	G0 = 0
	open(unit=999, file="input/adjacency_list.txt", status="old", action="read")	
	do
            read(999, *, end=100) ii, jj, nedge
            G0(ii,jj) = nedge
			G(ii,jj) = nedge*w0
	end do	
100 	continue
	close(999)
!dec$endif

    do i = 1, n; do j = 1, n
        if (G0(i,j) == 0) cycle
        write(11,*) i, j, G0(i,j)
        write(12,*) i, j, G(i,j)
    enddo; enddo

    nedge = 0
	do i = 1, N; do j = 1, N
	   nedge = nedge + G0(i,j)
	enddo; enddo
	print*, 'connection probability: ', 1.*nedge/(N*(N-1))
        
    indeg = 0; outdeg = 0
    inweightdeg = 0d0; outweightdeg = 0d0

    do i = 1, N
        indeg(i) = sum(G0(i,:))
        outdeg(i) = sum(G0(:,i))
        inweightdeg(i) = sum(G(i,:))
        outweightdeg(i) = sum(G(:,i))
    enddo

    a = 0; r = 0; noise = 0
    do i = 1, N
        a(i) = i*1d0
        r(i) = r0*(ran3(iseed)+0.5d0)
        noise(i,i) = dd
    enddo

    do i = 1, N
        write(14,*) i, r(i)
        write(13,*) i, noise(i,i)
    enddo

    write(6,*) repeat("=", 40)
    write(6,*) 'number of node = ', N
	write(6,*) 'directed = ', idir
	write(6,*) 'steady tstart, tend = ', 0d0, period
	write(6,*) 'number of time steps = ', nstep
    write(6,*) repeat("=", 40)

    call cpu_time(time1)

    ! initialize dynamic variables
    !xt = 0; ft = 0
    x_avg = 0; f_avg = 0; K_avg = 0
    dxi = 0; dqi = 0; dsi = 0; dwi = 0
    do k = 1, nreal
        if (mod(k, int(nreal*1./10)) == 1) call progress(k, nreal)
        protocol = 0; x1 = 0; f1 = 0; gau = 0
        do ii = -npast, NSTEP
            tt = ii*dt
            do i = 1, N
            gau(i) = sqrt(noise(i,i))*GAUDEV(0d0, 1.d0, iseed)
            enddo
        
            if (tt < 0) then
                protocol = 0
            elseif (tt >= 0 .and. tt < tau) then
                protocol = tt/tau*bb
            else
                protocol = bb
            end if

            if (k == 1) then
                if (mod(ii, 100) .eq. 1) then                    
                    write(17,*) tt, protocol
                endif
            endif

            x_real(:, ii) = x1
            call steady(N, r, a, G, protocol, f1, x1, gau)
            f_real(:, ii) = f1
        enddo
        do jj = -npast, nstep
            x_avg(:, jj) = x_avg(:, jj) + x_real(:, jj)
	   	    f_avg(:, jj) = f_avg(:, jj) + f_real(:, jj)	
        enddo
        do i = 1, N; do j = 1, N
            do jj = -npast, nstep
                K_avg(i,j,jj) = K_avg(i,j,jj) + x_real(i, jj)*x_real(j, jj)
            enddo
        enddo; enddo
    enddo
    call cpu_time(time2)

    write(*,*) ""
    write(*,*) 'CPU time for simulation: ', time2 - time1, ' seconds.'

    deallocate(x_real); deallocate(f_real)

    x_avg = x_avg/nreal; f_avg = f_avg/nreal; K_avg = K_avg/nreal

    do i = 1, n; do j = 1, n
        K_avg(i,j,:) = K_avg(i,j,:) - x_avg(i,:)*x_avg(j,:)
    enddo; enddo

    do i = -npast, nstep
        if (mod(i, 10) .eq. 1) then  
            write(15,*) i*dt, x_avg(1,i)
            write(16,*) i*dt, f_avg(1,i)
            write(18,*) i*dt, K_avg(1,1,i)
        endif
    enddo

contains

subroutine progress(i, nprog)
	implicit double precision(a-h,o-z)
	integer, intent(in) :: i, nprog
		
	write(*,*) 
	write(*, '(a17)', advance='no') ' Progress : 10.0%'
	
	val = i*1./nprog+0.1
	
	write(*,'(a, a, f5.1, a2, a1)', advance="no") char(13), ' progress: ', val*100, '%'		
	
end subroutine progress

subroutine read_variables_from_file()
    implicit none

    open(unit=10, file='input/input.dat', status='old')
    read(10, *) dd, dt, iseed ! noise, time step, seed
    read(10, *) idir, prob, ntau ! directed, connection probability, number of tau
    read(10, *) period, nreal, tau, bb ! nonsteady period, number of realizations, tau, bb
    read(10, *) W_mean, W_sig, r0, w0 ! mean and sig of weighted random graph, mean r, mean w
    read(10, *) na, nb ! na, nb
    read(10, *) NN ! number of nodes in dynamic network
    close(10)
end subroutine read_variables_from_file

subroutine steady(n, r, a, G, protocol, f1, x1, gau)
    implicit none
    integer, intent(in) :: n
    integer :: j, jn
    real(8), intent(in) :: protocol
    real(8), intent(in) :: G(n,n), r(n), a(n), gau(n)
    real(8), intent(inout) :: x1(n)
    real(8), intent(out) :: f1(n)
    real(8) :: Gt(n,n), rt(n), at(n), Gi(n), sumg

    Gt = G
    rt = r
    at = a

    rt(na) = rt(na) + protocol

    Gi = 0d0
    do jn = 1, n
        sumg = 0d0
        do j = 1, n
            sumg = sumg + Gt(jn,j)*(x1(j)-x1(jn)) !(x2-x1)
        enddo
        Gi(jn) = sumg
    enddo

    f1 = -rt*(x1-at) + Gi
    ! f1 = -rt*x1*(x1-1) + Gi
    x1 = x1 + f1 * dt + gau * sqrt(dt) ! x(t_{i+1})
end subroutine steady

real(8) function gaudev(rmu, sig, iseed)
    implicit none
    real(8), intent(in) :: rmu, sig
    integer, intent(inout) :: iseed
    real(8) :: tpi, rand1, rand2

    tpi = 6.28315307d0
    rand1 = ran3(iseed)
    do while (rand1 == 0d0)
        rand1 = ran3(iseed)
    enddo

    rand2 = ran3(iseed)
    do while (rand2 == 0d0)
        rand2 = ran3(iseed)
    enddo

    gaudev = rmu + sig*sqrt(-2.d0*log(rand1))*cos(tpi*rand2)
end function gaudev

real(8) function ran3(idum)
    implicit none
    integer, intent(inout) :: idum
    integer, parameter :: MBIG = 1000000000, MSEED = 161803398, MZ = 0
    real(8), parameter :: FAC = .999999d-9
    integer :: ma(55), iff, inext, inextp, mj, mk, i, ii, k
    save iff, inext, inextp, ma
    data iff /0/

    if (idum.lt.0 .or. iff.eq.0) then
        iff = 1
        mj = MSEED - iabs(idum)
        mj = mod(mj, MBIG)
        ma(55) = mj
        mk = 1

        do i = 1, 54
            ii = mod(21*i, 55)
            ma(ii) = mk
            mk = mj - mk
            if (mk.lt.MZ) mk = mk + MBIG
            mj = ma(ii)
        enddo

        do k = 1, 4
            do i = 1, 55
                ma(i) = ma(i) - ma(1 + mod(i + 30, 55))
                if (ma(i).lt.MZ) ma(i) = ma(i) + MBIG
            enddo
        enddo

        inext = 0
        inextp = 31
        idum = 1
    endif

    inext = inext + 1
    if (inext.eq.56) inext = 1
    inextp = inextp + 1
    if (inextp.eq.56) inextp = 1
    mj = ma(inext) - ma(inextp)
    if (mj.lt.MZ) mj = mj + MBIG
    ma(inext) = mj
    ran3 = mj * FAC
end function ran3

end program HS_Network
