program HS_Network
    implicit none
    integer :: i, j, k, ii, jj
    integer :: N, nstep, npast, NSTEP_PROTOCOL, nedge, nskip_edge
    integer :: LWORK, LWORK2, LDA, LDVL, LDVR
    integer :: nreal, iseed, idir, ntau, na, nb, NN
    real(8) :: xKx1, xKx2, HS_Y, HS_Y_avg, EXP_HS_Y_avg, Delta_phi
    real(8) :: phi_start, phi_end, Delta_phi_avg, entropy_ex, entropy_hk
    real(8) :: heat_total, heat_total_avg, entropy_total, entropy_total_avg
    real(8) :: tt, tpast, protocol, protocol_next    
    real(8) :: dum, time1, time2, time1_HS, time2_HS
    real(8) :: DD, dt, prob, W_mean, W_sig, r0, w0, period, tau, bb
    real(8), parameter :: pi = 4d0*atan(1d0)

    ! network variables
    integer, allocatable :: G0(:,:), indeg(:), outdeg(:)
    real(8), allocatable :: a(:), r(:), noise(:,:), G(:,:), Q(:,:)
    real(8), allocatable :: inweightdeg(:), outweightdeg(:)

    ! dynamic variables
    real(8), allocatable :: f1(:), x1(:), noise_increment(:), fluct(:), fluct_prev(:), x_prev(:), f2(:), dq_step(:)!, xt(:,:), ft(:,:)
    real(8), allocatable :: f_avg(:,:), x_avg(:,:), K_avg(:,:,:)
    real(8), allocatable :: f_real(:,:), x_real(:,:)
    real(8), allocatable :: dxi(:), dqi(:), dsi(:), dwi(:)
    ! real(8), allocatable :: protocol(:)

    ! schur variables
    real(8), allocatable :: QschurT(:,:), QschurVS(:,:), QschurWR(:), QschurWI(:), QschurWORK(:)
    real(8), allocatable :: QschurTMP(:,:), QschurC(:,:), QschurSol(:,:), ans(:,:)
    logical, allocatable :: QschurBWORK(:)

    ! energetics
    real(8), allocatable :: heat_sim(:), entropy_sim(:), work_sim(:), energy_sim(:)

    real(8), allocatable :: KoT(:,:), Ko_inv_protocol(:,:,:), logdet_K(:)

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
    NSTEP_PROTOCOL = INT(tau/dt)

    if (N <= 0) stop "N must be > 0."
    allocate(G0(N,N), indeg(N), outdeg(N))
    allocate(a(N), r(N), noise(N,N), G(N,N), Q(N,N))
    allocate(inweightdeg(N), outweightdeg(N))

    allocate(f1(N), x1(N), noise_increment(N), fluct(N), fluct_prev(N), x_prev(N), f2(N), dq_step(N))
    allocate(dxi(N), dqi(N), dsi(N), dwi(N))
    allocate(f_avg(N, -NPAST:NSTEP), x_avg(N, -NPAST:NSTEP), K_avg(N, N, -NPAST:NSTEP))
    allocate(f_real(N, -NPAST:NSTEP), x_real(N, -NPAST:NSTEP))

    allocate(QschurT(N,N), QschurVS(N,N), QschurWR(N), QschurWI(N), QschurWORK(LWORK))
    allocate(QschurBWORK(LWORK2))
    allocate(QschurTMP(N,N), QschurC(N,N), QschurSol(N,N), KoT(N,N), ans(N,N))
    allocate(Ko_inv_protocol(N,N,0:NSTEP_PROTOCOL), logdet_K(0:NSTEP_PROTOCOL))

    allocate(heat_sim(N), entropy_sim(N), work_sim(N), energy_sim(N))

    OPEN(UNIT=11, FILE='config/G0_connect.csv', STATUS='REPLACE')
	OPEN(UNIT=12, FILE='config/G_connect.csv', STATUS='REPLACE')	     
	OPEN(UNIT=13, FILE='config/noise.csv', STATUS='REPLACE')
	OPEN(UNIT=14, FILE='config/parameter.csv', STATUS='REPLACE')

    OPEN(UNIT=15, FILE='output/x.dat', STATUS='REPLACE')
    OPEN(UNIT=16, FILE='output/f.dat', STATUS='REPLACE')
    OPEN(UNIT=17, FILE='output/protocol.dat', STATUS='REPLACE')
    OPEN(UNIT=18, FILE='output/K.dat', STATUS='REPLACE')
    OPEN(UNIT=19, FILE='output/Ko.dat', STATUS='REPLACE')
    OPEN(UNIT=20, FILE='output/solution.dat', STATUS='REPLACE')
    OPEN(UNIT=21, FILE='output/heat.dat', STATUS='REPLACE')
    OPEN(UNIT=22, FILE='output/entropy_env.dat', STATUS='REPLACE')

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

!dec$if(.false.)
	G0 = 0
    nskip_edge = 0
	open(unit=999, file="input/adjacency_list.txt", status="old", action="read")	
	do
            read(999, *, end=100) ii, jj, nedge
            if (ii < 1 .or. ii > N .or. jj < 1 .or. jj > N) then
                nskip_edge = nskip_edge + 1
                cycle
            endif
            G0(ii,jj) = nedge
			G(ii,jj) = nedge*w0
	end do	
100 	continue
	close(999)
    if (nskip_edge > 0) write(6,*) 'skipped adjacency edges outside 1..N = ', nskip_edge
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
        ! noise(i,i) = dd
        noise(i,i) = dd*(ran3(iseed)+0.5d0)
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
    heat_sim = 0d0; entropy_sim = 0d0; work_sim = 0d0; energy_sim = 0d0
    do k = 1, nreal
        if (mod(k, int(nreal*1./10)) == 1) call progress(k, nreal)
        protocol = 0; protocol_next = 0; x1 = 0; f1 = 0; noise_increment = 0

        do ii = -npast, NSTEP
            tt = ii*dt
            do i = 1, N
                noise_increment(i) = sqrt(noise(i,i)*dt)*GAUDEV(0d0, 1.d0, iseed)
            enddo
            call protocol_at_time(tt, tau, bb, protocol)
            call protocol_at_time(tt + dt, tau, bb, protocol_next)

            if (k == 1) then
                if (mod(ii, 100) .eq. 1) then                    
                    write(17,*) tt, protocol
                endif
            endif

            x_real(:, ii) = x1
            call advance_hs_state_heun(N, r, a, G, protocol, protocol_next, noise_increment, f1, x1)
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

    x_avg = x_avg/nreal; f_avg = f_avg/nreal; K_avg = K_avg/nreal

    do i = 1, n; do j = 1, n
        K_avg(i,j,:) = K_avg(i,j,:) - x_avg(i,:)*x_avg(j,:)
    enddo; enddo

    do ii = 0, NSTEP_PROTOCOL
        tt = ii*dt
        call protocol_at_time(tt, tau, bb, protocol)
        call build_jacobian_from_protocol(N, r, G, protocol, Q)
        call schur_factorize_real_matrix(N, Q, QschurT, QschurVS, QschurWR, QschurWI, QschurWORK, QschurBWORK)
        call solve_qx_xqt_from_schur_real_ws(N, QschurT, QschurVS, -noise, QschurTMP, QschurC, QschurSol)
        KoT = QschurSol
        call matrix_logdet_spd(KoT, N, logdet_K(ii))
        call inverse(N, KoT, Ko_inv_protocol(:,:,ii))
    enddo

    call cpu_time(time1_HS)
    ! HATANO-SASA
    xKx1 = 0d0; xKx2 = 0d0
    HS_Y_avg = 0d0; EXP_HS_Y_avg = 0d0; Delta_phi_avg = 0d0
    heat_total_avg = 0d0; entropy_total_avg = 0d0
    do k = 1, nreal
        HS_Y = 0d0; phi_start = 0d0; phi_end = 0d0; Delta_phi = 0d0
        heat_total = 0d0; entropy_total = 0d0
        protocol = 0; protocol_next = 0; x1 = 0; f1 = 0; noise_increment = 0
        dqi = 0d0; dsi = 0d0

        do ii = -npast, NSTEP_PROTOCOL
            if (ii == 0) then
                fluct = x1 - x_avg(:, 0)
                xKx1 = DOT_PRODUCT(fluct, MATMUL(Ko_inv_protocol(:,:,0), fluct))
                phi_start = 0.5d0*logdet_K(0) + 0.5d0*xKx1
            endif

            if (ii >= 1) then
                fluct = x1 - x_avg(:, ii) ! fluctuation at time step ii
                fluct_prev = x1 - x_avg(:, ii-1)
                xkx1 = DOT_PRODUCT(fluct, MATMUL(Ko_inv_protocol(:,:,ii), fluct))
                xkx2 = DOT_PRODUCT(fluct_prev, MATMUL(Ko_inv_protocol(:,:,ii-1), fluct_prev))
                HS_Y = HS_Y + 0.5d0*(logdet_K(ii) - logdet_K(ii-1)) + 0.5d0*(xkx1 - xkx2)
            endif

            if (ii == NSTEP_PROTOCOL) then
                fluct = x1 - x_avg(:, NSTEP_PROTOCOL)
                xKx1 = DOT_PRODUCT(fluct, MATMUL(Ko_inv_protocol(:,:,NSTEP_PROTOCOL), fluct))
                phi_end = 0.5d0*logdet_K(NSTEP_PROTOCOL) + 0.5d0*xKx1
                Delta_phi = phi_end - phi_start
                exit ! no need to advance state at last protocol step
            endif

            tt = ii*dt
            do i = 1, N
                noise_increment(i) = sqrt(noise(i,i)*dt)*GAUDEV(0d0, 1.d0, iseed)
            enddo
            call protocol_at_time(tt, tau, bb, protocol)
            call protocol_at_time(tt + dt, tau, bb, protocol_next)
            x_prev = x1
            call advance_hs_state_heun(N, r, a, G, protocol, protocol_next, noise_increment, f1, x1)
            if (ii >= 0 .and. ii < NSTEP_PROTOCOL) then
                call calculate_hs_force(N, r, a, G, protocol_next, x1, f2)
                dq_step = -0.5d0*(f1 + f2)*(x1 - x_prev)
                dqi = dqi + dq_step
                do i = 1, N
                    dsi(i) = dsi(i) + 2d0/noise(i,i)*dq_step(i)
                enddo
            endif
        enddo
        HS_Y_avg = HS_Y_avg + HS_Y
        EXP_HS_Y_avg = EXP_HS_Y_avg + exp(-HS_Y)
        Delta_phi_avg = Delta_phi_avg + Delta_phi
        heat_sim = heat_sim + dqi
        entropy_sim = entropy_sim + dsi
        heat_total = sum(dqi)
        entropy_total = sum(dsi)
        heat_total_avg = heat_total_avg + heat_total
        entropy_total_avg = entropy_total_avg + entropy_total
    enddo
    call cpu_time(time2_HS)
    HS_Y_avg = HS_Y_avg/nreal
    EXP_HS_Y_avg = EXP_HS_Y_avg/nreal
    Delta_phi_avg = Delta_phi_avg/nreal
    heat_sim = heat_sim/nreal
    heat_total_avg = heat_total_avg/nreal
    entropy_sim = entropy_sim/nreal
    entropy_total_avg = entropy_total_avg/nreal
    write(*,*) 'CPU time for Hatano-Sasa: ', time2_HS - time1_HS, ' seconds.'
    write(6,*) 'Hatano-Sasa Y average = ', HS_Y_avg
    write(6,*) 'Hatano-Sasa Y exponential average = ', EXP_HS_Y_avg
    write(6,*) 'Delta phi average = ', Delta_phi_avg
    write(6,*) 'Excess Entropy average = ', HS_Y_avg - Delta_phi_avg
    ! write(6,*) 'Housekeeping Entropy average = ', Delta_phi_avg - HS_Y_avg 
    write(6,*) 'Heat total average = ', heat_total_avg
    write(6,*) 'Environment entropy average = ', entropy_total_avg

    do i = 1, N
        write(21,*) i, heat_sim(i)
        write(22,*) i, entropy_sim(i)
    enddo

    call protocol_at_time(0d0, tau, bb, protocol)
    call build_jacobian_from_protocol(N, r, G, protocol, Q)
    call schur_factorize_real_matrix(N, Q, QschurT, QschurVS, QschurWR, QschurWI, QschurWORK, QschurBWORK)
    call solve_qx_xqt_from_schur_real_ws(N, QschurT, QschurVS, -noise, QschurTMP, QschurC, QschurSol)
    KoT = QschurSol

    ans = MATMUL(Q, KoT) + MATMUL(KoT, TRANSPOSE(Q)) + noise

    deallocate(x_real); deallocate(f_real)

    do i = -npast/5, nstep
        if (mod(i+npast, 10) .eq. 1) then  
            write(15,*) i*dt, x_avg(1,i), x_avg(2,i)
            write(16,*) i*dt, f_avg(1,i), f_avg(2,i)
            write(18,*) i*dt, K_avg(1,1,i), K_avg(1,2,i)            
        endif
    enddo

    do i = 1, N
        do j = i, N
            write(19,*) K_avg(i,j,0), real(KoT(i,j))
            write(20,*) ans(i,j)
        enddo
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
    read(10, *) dd, dt, iseed ! noise strength, time step, seed
    read(10, *) idir, prob, ntau ! directed, connection probability, number of tau
    read(10, *) period, nreal, tau, bb ! nonsteady period, number of realizations, tau, bb
    read(10, *) W_mean, W_sig, r0, w0 ! mean and sig of weighted random graph, mean r, mean w
    read(10, *) na, nb ! na, nb
    read(10, *) NN ! number of nodes in dynamic network
    close(10)
end subroutine read_variables_from_file

subroutine protocol_at_time(tt, tau, bb, protocol)
    implicit none
    real(8), intent(in) :: tt, tau, bb
    real(8), intent(out) :: protocol

    if (tt < 0d0) then
        protocol = 0d0
    elseif (tt < tau) then
        protocol = tt/tau*bb
    else
        protocol = bb
    endif
end subroutine protocol_at_time

subroutine calculate_hs_force(n, r, a, G, protocol, x, f)
    implicit none
    integer, intent(in) :: n
    integer :: j, jn
    real(8), intent(in) :: protocol
    real(8), intent(in) :: G(n,n), r(n), a(n), x(n)
    real(8), intent(out) :: f(n)
    real(8) :: Gt(n,n), rt(n), at(n), Gi(n), sumg

    Gt = G
    rt = r
    at = a

    rt(na) = rt(na) + protocol

    Gi = 0d0
    do jn = 1, n
        sumg = 0d0
        do j = 1, n
            sumg = sumg + Gt(jn,j)*(x(j)-x(jn)) !(x2-x1)
        enddo
        Gi(jn) = sumg
    enddo

    f = -rt*(x-at) + Gi
    ! f = -rt*x*(x-1) + Gi
end subroutine calculate_hs_force

subroutine advance_hs_state_heun(n, r, a, G, protocol, protocol_next, noise_increment, f1, x1)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: protocol, protocol_next
    real(8), intent(in) :: G(n,n), r(n), a(n), noise_increment(n)
    real(8), intent(inout) :: x1(n)
    real(8), intent(out) :: f1(n)
    real(8) :: x_predictor(n), f_predictor(n)

    call calculate_hs_force(n, r, a, G, protocol, x1, f1)
    x_predictor = x1 + f1*dt + noise_increment
    call calculate_hs_force(n, r, a, G, protocol_next, x_predictor, f_predictor)
    x1 = x1 + 0.5d0*(f1 + f_predictor)*dt + noise_increment ! x(t_{i+1})
end subroutine advance_hs_state_heun

subroutine build_jacobian_from_protocol(N, r, G, protocol, Q)
		implicit none
		integer, intent(in) :: N
		real(8), intent(in) :: r(N), G(N,N), protocol
		real(8), intent(out) :: Q(N,N)
		integer :: i
		real(8) :: inweightdegree_i, rt(N)

		! HS local dynamics:
		!   f_i(x_i) = -r_i*(x_i-a_i)
		! gives f_i'(x_i) = -r_i.
        ! The protocol shifts the relaxation rate of node na.
		! Diffusive coupling h(x_i, x_j) = x_j - x_i contributes:
		!   off-diagonal:  +G(i,j)
		!   diagonal:      -sum_j G(i,j)
        rt = r
        rt(na) = rt(na) + protocol
		Q = G
		do i = 1, N
			inweightdegree_i = sum(G(i,:))
			Q(i,i) = Q(i,i) - rt(i) - inweightdegree_i
		enddo
	end subroutine build_jacobian_from_protocol

logical function select_eig(wr, wi)
	implicit none
	real(8), intent(in) :: wr, wi
	! Callback for DGEES; no reordering needed in this code path.
	select_eig = .false.
end function select_eig

subroutine schur_factorize_real_matrix(N, Q, Tschur, VS, WR, WI, WORK, BWORK)
	implicit none
	integer, intent(in) :: N
	integer :: INFO, SDIM
	real(8), intent(in) :: Q(N,N)
	real(8), intent(out) :: Tschur(N,N), VS(N,N), WR(N), WI(N)
	real(8), intent(inout) :: WORK(*)
	logical, intent(inout) :: BWORK(N)
	! logical, external :: select_eig

	Tschur = Q
	call DGEES('V', 'N', select_eig, N, Tschur, N, SDIM, WR, WI, VS, N, WORK, 8*N, BWORK, INFO)
	if (INFO .ne. 0) then
		write(6,*) "Warning: DGEES failed in schur_factorize_real_matrix, INFO =", INFO
	endif
end subroutine schur_factorize_real_matrix

subroutine solve_qx_xqt_from_schur_real_ws(N, Tschur, VS, RHS, TMPR, C, X)
		! Solves Q*X + X*Q^T = RHS for X, where Q is given in real Schur form as
		! Q = VS*Tschur*VS^T. The continuous Lyapunov equation is the special
		! case RHS = -Dnoise.
		implicit none
		integer, intent(in) :: N
		integer :: i, j, INFO
		real(8), intent(in) :: Tschur(N,N), VS(N,N), RHS(N,N)
		real(8), intent(inout) :: TMPR(N,N), C(N,N)
		real(8), intent(out) :: X(N,N)
		real(8) :: SCALE, symv
	
		! Solve Q*X + X*Q^T = RHS using the cached real Schur form Q = VS*Tschur*VS^T.
		call dgemm('N', 'N', N, N, N, 1d0, RHS, N, VS, N, 0d0, TMPR, N)
		call dgemm('T', 'N', N, N, N, 1d0, VS, N, TMPR, N, 0d0, C, N)
		call DTRSYL('N', 'T', 1, N, N, Tschur, N, Tschur, N, C, N, SCALE, INFO)
		if (INFO .ne. 0) then
			write(6,*) "Warning: DTRSYL returned INFO in solve_qx_xqt_from_schur_real_ws =", INFO
		endif
		if (SCALE .gt. 0d0 .and. abs(SCALE-1d0) .gt. 1d-14) C = C / SCALE
		call dgemm('N', 'T', N, N, N, 1d0, C, N, VS, N, 0d0, TMPR, N)
		call dgemm('N', 'N', N, N, N, 1d0, VS, N, TMPR, N, 0d0, X, N)
	
		do i = 1, N
			do j = i+1, N
				symv = 0.5d0*(X(i,j) + X(j,i))
				X(i,j) = symv
				X(j,i) = symv
			enddo
		enddo
end subroutine solve_qx_xqt_from_schur_real_ws

subroutine matrix_trace(A, n, trace)
	  	implicit none
	  	integer, intent(in) :: n
	  	real(8), intent(in) :: A(n,n)
	  	real(8), intent(out) :: trace

  	integer :: i
	! trace(A) = sum diagonal elements.

  	trace = 0.0d0
  	do i = 1, n
	     	trace = trace + A(i,i)
	  	end do
	end subroutine matrix_trace

subroutine matrix_logdet_spd(A, n, logdet)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: A(n,n)
        real(8), intent(out) :: logdet
        integer :: i, info
        real(8) :: Acopy(n,n)

        Acopy = A
        call dpotrf('U', n, Acopy, n, info)
        if (info .ne. 0) then
            write(6,*) "Warning: DPOTRF failed in matrix_logdet_spd, INFO =", info
            logdet = 0d0
            return
        endif

        logdet = 0d0
        do i = 1, n
            logdet = logdet + 2d0*log(Acopy(i,i))
        enddo
    end subroutine matrix_logdet_spd

subroutine inverse(n, A, Ainv)
		implicit double precision(a-h, o-z)
	integer :: n, lda, info
	integer :: ipiv(n)
	real(8) :: A(n,n), Ainv(n,n)
	real(8) :: work(n)
	! Utility inverse through LU (DGETRF/DGETRI); used for diagnostics.
	
	lda = n; Ainv = A	
	call dgetrf ( n, n, Ainv, n, ipiv, info )
	if (INFO .ne. 0) write(6,*) "Error in dgetrf: INFO =", INFO
	call dgetri ( n, Ainv, n, ipiv, work, n, info )
	if (INFO .ne. 0) write(6,*) "Error in dgetrf: INFO =", INFO	
end subroutine inverse

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
