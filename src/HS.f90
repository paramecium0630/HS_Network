program HS_Network
    implicit none
    integer :: i, j, k, ii, jj, iread, progress_stride
    integer :: N, nstep, npast, NSTEP_PROTOCOL, nedge, nskip_edge
    integer :: LWORK, LWORK2, LDA, LDVL, LDVR
    integer :: nreal, iseed, idir, ntau, na, nb, NN
    integer, parameter :: small_blas_threshold = 16
    real(8) :: xKx1, xKx2, hs_functional, hs_functional_avg, exp_minus_hs_avg
    real(8) :: entropy_initial, entropy_final, entropy_change, entropy_change_avg
    real(8) :: logdet_K_start, logdet_K_end, entropy_change_logdet, entropy_change_quadratic_avg
    real(8) :: sum_heat, sum_heat_avg, sum_entropy, sum_entropy_avg, sum_work, sum_work_avg
    real(8) :: tt, tpast, protocol, protocol_next    
    real(8) :: dum, time1, time2, time1_HS, time2_HS
    real(8) :: DD, dt, prob, W_mean, W_sig, r0, w0, period, tau, bb
    real(8), parameter :: pi = 4d0*atan(1d0)
    logical :: use_blas_level2

    ! network variables
    integer, allocatable :: G0(:,:), indeg(:), outdeg(:)
    real(8), allocatable :: a(:), r(:), ar(:), noise(:,:), G(:,:), Q(:,:), Q_inv(:,:), noise_scale(:)
    real(8), allocatable :: inweightdeg(:), outweightdeg(:)

    ! dynamic variables
    real(8), allocatable :: f1(:), x1(:), noise_increment(:), fluct(:), fluct_prev(:), x_prev(:), f2(:)
    real(8), allocatable :: node_heat_step(:)!, xt(:,:), ft(:,:)
    real(8), allocatable :: f_avg(:,:), x_avg(:,:), K_avg(:,:,:), fixed_point(:,:)
    real(8), allocatable :: node_heat(:), node_entropy(:), node_work(:)
    ! real(8), allocatable :: protocol(:)

    ! schur variables
    real(8), allocatable :: QschurT(:,:), QschurVS(:,:), QschurWR(:), QschurWI(:), QschurWORK(:)
    real(8), allocatable :: QschurTMP(:,:), QschurC(:,:), QschurSol(:,:)
    real(8), allocatable :: Q_protocol_start(:,:), Q_protocol_end(:,:)
    real(8), allocatable :: lyapunov_residual_start(:,:), lyapunov_residual_end(:,:)
    logical, allocatable :: QschurBWORK(:)

    ! energetics
    real(8), allocatable :: node_heat_avg(:), node_entropy_avg(:), node_work_avg(:)
    real(8), allocatable :: Ko_protocol_start(:,:), Ko_protocol_end(:,:), Ko_inv_protocol(:,:,:), logdet_K(:)

    call read_variables_from_file()
    if (iread .eq. 1) then
        call infer_network_size_from_adjacency("input/adjacency_list.txt", N)
    elseif (iread .eq. 0) then
        N = NN
    else
        stop "iread must be 0 (random graph) or 1 (read adjacency_list.txt)."
    endif
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
    if (dt <= 0d0) stop "dt must be > 0."
    if (nreal <= 0) stop "nreal must be > 0."
    if (na < 1 .or. na > N) stop "na must be within 1..N."
    use_blas_level2 = (N > small_blas_threshold)
    allocate(G0(N,N), indeg(N), outdeg(N))
    allocate(a(N), r(N), ar(N), noise(N,N), G(N,N), Q(N,N), Q_inv(N,N), noise_scale(N))
    allocate(inweightdeg(N), outweightdeg(N))

    allocate(f1(N), x1(N), noise_increment(N), fluct(N), fluct_prev(N), x_prev(N), f2(N), node_heat_step(N))
    allocate(node_heat(N), node_entropy(N))
    allocate(f_avg(N, -NPAST:NSTEP), x_avg(N, -NPAST:NSTEP), K_avg(N, N, -NPAST:NSTEP))

    allocate(QschurT(N,N), QschurVS(N,N), QschurWR(N), QschurWI(N), QschurWORK(LWORK))
    allocate(QschurBWORK(LWORK2))
    allocate(QschurTMP(N,N), QschurC(N,N), QschurSol(N,N))
    allocate(Q_protocol_start(N,N), Q_protocol_end(N,N))
    allocate(lyapunov_residual_start(N,N), lyapunov_residual_end(N,N))
    allocate(Ko_protocol_start(N,N), Ko_protocol_end(N,N))
    allocate(Ko_inv_protocol(N,N,0:NSTEP_PROTOCOL), logdet_K(0:NSTEP_PROTOCOL), fixed_point(N,0:NSTEP_PROTOCOL))

    allocate(node_heat_avg(N), node_entropy_avg(N), node_work_avg(N))

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
    OPEN(UNIT=21, FILE='output/node_heat.dat', STATUS='REPLACE')
    OPEN(UNIT=22, FILE='output/node_entropy.dat', STATUS='REPLACE')
    OPEN(UNIT=24, FILE='output/fixed_point.dat', STATUS='REPLACE')

    G = 0d0; G0 = 0
    if (iread .eq. 1) then
        call read_adjacency_list("input/adjacency_list.txt", N, w0, G0, G, nskip_edge)
        if (nskip_edge > 0) then
            write(6,*) 'skipped adjacency edges outside 1..N = ', nskip_edge
        endif
    elseif (iread .eq. 0) then
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
    else
        stop "iread must be 0 (random graph) or 1 (read adjacency_list.txt)."
    endif

    do i = 1, n; do j = 1, n
        if (G0(i,j) == 0) cycle
        write(11,*) i, j, G0(i,j)
        write(12,*) i, j, G(i,j)
    enddo; enddo

    nedge = 0
	do i = 1, N; do j = 1, N
	   nedge = nedge + G0(i,j)
	enddo; enddo
    if (N > 1) then
	    print*, 'connection probability: ', 1.d0*nedge/(N*(N-1))
    else
        print*, 'single-node network: connection probability is undefined.'
    endif
        
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
        ! r(i) = r0*(ran3(iseed)+0.5d0)
        r(i) = r0
        noise(i,i) = dd*i
        ! noise(i,i) = dd*(ran3(iseed)+0.5d0)
    enddo

    ar = r*a

    do i = 1, N
        if (noise(i,i) < 0d0) stop "noise diagonal must be >= 0."
        if (noise(i,i) > 0d0) then
            noise_scale(i) = sqrt(noise(i,i)*dt)
        else
            noise_scale(i) = 0d0
        endif
    enddo

    do i = 1, N
        write(14,*) i, r(i), a(i)
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
    x_avg = 0; f_avg = 0; K_avg = 0
    node_heat = 0d0; node_entropy = 0d0; node_work = 0d0
    node_heat_avg = 0d0; node_entropy_avg = 0d0; node_work_avg = 0d0
    progress_stride = max(1, nreal/10)
    do k = 1, nreal
        if (mod(k-1, progress_stride) == 0) call progress(k, nreal)
        protocol = 0; protocol_next = 0; x1 = 0; f1 = 0; noise_increment = 0

        do ii = -npast, NSTEP
            tt = ii*dt
            do i = 1, N
                if (noise_scale(i) > 0d0) then
                    noise_increment(i) = noise_scale(i)*GAUDEV(0d0, 1.d0, iseed)
                else
                    noise_increment(i) = 0d0
                endif
            enddo
            call protocol_at_time(tt, tau, bb, protocol)
            call protocol_at_time(tt + dt, tau, bb, protocol_next)

            if (k == 1) then
                if (mod(abs(ii), 100) .eq. 1) then                    
                    write(17,*) tt, protocol
                endif
            endif

            x_avg(:, ii) = x_avg(:, ii) + x1
            if (use_blas_level2) then
                call dger(N, N, 1d0, x1, 1, x1, 1, K_avg(:,:,ii), N)
            else
                do j = 1, N
                    do i = 1, N
                        K_avg(i,j,ii) = K_avg(i,j,ii) + x1(i)*x1(j)
                    enddo
                enddo
            endif
            call advance_hs_state_heun(N, r, a, G, protocol, protocol_next, noise_increment, f1, x1)
            f_avg(:, ii) = f_avg(:, ii) + f1
        enddo
    enddo
    call cpu_time(time2)

    write(*,*) ""
    write(*,*) 'CPU time for simulation: ', time2 - time1, ' seconds.'

    x_avg = x_avg/dble(nreal)
    f_avg = f_avg/dble(nreal)
    K_avg = K_avg/dble(nreal)

    if (use_blas_level2) then
        do jj = -npast, nstep
            call dger(N, N, -1d0, x_avg(:,jj), 1, x_avg(:,jj), 1, K_avg(:,:,jj), N)
        enddo
    else
        do jj = -npast, nstep
            do j = 1, N
                do i = 1, N
                    K_avg(i,j,jj) = K_avg(i,j,jj) - x_avg(i,jj)*x_avg(j,jj)
                enddo
            enddo
        enddo
    endif

    do ii = 0, NSTEP_PROTOCOL
        tt = ii*dt
        call protocol_at_time(tt, tau, bb, protocol)
        call solve_fixed_point_from_protocol(N, r, a, ar, G, protocol, fixed_point(:,ii), Q, Q_inv)
        call schur_factorize_real_matrix(N, Q, QschurT, QschurVS, QschurWR, QschurWI, QschurWORK, QschurBWORK)
        call solve_qx_xqt_from_schur_real_ws(N, QschurT, QschurVS, -noise, QschurTMP, QschurC, QschurSol)
        call matrix_logdet_spd(QschurSol, N, logdet_K(ii)) ! logdet of Ko at protocol step ii
        call inverse(N, QschurSol, Ko_inv_protocol(:,:,ii))
        if (ii == 0) then
            Q_protocol_start = Q
            Ko_protocol_start = QschurSol
        endif
        if (ii == NSTEP_PROTOCOL) then
            Q_protocol_end = Q
            Ko_protocol_end = QschurSol
        endif
    enddo

    call cpu_time(time1_HS)
    ! HATANO-SASA
    xKx1 = 0d0; xKx2 = 0d0
    hs_functional_avg = 0d0; exp_minus_hs_avg = 0d0; entropy_change_avg = 0d0
    sum_heat_avg = 0d0; sum_entropy_avg = 0d0
    do k = 1, nreal
        if (mod(k-1, progress_stride) == 0) call progress(k, nreal)
        hs_functional = 0d0; entropy_initial = 0d0; entropy_final = 0d0; entropy_change = 0d0
        sum_heat = 0d0; sum_entropy = 0d0
        protocol = 0; protocol_next = 0
        x1 = 0; f1 = 0; noise_increment = 0
        node_heat = 0d0; node_entropy = 0d0

        do ii = -npast, NSTEP_PROTOCOL
            if (ii == 0) then
                fluct = x1 - x_avg(:, 0) ! fluctuation at initial time
                xKx1 = DOT_PRODUCT(fluct, MATMUL(Ko_inv_protocol(:,:,0), fluct)) ! delta x^T * Ko^{-1} * delta x
                entropy_initial = 0.5d0*logdet_K(0) + 0.5d0*xKx1
            endif

            if (ii >= 1) then
                fluct = x1 - x_avg(:, ii) ! fluctuation at time step ii
                fluct_prev = x1 - x_avg(:, ii-1)
                xkx1 = DOT_PRODUCT(fluct, MATMUL(Ko_inv_protocol(:,:,ii), fluct))
                xkx2 = DOT_PRODUCT(fluct_prev, MATMUL(Ko_inv_protocol(:,:,ii-1), fluct_prev))
                hs_functional = hs_functional + 0.5d0*(logdet_K(ii) - logdet_K(ii-1)) + 0.5d0*(xkx1 - xkx2)
            endif

            if (ii == NSTEP_PROTOCOL) then
                fluct = x1 - x_avg(:, NSTEP_PROTOCOL) ! fluctuation at final time
                xKx1 = DOT_PRODUCT(fluct, MATMUL(Ko_inv_protocol(:,:,NSTEP_PROTOCOL), fluct))
                entropy_final = 0.5d0*logdet_K(NSTEP_PROTOCOL) + 0.5d0*xKx1
                entropy_change = entropy_final - entropy_initial ! Delta \phi
                exit ! no need to advance state at last protocol step
            endif

            tt = ii*dt
            do i = 1, N
                if (noise_scale(i) > 0d0) then
                    noise_increment(i) = noise_scale(i)*GAUDEV(0d0, 1.d0, iseed)
                else
                    noise_increment(i) = 0d0
                endif
            enddo
            call protocol_at_time(tt, tau, bb, protocol)
            call protocol_at_time(tt + dt, tau, bb, protocol_next)
            x_prev = x1
            call advance_hs_state_heun(N, r, a, G, protocol, protocol_next, noise_increment, f1, x1)
            if (ii >= 0 .and. ii < NSTEP_PROTOCOL) then
                call calculate_hs_force(N, r, a, G, protocol_next, x1, f2)
                node_heat_step = -0.5d0*(f1 + f2)*(x1 - x_prev) ! Stratonovich heat increment at each node
                node_heat = node_heat + node_heat_step
                do i = 1, N
                    if (noise(i,i) > 0d0) then
                        node_entropy(i) = node_entropy(i) - 2d0/noise(i,i)*node_heat_step(i)
                    endif
                enddo
            endif
        enddo
        hs_functional_avg = hs_functional_avg + hs_functional
        exp_minus_hs_avg = exp_minus_hs_avg + exp(-hs_functional)
        entropy_change_avg = entropy_change_avg + entropy_change
        node_heat_avg = node_heat_avg + node_heat
        node_entropy_avg = node_entropy_avg + node_entropy
        sum_heat = sum(node_heat)
        sum_entropy = sum(node_entropy)
        sum_heat_avg = sum_heat_avg + sum_heat
        sum_entropy_avg = sum_entropy_avg + sum_entropy
    enddo
    call cpu_time(time2_HS)
    hs_functional_avg = hs_functional_avg/nreal ! Y-value average
    exp_minus_hs_avg = exp_minus_hs_avg/nreal ! exp(-Y) average, should be close to 1 by Hatano-Sasa equality
    entropy_change_avg = entropy_change_avg/nreal ! Delta \phi average
    node_heat_avg = node_heat_avg/nreal ! average heat at each node
    sum_heat_avg = sum_heat_avg/nreal ! average sum over node_heat
    node_entropy_avg = node_entropy_avg/nreal ! average entropy at each node
    sum_entropy_avg = sum_entropy_avg/nreal ! average sum over node_entropy

    call matrix_logdet_spd(Ko_protocol_start, N, logdet_K_start)
    call matrix_logdet_spd(Ko_protocol_end, N, logdet_K_end)
    entropy_change_logdet = 0.5d0*(logdet_K_end - logdet_K_start)
    entropy_change_quadratic_avg = entropy_change_avg - entropy_change_logdet
    write(*,*) 'CPU time for Hatano-Sasa: ', time2_HS - time1_HS, ' seconds.'
    write(6,*) 'Hatano-Sasa Y average = ', hs_functional_avg
    write(6,*) 'Hatano-Sasa exp(-Y) average = ', exp_minus_hs_avg
    write(6,*) 'Entropy change average = ', entropy_change_avg
    write(6,*) 'Logdet Ko start = ', logdet_K_start
    write(6,*) 'Logdet Ko end = ', logdet_K_end
    ! write(6,*) 'Entropy change logdet part = ', entropy_change_logdet
    ! write(6,*) 'Entropy change quadratic average = ', entropy_change_quadratic_avg
    write(6,*) 'Excess entropy average = ', hs_functional_avg - entropy_change_avg
    write(6,*) 'Sum heat average = ', sum_heat_avg
    write(6,*) 'Sum entropy average = ', sum_entropy_avg

    do i = 1, N
        write(21,*) i, node_heat_avg(i)
        write(22,*) i, node_entropy_avg(i)
    enddo

    lyapunov_residual_start = MATMUL(Q_protocol_start, Ko_protocol_start) &
        + MATMUL(Ko_protocol_start, TRANSPOSE(Q_protocol_start)) + noise
    lyapunov_residual_end = MATMUL(Q_protocol_end, Ko_protocol_end) &
        + MATMUL(Ko_protocol_end, TRANSPOSE(Q_protocol_end)) + noise

    do i = -npast/5, nstep
        if (mod(i+npast, 10) .eq. 1) then  
            write(15,*) i*dt, x_avg(1,i), x_avg(2,i)            
            write(16,*) i*dt, f_avg(1,i), f_avg(2,i)
            write(18,*) i*dt, K_avg(1,1,i), K_avg(1,2,i)
            if (i >= 0 .and. i < NSTEP_PROTOCOL) then
                write(24,*) i*dt, fixed_point(1,i), fixed_point(2,i)
            endif
        endif
    enddo

    do i = 1, N
        do j = i, N
            write(19,'(*(G0,1X))') K_avg(i,j,0), real(Ko_protocol_start(i,j)), K_avg(i,j,NSTEP_PROTOCOL), &
                real(Ko_protocol_end(i,j))
            write(20,'(*(G0,1X))') lyapunov_residual_start(i,j), lyapunov_residual_end(i,j)
        enddo
    enddo

contains

subroutine progress(i, nprog)
	implicit none
	integer, intent(in) :: i, nprog
    real(8) :: val
		
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
    read(10, *) NN, iread ! NN is used for random graph; iread=1 infers N from adjacency_list.txt
    close(10)
end subroutine read_variables_from_file

subroutine infer_network_size_from_adjacency(filename, n_inferred)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(out) :: n_inferred
    integer :: row_node, col_node, edge_weight

    n_inferred = 0
    open(unit=998, file=filename, status="old", action="read")
    do
        read(998, *, end=110) row_node, col_node, edge_weight
        if (row_node < 1 .or. col_node < 1) cycle
        n_inferred = max(n_inferred, row_node, col_node)
    enddo
110 continue
    close(998)

    if (n_inferred <= 0) stop "Could not infer N from adjacency_list.txt."
end subroutine infer_network_size_from_adjacency

subroutine read_adjacency_list(filename, n, w0, G0, G, nskip_edge)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: n
    real(8), intent(in) :: w0
    integer, intent(out) :: G0(n,n), nskip_edge
    real(8), intent(out) :: G(n,n)
    integer :: row_node, col_node, edge_weight

    G0 = 0
    G = 0d0
    nskip_edge = 0
    open(unit=999, file=filename, status="old", action="read")
    do
        read(999, *, end=120) row_node, col_node, edge_weight
        if (row_node < 1 .or. row_node > n .or. col_node < 1 .or. col_node > n) then
            nskip_edge = nskip_edge + 1
            cycle
        endif
        G0(row_node,col_node) = edge_weight
        G(row_node,col_node) = edge_weight*w0
    enddo
120 continue
    close(999)
end subroutine read_adjacency_list

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
    real(8) :: coupling_sum, local_drift, xjn

    if (use_blas_level2) then
        call dgemv('N', n, n, 1d0, G, n, x, 1, 0d0, f, 1)
        do jn = 1, n
            local_drift = -r(jn)*(x(jn)-a(jn))
            f(jn) = f(jn) - inweightdeg(jn)*x(jn) + local_drift
        enddo
    else
        do jn = 1, n
            xjn = x(jn)
            coupling_sum = 0d0
            do j = 1, n
                coupling_sum = coupling_sum + G(jn,j)*(x(j)-xjn)
            enddo
            f(jn) = -r(jn)*(xjn-a(jn)) + coupling_sum
        enddo
    endif
    f(na) = f(na) - protocol*(x(na)-a(na))

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
		real(8) :: inweightdegree_i

		! HS local dynamics:
		!   f_i(x_i) = -r_i*(x_i-a_i)
		! gives f_i'(x_i) = -r_i.
        ! The protocol shifts the relaxation rate of node na.
		! Diffusive coupling h(x_i, x_j) = x_j - x_i contributes:
		!   off-diagonal:  +G(i,j)
		!   diagonal:      -sum_j G(i,j)
		Q = G
		do i = 1, N
			inweightdegree_i = inweightdeg(i)
			Q(i,i) = Q(i,i) - r(i) - inweightdegree_i
		enddo
        Q(na,na) = Q(na,na) - protocol
	end subroutine build_jacobian_from_protocol

subroutine solve_fixed_point_from_protocol(N, r, a, ar, G, protocol, x_fixed, Q, Q_inv)
        implicit none
        integer, intent(in) :: N
        real(8), intent(in) :: r(N), a(N), ar(N), G(N,N), protocol
        real(8), intent(out) :: x_fixed(N), Q(N,N), Q_inv(N,N)
        real(8) :: rhs_fp(N)

        call build_jacobian_from_protocol(N, r, G, protocol, Q)
        rhs_fp = ar
        rhs_fp(na) = rhs_fp(na) + protocol*a(na)
        call inverse(N, Q, Q_inv)
        x_fixed = -MATMUL(Q_inv, rhs_fp)
end subroutine solve_fixed_point_from_protocol

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
