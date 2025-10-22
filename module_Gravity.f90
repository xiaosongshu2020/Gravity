module module_Gravity
    
    implicit none
    
    ! 引力场对象类型
    type :: GravityField
        character(len=50) :: name
        integer :: ID
        character(len=50) :: centralBody
        real*8 :: GM ! m^3/s^2
        real*8 :: referenceRadius ! m
        integer :: maxDegree
        integer :: maxOrder
        integer :: normalizationState
        real*8 :: referenceLongitude
        real*8 :: referenceLatitude
        real*8, allocatable :: C(:,:)  ! 球谐系数 C_lm
        real*8, allocatable :: S(:,:)  ! 球谐系数 S_lm
        logical :: coefficientsLoaded
    contains
        procedure :: init => initialize_gravity_field
        procedure :: read_coeff => read_gravity_coefficients
        procedure :: acceleration => compute_gravity_acceleration
        procedure :: show => get_gravity_field_info
    end type GravityField
    
    ! 模块参数
    integer, parameter :: MAX_DEGREE = 80
    integer, parameter :: MAX_ORDER = 80
    
    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 初始化引力场对象
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine initialize_gravity_field(this, name, central_body, max_degree, max_order)
        implicit none
        class(GravityField), intent(inout) :: this
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: central_body
        integer, intent(in) :: max_degree, max_order
        
        this%name = name
        this%centralBody = central_body
        this%maxDegree = max_degree
        this%maxOrder = max_order
        this%coefficientsLoaded = .false.
        
        ! 分配球谐系数数组
        if (allocated(this%C)) deallocate(this%C)
        if (allocated(this%S)) deallocate(this%S)
        
        allocate(this%C(0:max_degree, 0:max_order))
        allocate(this%S(0:max_degree, 0:max_order))
        
        this%C = 0.0d0
        this%S = 0.0d0
        
    end subroutine initialize_gravity_field
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 读取引力场球谐系数
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_gravity_coefficients(this, filename)
        implicit none
        class(GravityField), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer :: file_unit, iostat_val
        integer :: l, m
        real*8 :: C_val, S_val, C_uncertainty, S_uncertainty
        character(len=200) :: header_line
        
        ! 打开文件
        open(newunit=file_unit, file=trim(filename), status='old', action='read', iostat=iostat_val)
        if (iostat_val /= 0) then
            write(*,*) "Error opening file: ", trim(filename)
            return
        endif
        
        ! 读取第一行头部信息
        read(file_unit, '(A)', iostat=iostat_val) header_line
        if (iostat_val == 0) then
            ! 解析头部信息
            read(header_line, *) this%referenceRadius, &
                                this%GM, &
                                this%referenceLongitude, &
                                this%referenceLatitude
        endif
        
        ! 读取球谐系数
        do
            read(file_unit, *, iostat=iostat_val) l, m, C_val, S_val, C_uncertainty, S_uncertainty
            if (iostat_val /= 0) exit
            
            if (l <= this%maxDegree .and. m <= this%maxOrder) then
                this%C(l, m) = C_val
                this%S(l, m) = S_val
            endif
        end do
        
        close(file_unit)
        
        this%coefficientsLoaded = .true.
        write(*,*) "Gravity coefficients loaded successfully for ", trim(this%name)
        write(*,*) "Reference Radius: ", this%referenceRadius
        write(*,*) "GM: ", this%GM
        
    end subroutine read_gravity_coefficients
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 计算引力加速度（基于参考模块算法）
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine compute_gravity_acceleration(this, position, degree, order, acceleration)
    implicit none
    class(GravityField), intent(in) :: this
    real*8, intent(in) :: position(3)
    integer, intent(in) :: degree, order
    real*8, intent(out) :: acceleration(3)
        
    real*8 :: r, phi, lambda, sin_phi, cos_phi
    real*8 :: P(0:degree+1), Pd(0:degree+1, 0:order+1), Pp(0:degree+1), Pdp(0:degree+1, 0:order+1)
    real*8 :: k(3), G(3)
    real*8 :: F_zonal(3), F_tesseral(3), F(3)
    real*8 :: Rp
    integer :: l, m, n
    real*8 :: l_dble, m_dble
    real*8 :: cos_mlambda, sin_mlambda
    real*8 :: r_over_R
    real*8 :: rs_f(3), rs_f_norm, i_dble, j_dble
    integer :: i, j
        
    ! 检查系数是否已加载
    if (.not. this%coefficientsLoaded) then
        write(*,*) "Error: Gravity coefficients not loaded for ", trim(this%name)
        acceleration = 0.0d0
        return
    endif
        
    !--------------------------------------------------

    P = 0d0
    Pd = 0d0
    Pp = 0d0
    Pdp = 0d0

    !---------------------------------------------------

    ! m -> length_unit
    rs_f = position/this%referenceRadius

    !get longitude(lambda) and latitude(phi) in fixed coordinate
    rs_f_norm = norm2(rs_f)
    phi = dasin(rs_f(3)/rs_f_norm)
    call get_angle1(rs_f(2)/norm2(rs_f(1:2)), rs_f(1)/norm2(rs_f(1:2)), lambda)

    !get legendre polynomial(P) and its derivative(Pp)
    sin_phi = dsin(phi)

    P(0) = 1d0
    P(1) = dsqrt(3d0)*sin_phi
    n = degree
    do i = 2, n
        i_dble = dble(i)
        P(i) = dsqrt((2d0*i_dble + 1d0)/(2d0*i_dble - 1d0))*((2d0 - 1d0/i_dble)*sin_phi*P(i - 1) - &
        & dsqrt((2d0*i_dble - 1d0)/(2d0*i_dble - 3d0))*(1d0 - 1d0/i_dble)*P(i - 2))

         Pp(i) = i_dble/(1d0 - sin_phi**2)*(dsqrt((2d0*i_dble + 1d0)/(2d0*i_dble - 1d0))*P(i - 1) - sin_phi*P(i))

    end do

    !get zonal F
    k = (/0d0, 0d0, 1d0/)
    F_zonal = 0d0
    do i = 2, n
        i_dble = dble(i)
        F_zonal = F_zonal - this%C(i, 0)/rs_f_norm**(i + 3)*(((i_dble + 1d0)*P(i) + &
        & sin_phi*Pp(i))*rs_f - rs_f_norm*Pp(i)*k)

    end do

    !get the associative legendre polynomial(Pd) and its derivative(Pdp)
    Pd(1, 1) = dsqrt(3d0*(1d0 - sin_phi**2))
    do i = 2, n
        i_dble = dble(i)
        do j = 1, i - 1
            j_dble = dble(j)
            Pd(i, j) = dsqrt((2d0*i_dble + 1d0)*(2d0*i_dble - 1d0)/(i_dble + j_dble)/(i_dble - j_dble))*sin_phi*Pd(i - 1, j) - &
            & dsqrt((2d0*i_dble + 1d0)*(i_dble - 1d0 + j_dble)*(i_dble - 1d0 - j_dble)/(2d0*i_dble - 3d0)/(i_dble + j_dble)/(i_dble - j_dble))*Pd(i - 2, j)

        end do
        Pd(i, i) = dsqrt((2d0*i_dble + 1d0)/2d0/i_dble)*dsqrt(1d0 - sin_phi**2)*Pd(i - 1, i- 1)
    end do

    do i = 2, n
        i_dble = dble(i)
        do j = 1, i
            j_dble = dble(j)
            Pdp(i, j) = 1d0/dsqrt(1d0 - sin_phi**2)*(dsqrt((i_dble + j_dble + 1d0)*(i_dble - j_dble))*Pd(i, j + 1) - &
            & j_dble*sin_phi/dsqrt(1d0 - sin_phi**2)*Pd(i, j))

        end do
    end do

    !get tesseral F
    F_tesseral = 0d0
    G = (/-dsin(lambda), dcos(lambda), 0d0/)
    Rp = norm2(rs_f(1:2))
    do i = 2, n
        i_dble = dble(i)
        do j = 1, i
            j_dble = dble(j)
            !get sin_mlambda and cos_mlambda
            sin_mlambda = dsin(j_dble*lambda)
            cos_mlambda = dcos(j_dble*lambda)

            F_tesseral = F_tesseral - (1d0/rs_f_norm**(i_dble + 3d0)*(((i_dble + 1d0)*Pd(i, j) + sin_phi*Pdp(i, j))*rs_f - &
            & rs_f_norm*Pdp(i, j)*k)*(this%C(i, j)*cos_mlambda + this%S(i, j)*sin_mlambda) + &
            & j_dble/Rp/rs_f_norm**(i_dble + 1d0)*Pd(i, j)*(this%C(i, j)*sin_mlambda - this%S(i, j)*cos_mlambda)*G)

        end do
    end do

    F = F_zonal + F_tesseral
    
    ! -> unit: m/s^2
    acceleration = this%GM/(this%referenceRadius)**2 * F
        
    end subroutine compute_gravity_acceleration
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 获取引力场信息
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_gravity_field_info(this)
        implicit none
        class(GravityField), intent(in) :: this
        
        write(*,*) "=== Gravity Field Information ==="
        write(*,*) "Name: ", trim(this%name)
        write(*,*) "Central Body: ", trim(this%centralBody)
        write(*,*) "GM: ", this%GM
        write(*,*) "Reference Radius: ", this%referenceRadius
        write(*,*) "Max Degree: ", this%maxDegree
        write(*,*) "Max Order: ", this%maxOrder
        write(*,*) "Coefficients Loaded: ", this%coefficientsLoaded
        write(*,*) "=================================="
        
    end subroutine get_gravity_field_info
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This subroutine is written to get the value of the angle while 
    ! we know its sine and cosine.
    !The returned value of function 'dacos' is between 0 and pi, 
    ! however the value of the angle
    ! is between 0 and 2*pi
    subroutine get_angle1(sin_phi, cos_phi, phi)
    implicit none
    real*8 sin_phi, cos_phi, phi
    real*8 sin_phi2, cos_phi2
    !----------------------------------------------
    ! Normalize sin_phi and cos_phi
    real*8 norm
    norm = sqrt(sin_phi**2 + cos_phi**2)
    if (norm > 0d0) then
        sin_phi2 = sin_phi / norm
        cos_phi2 = cos_phi / norm
    end if
    phi = dacos(cos_phi2)
    if (sin_phi2 < 0) then
        phi = 2d0*3.141592653589793d0 - phi
    end if
    end subroutine

end module module_Gravity
