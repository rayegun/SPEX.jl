const mp_limb_t = Culong

mutable struct __mpz_struct
    _mp_alloc::Cint
    _mp_size::Cint
    _mp_d::Ptr{mp_limb_t}
    __mpz_struct() = new()
end

const mpz_t = NTuple{1, __mpz_struct}

mutable struct __mpq_struct
    _mp_num::__mpz_struct
    _mp_den::__mpz_struct
    __mpq_struct() = new()
end

const mpq_t = NTuple{1, __mpq_struct}

@enum mpfr_rnd_t::Int32 begin
    MPFR_RNDN = 0
    MPFR_RNDZ = 1
    MPFR_RNDU = 2
    MPFR_RNDD = 3
    MPFR_RNDA = 4
    MPFR_RNDF = 5
    MPFR_RNDNA = -1
end

const mpfr_prec_t = Clong

const mpfr_sign_t = Cint

const mpfr_exp_t = Clong

mutable struct __mpfr_struct
    _mpfr_prec::mpfr_prec_t
    _mpfr_sign::mpfr_sign_t
    _mpfr_exp::mpfr_exp_t
    _mpfr_d::Ptr{mp_limb_t}
    __mpfr_struct() = new()
end

const mpfr_t = NTuple{1, __mpfr_struct}

function SPEX_free(p)
    @ccall libspexutil.SPEX_free(p::Ptr{Cvoid})::Cvoid
end

mutable struct SuiteSparse_config_struct
    malloc_func::Ptr{Cvoid}
    calloc_func::Ptr{Cvoid}
    realloc_func::Ptr{Cvoid}
    free_func::Ptr{Cvoid}
    printf_func::Ptr{Cvoid}
    hypot_func::Ptr{Cvoid}
    divcomplex_func::Ptr{Cvoid}
    SuiteSparse_config_struct() = new()
end

function SuiteSparse_start()
    @ccall libxxx.SuiteSparse_start()::Cvoid
end

function SuiteSparse_finish()
    @ccall libxxx.SuiteSparse_finish()::Cvoid
end

function SuiteSparse_malloc(nitems, size_of_item)
    @ccall libxxx.SuiteSparse_malloc(nitems::Csize_t, size_of_item::Csize_t)::Ptr{Cvoid}
end

function SuiteSparse_calloc(nitems, size_of_item)
    @ccall libxxx.SuiteSparse_calloc(nitems::Csize_t, size_of_item::Csize_t)::Ptr{Cvoid}
end

function SuiteSparse_realloc(nitems_new, nitems_old, size_of_item, p, ok)
    @ccall libxxx.SuiteSparse_realloc(nitems_new::Csize_t, nitems_old::Csize_t, size_of_item::Csize_t, p::Ptr{Cvoid}, ok::Ptr{Cint})::Ptr{Cvoid}
end

function SuiteSparse_free(p)
    @ccall libxxx.SuiteSparse_free(p::Ptr{Cvoid})::Ptr{Cvoid}
end

function SuiteSparse_tic(tic)
    @ccall libxxx.SuiteSparse_tic(tic::Ptr{Cdouble})::Cvoid
end

function SuiteSparse_toc(tic)
    @ccall libxxx.SuiteSparse_toc(tic::Ptr{Cdouble})::Cdouble
end

function SuiteSparse_time()
    @ccall libxxx.SuiteSparse_time()::Cdouble
end

function SuiteSparse_hypot(x, y)
    @ccall libxxx.SuiteSparse_hypot(x::Cdouble, y::Cdouble)::Cdouble
end

function SuiteSparse_divcomplex(ar, ai, br, bi, cr, ci)
    @ccall libxxx.SuiteSparse_divcomplex(ar::Cdouble, ai::Cdouble, br::Cdouble, bi::Cdouble, cr::Ptr{Cdouble}, ci::Ptr{Cdouble})::Cint
end

function SuiteSparse_version(version)
    @ccall libxxx.SuiteSparse_version(version::Ptr{Cint})::Cint
end

@enum SPEX_info::Int32 begin
    SPEX_OK = 0
    SPEX_OUT_OF_MEMORY = -1
    SPEX_SINGULAR = -2
    SPEX_INCORRECT_INPUT = -3
    SPEX_INCORRECT = -4
    SPEX_UNSYMMETRIC = -5
    SPEX_PANIC = -6
end

@enum SPEX_pivot::UInt32 begin
    SPEX_SMALLEST = 0
    SPEX_DIAGONAL = 1
    SPEX_FIRST_NONZERO = 2
    SPEX_TOL_SMALLEST = 3
    SPEX_TOL_LARGEST = 4
    SPEX_LARGEST = 5
end

@enum SPEX_col_order::UInt32 begin
    SPEX_NO_ORDERING = 0
    SPEX_COLAMD = 1
    SPEX_AMD = 2
end

struct SPEX_options
    pivot::SPEX_pivot
    order::SPEX_col_order
    tol::Cdouble
    print_level::Cint
    prec::Int32
    round::mpfr_rnd_t
    check::Bool
end

function SPEX_create_default_options(option)
    @ccall libspexutil.SPEX_create_default_options(option::Ptr{Ptr{SPEX_options}})::SPEX_info
end

@enum SPEX_kind::UInt32 begin
    SPEX_CSC = 0
    SPEX_TRIPLET = 1
    SPEX_DENSE = 2
end

@enum SPEX_type::UInt32 begin
    SPEX_MPZ = 0
    SPEX_MPQ = 1
    SPEX_MPFR = 2
    SPEX_INT64 = 3
    SPEX_FP64 = 4
end

struct __JL_Ctag_44
    data::NTuple{8, UInt8}
end

function Base.getproperty(x::Ptr{__JL_Ctag_44}, f::Symbol)
    f === :mpz && return Ptr{Ptr{mpz_t}}(x + 0)
    f === :mpq && return Ptr{Ptr{mpq_t}}(x + 0)
    f === :mpfr && return Ptr{Ptr{mpfr_t}}(x + 0)
    f === :int64 && return Ptr{Ptr{Int64}}(x + 0)
    f === :fp64 && return Ptr{Ptr{Cdouble}}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::__JL_Ctag_44, f::Symbol)
    r = Ref{__JL_Ctag_44}(x)
    ptr = Base.unsafe_convert(Ptr{__JL_Ctag_44}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{__JL_Ctag_44}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

mutable struct SPEX_matrix
    m::Int64
    n::Int64
    nzmax::Int64
    nz::Int64
    kind::SPEX_kind
    type::SPEX_type
    p::Ptr{Int64}
    p_shallow::Bool
    i::Ptr{Int64}
    i_shallow::Bool
    j::Ptr{Int64}
    j_shallow::Bool
    x::__JL_Ctag_44
    x_shallow::Bool
    scale::mpq_t
    SPEX_matrix() = new()
end

function SPEX_matrix_allocate(A_handle, kind, type, m, n, nzmax, shallow, init, option)
    @ccall libspexutil.SPEX_matrix_allocate(A_handle::Ptr{Ptr{SPEX_matrix}}, kind::SPEX_kind, type::SPEX_type, m::Int64, n::Int64, nzmax::Int64, shallow::Bool, init::Bool, option::Ptr{SPEX_options})::SPEX_info
end

function SPEX_matrix_free(A_handle, option)
    @ccall libspexutil.SPEX_matrix_free(A_handle::Ptr{Ptr{SPEX_matrix}}, option::Ptr{SPEX_options})::SPEX_info
end

function SPEX_matrix_nnz(nnz, A, option)
    @ccall libspexutil.SPEX_matrix_nnz(nnz::Ptr{Int64}, A::Ptr{SPEX_matrix}, option::Ptr{SPEX_options})::SPEX_info
end

function SPEX_matrix_copy(C_handle, C_kind, C_type, A, option)
    @ccall libspexutil.SPEX_matrix_copy(C_handle::Ptr{Ptr{SPEX_matrix}}, C_kind::SPEX_kind, C_type::SPEX_type, A::Ptr{SPEX_matrix}, option::Ptr{SPEX_options})::SPEX_info
end

mutable struct SPEX_LU_analysis
    q::Ptr{Int64}
    lnz::Int64
    unz::Int64
    SPEX_LU_analysis() = new()
end

function SPEX_LU_analysis_free(S, option)
    @ccall libspexutil.SPEX_LU_analysis_free(S::Ptr{Ptr{SPEX_LU_analysis}}, option::Ptr{SPEX_options})::SPEX_info
end

function SPEX_LU_analyze(S, A, option)
    @ccall libspexutil.SPEX_LU_analyze(S::Ptr{Ptr{SPEX_LU_analysis}}, A::Ptr{SPEX_matrix}, option::Ptr{SPEX_options})::SPEX_info
end

function SPEX_calloc(nitems, size)
    @ccall libspexutil.SPEX_calloc(nitems::Csize_t, size::Csize_t)::Ptr{Cvoid}
end

function SPEX_malloc(size)
    @ccall libspexutil.SPEX_malloc(size::Csize_t)::Ptr{Cvoid}
end

function SPEX_realloc(nitems_new, nitems_old, size_of_item, p, ok)
    @ccall libspexutil.SPEX_realloc(nitems_new::Int64, nitems_old::Int64, size_of_item::Csize_t, p::Ptr{Cvoid}, ok::Ptr{Bool})::Ptr{Cvoid}
end

function SPEX_initialize()
    @ccall libspexutil.SPEX_initialize()::SPEX_info
end

function SPEX_initialize_expert(MyMalloc, MyCalloc, MyRealloc, MyFree)
    @ccall libspexutil.SPEX_initialize_expert(MyMalloc::Ptr{Cvoid}, MyCalloc::Ptr{Cvoid}, MyRealloc::Ptr{Cvoid}, MyFree::Ptr{Cvoid})::SPEX_info
end

function SPEX_finalize()
    @ccall libspexutil.SPEX_finalize()::SPEX_info
end

function SPEX_matrix_check(A, option)
    @ccall libspexutil.SPEX_matrix_check(A::Ptr{SPEX_matrix}, option::Ptr{SPEX_options})::SPEX_info
end

function SPEX_matrix_div(x2_handle, x, scalar, option)
    @ccall libspexutil.SPEX_matrix_div(x2_handle::Ptr{Ptr{SPEX_matrix}}, x::Ptr{SPEX_matrix}, scalar::Ptr{Cvoid}, option::Ptr{SPEX_options})::SPEX_info
end

function SPEX_matrix_mul(x, scalar)
    @ccall libspexutil.SPEX_matrix_mul(x::Ptr{SPEX_matrix}, scalar::Ptr{Cvoid})::SPEX_info
end

function SPEX_check_solution(A, x, b, option)
    @ccall libspexutil.SPEX_check_solution(A::Ptr{SPEX_matrix}, x::Ptr{SPEX_matrix}, b::Ptr{SPEX_matrix}, option::Ptr{SPEX_options})::SPEX_info
end

function SPEX_cumsum(p, c, n)
    @ccall libspexutil.SPEX_cumsum(p::Ptr{Int64}, c::Ptr{Int64}, n::Int64)::SPEX_info
end

function SPEX_mpz_init(x)
    @ccall libspexutil.SPEX_mpz_init(x::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpz_init2(x, size)
    @ccall libspexutil.SPEX_mpz_init2(x::Ptr{Cvoid}, size::Csize_t)::SPEX_info
end

function SPEX_mpz_set(x, y)
    @ccall libspexutil.SPEX_mpz_set(x::Ptr{Cvoid}, y::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpz_set_ui(x, y)
    @ccall libspexutil.SPEX_mpz_set_ui(x::Ptr{Cvoid}, y::UInt64)::SPEX_info
end

function SPEX_mpz_set_si(x, y)
    @ccall libspexutil.SPEX_mpz_set_si(x::Ptr{Cvoid}, y::Int64)::SPEX_info
end

function SPEX_mpz_get_d(x, y)
    @ccall libspexutil.SPEX_mpz_get_d(x::Ptr{Cdouble}, y::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpz_get_si(x, y)
    @ccall libspexutil.SPEX_mpz_get_si(x::Ptr{Int64}, y::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpz_set_q(x, y)
    @ccall libspexutil.SPEX_mpz_set_q(x::Ptr{Cvoid}, y::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpz_mul(a, b, c)
    @ccall libspexutil.SPEX_mpz_mul(a::Ptr{Cvoid}, b::Ptr{Cvoid}, c::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpz_submul(x, y, z)
    @ccall libspexutil.SPEX_mpz_submul(x::Ptr{Cvoid}, y::Ptr{Cvoid}, z::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpz_divexact(x, y, z)
    @ccall libspexutil.SPEX_mpz_divexact(x::Ptr{Cvoid}, y::Ptr{Cvoid}, z::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpz_gcd(x, y, z)
    @ccall libspexutil.SPEX_mpz_gcd(x::Ptr{Cvoid}, y::Ptr{Cvoid}, z::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpz_lcm(lcm, x, y)
    @ccall libspexutil.SPEX_mpz_lcm(lcm::Ptr{Cvoid}, x::Ptr{Cvoid}, y::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpz_abs(x, y)
    @ccall libspexutil.SPEX_mpz_abs(x::Ptr{Cvoid}, y::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpz_cmp(r, x, y)
    @ccall libspexutil.SPEX_mpz_cmp(r::Ptr{Cint}, x::Ptr{Cvoid}, y::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpz_cmpabs(r, x, y)
    @ccall libspexutil.SPEX_mpz_cmpabs(r::Ptr{Cint}, x::Ptr{Cvoid}, y::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpz_cmp_ui(r, x, y)
    @ccall libspexutil.SPEX_mpz_cmp_ui(r::Ptr{Cint}, x::Ptr{Cvoid}, y::UInt64)::SPEX_info
end

function SPEX_mpz_sgn(sgn, x)
    @ccall libspexutil.SPEX_mpz_sgn(sgn::Ptr{Cint}, x::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpz_sizeinbase(size, x, base)
    @ccall libspexutil.SPEX_mpz_sizeinbase(size::Ptr{Csize_t}, x::Ptr{Cvoid}, base::Int64)::SPEX_info
end

function SPEX_mpq_init(x)
    @ccall libspexutil.SPEX_mpq_init(x::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpq_set(x, y)
    @ccall libspexutil.SPEX_mpq_set(x::Ptr{Cvoid}, y::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpq_set_z(x, y)
    @ccall libspexutil.SPEX_mpq_set_z(x::Ptr{Cvoid}, y::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpq_set_d(x, y)
    @ccall libspexutil.SPEX_mpq_set_d(x::Ptr{Cvoid}, y::Cdouble)::SPEX_info
end

function SPEX_mpq_set_ui(x, y, z)
    @ccall libspexutil.SPEX_mpq_set_ui(x::Ptr{Cvoid}, y::UInt64, z::UInt64)::SPEX_info
end

function SPEX_mpq_set_si(x, y, z)
    @ccall libspexutil.SPEX_mpq_set_si(x::Ptr{Cvoid}, y::Int64, z::UInt64)::SPEX_info
end

function SPEX_mpq_set_num(x, y)
    @ccall libspexutil.SPEX_mpq_set_num(x::Ptr{Cvoid}, y::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpq_set_den(x, y)
    @ccall libspexutil.SPEX_mpq_set_den(x::Ptr{Cvoid}, y::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpq_get_den(x, y)
    @ccall libspexutil.SPEX_mpq_get_den(x::Ptr{Cvoid}, y::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpq_get_d(x, y)
    @ccall libspexutil.SPEX_mpq_get_d(x::Ptr{Cdouble}, y::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpq_abs(x, y)
    @ccall libspexutil.SPEX_mpq_abs(x::Ptr{Cvoid}, y::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpq_add(x, y, z)
    @ccall libspexutil.SPEX_mpq_add(x::Ptr{Cvoid}, y::Ptr{Cvoid}, z::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpq_mul(x, y, z)
    @ccall libspexutil.SPEX_mpq_mul(x::Ptr{Cvoid}, y::Ptr{Cvoid}, z::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpq_div(x, y, z)
    @ccall libspexutil.SPEX_mpq_div(x::Ptr{Cvoid}, y::Ptr{Cvoid}, z::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpq_cmp(r, x, y)
    @ccall libspexutil.SPEX_mpq_cmp(r::Ptr{Cint}, x::Ptr{Cvoid}, y::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpq_cmp_ui(r, x, num, den)
    @ccall libspexutil.SPEX_mpq_cmp_ui(r::Ptr{Cint}, x::Ptr{Cvoid}, num::UInt64, den::UInt64)::SPEX_info
end

function SPEX_mpq_sgn(sgn, x)
    @ccall libspexutil.SPEX_mpq_sgn(sgn::Ptr{Cint}, x::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpq_equal(r, x, y)
    @ccall libspexutil.SPEX_mpq_equal(r::Ptr{Cint}, x::Ptr{Cvoid}, y::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpfr_init2(x, size)
    @ccall libspexutil.SPEX_mpfr_init2(x::Ptr{Cvoid}, size::UInt64)::SPEX_info
end

function SPEX_mpfr_set(x, y, rnd)
    @ccall libspexutil.SPEX_mpfr_set(x::Ptr{Cvoid}, y::Ptr{Cvoid}, rnd::mpfr_rnd_t)::SPEX_info
end

function SPEX_mpfr_set_d(x, y, rnd)
    @ccall libspexutil.SPEX_mpfr_set_d(x::Ptr{Cvoid}, y::Cdouble, rnd::mpfr_rnd_t)::SPEX_info
end

function SPEX_mpfr_set_si(x, y, rnd)
    @ccall libspexutil.SPEX_mpfr_set_si(x::Ptr{Cvoid}, y::Int64, rnd::mpfr_rnd_t)::SPEX_info
end

function SPEX_mpfr_set_q(x, y, rnd)
    @ccall libspexutil.SPEX_mpfr_set_q(x::Ptr{Cvoid}, y::Ptr{Cvoid}, rnd::mpfr_rnd_t)::SPEX_info
end

function SPEX_mpfr_set_z(x, y, rnd)
    @ccall libspexutil.SPEX_mpfr_set_z(x::Ptr{Cvoid}, y::Ptr{Cvoid}, rnd::mpfr_rnd_t)::SPEX_info
end

function SPEX_mpfr_get_z(x, y, rnd)
    @ccall libspexutil.SPEX_mpfr_get_z(x::Ptr{Cvoid}, y::Ptr{Cvoid}, rnd::mpfr_rnd_t)::SPEX_info
end

function SPEX_mpfr_get_q(x, y, rnd)
    @ccall libspexutil.SPEX_mpfr_get_q(x::Ptr{Cvoid}, y::Ptr{Cvoid}, rnd::mpfr_rnd_t)::SPEX_info
end

function SPEX_mpfr_get_d(x, y, rnd)
    @ccall libspexutil.SPEX_mpfr_get_d(x::Ptr{Cdouble}, y::Ptr{Cvoid}, rnd::mpfr_rnd_t)::SPEX_info
end

function SPEX_mpfr_get_si(x, y, rnd)
    @ccall libspexutil.SPEX_mpfr_get_si(x::Ptr{Int64}, y::Ptr{Cvoid}, rnd::mpfr_rnd_t)::SPEX_info
end

function SPEX_mpfr_mul(x, y, z, rnd)
    @ccall libspexutil.SPEX_mpfr_mul(x::Ptr{Cvoid}, y::Ptr{Cvoid}, z::Ptr{Cvoid}, rnd::mpfr_rnd_t)::SPEX_info
end

function SPEX_mpfr_mul_d(x, y, z, rnd)
    @ccall libspexutil.SPEX_mpfr_mul_d(x::Ptr{Cvoid}, y::Ptr{Cvoid}, z::Cdouble, rnd::mpfr_rnd_t)::SPEX_info
end

function SPEX_mpfr_div_d(x, y, z, rnd)
    @ccall libspexutil.SPEX_mpfr_div_d(x::Ptr{Cvoid}, y::Ptr{Cvoid}, z::Cdouble, rnd::mpfr_rnd_t)::SPEX_info
end

function SPEX_mpfr_ui_pow_ui(x, y, z, rnd)
    @ccall libspexutil.SPEX_mpfr_ui_pow_ui(x::Ptr{Cvoid}, y::UInt64, z::UInt64, rnd::mpfr_rnd_t)::SPEX_info
end

function SPEX_mpfr_sgn(sgn, x)
    @ccall libspexutil.SPEX_mpfr_sgn(sgn::Ptr{Cint}, x::Ptr{Cvoid})::SPEX_info
end

function SPEX_mpfr_free_cache()
    @ccall libspexutil.SPEX_mpfr_free_cache()::SPEX_info
end

function SPEX_mpfr_free_str(str)
    @ccall libspexutil.SPEX_mpfr_free_str(str::Ptr{Cchar})::SPEX_info
end

function SPEX_Left_LU_backslash(X_handle, type, A, b, option)
    @ccall libspexleftlu.SPEX_Left_LU_backslash(X_handle::Ptr{Ptr{SPEX_matrix}}, type::SPEX_type, A::Ptr{SPEX_matrix}, b::Ptr{SPEX_matrix}, option::Ptr{SPEX_options})::SPEX_info
end

function SPEX_Left_LU_factorize(L_handle, U_handle, rhos_handle, pinv_handle, A, S, option)
    @ccall libspexleftlu.SPEX_Left_LU_factorize(L_handle::Ptr{Ptr{SPEX_matrix}}, U_handle::Ptr{Ptr{SPEX_matrix}}, rhos_handle::Ptr{Ptr{SPEX_matrix}}, pinv_handle::Ptr{Ptr{Int64}}, A::Ptr{SPEX_matrix}, S::Ptr{SPEX_LU_analysis}, option::Ptr{SPEX_options})::SPEX_info
end

function SPEX_Left_LU_solve(X_handle, b, A, L, U, rhos, S, pinv, option)
    @ccall libspexleftlu.SPEX_Left_LU_solve(X_handle::Ptr{Ptr{SPEX_matrix}}, b::Ptr{SPEX_matrix}, A::Ptr{SPEX_matrix}, L::Ptr{SPEX_matrix}, U::Ptr{SPEX_matrix}, rhos::Ptr{SPEX_matrix}, S::Ptr{SPEX_LU_analysis}, pinv::Ptr{Int64}, option::Ptr{SPEX_options})::SPEX_info
end

const SuiteSparse_long = Clong

const SuiteSparse_long_max = typemax(Clong)

const SuiteSparse_long_idd = "ld"

const SUITESPARSE_DATE = "Mar 3, 2021"

const SUITESPARSE_MAIN_VERSION = 5

const SUITESPARSE_SUB_VERSION = 9

const SUITESPARSE_SUBSUB_VERSION = 0

const SPEX_UTIL_VERSION = "1.1.1"

const SPEX_UTIL_VERSION_MAJOR = 1

const SPEX_UTIL_VERSION_MINOR = 1

const SPEX_UTIL_VERSION_SUB = 1

const SPEX_LEFT_LU_VERSION = "1.1.1"

const SPEX_LEFT_LU_VERSION_MAJOR = 1

const SPEX_LEFT_LU_VERSION_MINOR = 1

const SPEX_LEFT_LU_VERSION_SUB = 1

