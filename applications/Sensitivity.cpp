/****************************************************************************
 * Parameter sensitivity of the Finch heat solve, by forward-mode AD.
 *
 * Runs the single-layer heat transport problem once with OTI-valued material
 * and source parameters, obtaining d(QoI)/dp for every parameter from that one
 * solve, then validates each derivative against a central finite difference
 * built from two ordinary double solves per parameter.
 ****************************************************************************/

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <mpi.h>

#include <Kokkos_Core.hpp>

#include "Finch_Core.hpp"
#include "Finch_OTI.hpp"

namespace
{

constexpr int NP = Finch::Sensitivity::NumParameters;

// First-order jets: one variable per differentiated parameter, derivatives
// through order one. Raising the order to 2 here is the only change needed to
// obtain the full parameter Hessian of the same quantities.
using OTI = oti::otinum<NP, 1>;

// Quantities of interest evaluated on the final temperature field. Both are
// smooth functionals of the field, which makes them meaningful finite-
// difference references.
template <class Scalar>
struct QoI
{
    // Sum of temperature over all owned nodes.
    Scalar sum;
    // Temperature at a single fixed node (centre of the owned index space).
    Scalar probe;
};

// Run the transient solve to completion and evaluate the QoIs.
template <class Scalar, class MemorySpace, class ExecSpace>
QoI<Scalar> solve( MPI_Comm comm, Finch::Inputs db,
                   const Finch::MaterialProperties<Scalar>& props )
{
    std::array<std::string, 6> bc_types = { "adiabatic", "adiabatic",
                                            "adiabatic", "adiabatic",
                                            "adiabatic", "adiabatic" };

    Finch::Grid<MemorySpace, Scalar> grid(
        comm, db.space.cell_size, db.space.global_low_corner,
        db.space.global_high_corner, db.space.ranks_per_dim, bc_types,
        db.space.initial_temperature );

    auto fd = Finch::createSolver( db, grid, props );

    Finch::MovingBeam beam( db.source.scan_path_file );

    // Time loop. This mirrors Finch::Layer::step without the solidification
    // sampling and file output, which are not needed for the QoIs here.
    double time = db.time.start_time;
    const double dt = db.time.time_step;

    for ( int n = 0; n < db.time.num_steps; ++n )
    {
        time += dt;

        beam.move( time );
        double beam_power = beam.power();
        double beam_pos[3];
        for ( std::size_t d = 0; d < 3; ++d )
            beam_pos[d] = beam.position( d );

        auto T = grid.getTemperature();
        auto T0 = grid.getPreviousTemperature();

        Kokkos::deep_copy( T0, T );

        auto owned_space = grid.getIndexSpace();
        fd.solve( ExecSpace(), owned_space, T, T0, beam_power, beam_pos );

        grid.updateBoundaries();
        grid.gather();
    }

    // Evaluate the QoIs on the host. Copying the field to a host mirror keeps
    // this independent of Kokkos reducer support for compound scalar types,
    // which is not needed for a once-per-run diagnostic.
    auto T = grid.getTemperature();
    auto T_host = Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), T );

    auto owned = grid.getIndexSpace();

    QoI<Scalar> q;
    q.sum = Scalar( 0 );
    for ( int i = owned.min( 0 ); i < owned.max( 0 ); ++i )
        for ( int j = owned.min( 1 ); j < owned.max( 1 ); ++j )
            for ( int k = owned.min( 2 ); k < owned.max( 2 ); ++k )
                q.sum = q.sum + T_host( i, j, k, 0 );

    q.probe = T_host( ( owned.min( 0 ) + owned.max( 0 ) ) / 2,
                      ( owned.min( 1 ) + owned.max( 1 ) ) / 2,
                      ( owned.min( 2 ) + owned.max( 2 ) ) / 2, 0 );

    // Reduce the sum across ranks. An OTI jet is a contiguous block of
    // coefficients, and summing jets is summing coefficients elementwise, so a
    // plain MPI_DOUBLE reduction over ncoeffs entries is exact -- no derived
    // datatype or user-defined operator is required.
    int comm_size;
    MPI_Comm_size( comm, &comm_size );
    if ( comm_size > 1 )
    {
        if constexpr ( std::is_same<Scalar, double>::value )
        {
            double local = q.sum;
            MPI_Allreduce( &local, &q.sum, 1, MPI_DOUBLE, MPI_SUM, comm );
        }
        else
        {
            Scalar local = q.sum;
            MPI_Allreduce( &local[0], &q.sum[0], Scalar::ncoeffs, MPI_DOUBLE,
                           MPI_SUM, comm );
        }
    }

    return q;
}

void run( MPI_Comm comm, int argc, char* argv[] )
{
    using exec_space = Kokkos::DefaultExecutionSpace;
    using memory_space = exec_space::memory_space;

    int rank;
    MPI_Comm_rank( comm, &rank );

    Finch::Inputs db( comm, argc, argv );

    auto nominal = Finch::Sensitivity::nominal( db );

    // FD step is overridable so the derivative can be checked for step-size
    // independence -- the standard way to tell a genuine discrepancy from a
    // finite-difference artifact.
    double rel_step = 1e-6;
    if ( const char* s = std::getenv( "FINCH_FD_STEP" ) )
        rel_step = std::atof( s );

    // Warm up: the first solve of a process pays for page faults and OpenMP
    // thread start-up, which would otherwise be charged to whichever solve
    // happens to run first.
    {
        std::array<double, NP> warm;
        for ( int i = 0; i < NP; ++i )
            warm[i] = nominal[i];
        Finch::Inputs warm_db = db;
        warm_db.time.num_steps = 2;
        solve<double, memory_space, exec_space>(
            comm, warm_db, Finch::Sensitivity::build<double>( warm ) );
    }

    // ---- One OTI solve gives every derivative -----------------------------
    double t0 = MPI_Wtime();
    auto seeded = Finch::Sensitivity::seed<OTI>( db );
    auto oti_result = solve<OTI, memory_space, exec_space>(
        comm, db, Finch::Sensitivity::build<OTI>( seeded ) );
    double t_oti = MPI_Wtime() - t0;

    // ---- Reference: one plain double solve, then two per parameter --------
    t0 = MPI_Wtime();
    std::array<double, NP> base;
    for ( int i = 0; i < NP; ++i )
        base[i] = nominal[i];
    auto ref = solve<double, memory_space, exec_space>(
        comm, db, Finch::Sensitivity::build<double>( base ) );
    double t_double = MPI_Wtime() - t0;

    std::array<double, NP> fd_sum, fd_probe;
    t0 = MPI_Wtime();
    for ( int i = 0; i < NP; ++i )
    {
        double h = rel_step * std::fabs( nominal[i] );

        auto plus = base;
        plus[i] = nominal[i] + h;
        auto qp = solve<double, memory_space, exec_space>(
            comm, db, Finch::Sensitivity::build<double>( plus ) );

        auto minus = base;
        minus[i] = nominal[i] - h;
        auto qm = solve<double, memory_space, exec_space>(
            comm, db, Finch::Sensitivity::build<double>( minus ) );

        fd_sum[i] = ( qp.sum - qm.sum ) / ( 2.0 * h );
        fd_probe[i] = ( qp.probe - qm.probe ) / ( 2.0 * h );
    }
    double t_fd = MPI_Wtime() - t0;

    if ( rank != 0 )
        return;

    // ---- Report -----------------------------------------------------------
    std::printf( "\n" );
    std::printf( "================================================="
                 "=================================\n" );
    std::printf( "  Finch parameter sensitivity: OTI forward-mode AD vs "
                 "central finite differences\n" );
    std::printf( "================================================="
                 "=================================\n" );
    std::printf( "  algebra          : oti::otinum<%d,%d>  (%d coefficients "
                 "per node)\n",
                 OTI::nvars, OTI::order, OTI::ncoeffs );
    std::printf( "  time steps       : %d\n", db.time.num_steps );
    std::printf( "  FD relative step : %g\n", rel_step );
    std::printf( "\n" );
    std::printf( "  value check   T_sum   OTI %.10e   double %.10e\n",
                 oti_result.sum.real(), ref.sum );
    std::printf( "  value check   T_probe OTI %.10e   double %.10e\n",
                 oti_result.probe.real(), ref.probe );
    std::printf( "\n" );

    const char* qoi_name[2] = { "T_sum  [K]", "T_probe[K]" };

    for ( int q = 0; q < 2; ++q )
    {
        const OTI& jet = ( q == 0 ) ? oti_result.sum : oti_result.probe;
        const std::array<double, NP>& fd = ( q == 0 ) ? fd_sum : fd_probe;

        std::printf( "  d(%s)/dp\n", qoi_name[q] );
        std::printf( "  %-22s %16s %16s %11s %14s\n", "parameter", "OTI",
                     "central FD", "rel.diff", "p*dQ/dp [K]" );
        std::printf( "  %-22s %16s %16s %11s %14s\n",
                     "----------------------", "----------------",
                     "----------------", "-----------", "--------------" );

        // A relative difference is only meaningful if the derivative is itself
        // non-negligible. Compare each row against the largest normalized
        // sensitivity in the table so that a structurally-zero derivative is
        // reported as such rather than as a 100% disagreement.
        double ref = 0.0;
        for ( int i = 0; i < NP; ++i )
        {
            typename OTI::alpha_type alpha{};
            alpha[i] = 1;
            ref = std::max( ref, std::fabs( nominal[i] * jet.partial( alpha ) ) );
        }

        for ( int i = 0; i < NP; ++i )
        {
            typename OTI::alpha_type alpha{};
            alpha[i] = 1;
            double d_oti = jet.partial( alpha );

            double scale = std::max( std::fabs( d_oti ), std::fabs( fd[i] ) );
            bool negligible =
                nominal[i] != 0.0 &&
                std::fabs( nominal[i] ) * scale < 1e-8 * ref;

            char label[64];
            std::snprintf( label, sizeof( label ), "%s [%s]",
                           Finch::Sensitivity::name( i ),
                           Finch::Sensitivity::units( i ) );

            std::printf( "  %-22s %16.6e %16.6e ", label, d_oti, fd[i] );
            if ( negligible )
                std::printf( "%11s ", "negligible" );
            else if ( scale > 0.0 )
                std::printf( "%11.2e ",
                             std::fabs( d_oti - fd[i] ) / scale );
            else
                std::printf( "%11s ", "-" );
            std::printf( "%14.4e\n", nominal[i] * d_oti );
        }
        std::printf( "\n" );
    }

    std::printf( "  cost:  1 OTI solve %.3f s   |   1 double solve %.3f s"
                 "   |   %d FD solves %.3f s\n",
                 t_oti, t_double, 2 * NP, t_fd );
    std::printf( "  OTI overhead vs one double solve: %.2fx  "
                 "(all %d derivatives)\n",
                 t_oti / t_double, NP );
    std::printf( "  FD cost for the same %d derivatives: %.2fx\n", NP,
                 t_fd / t_double );
    std::printf( "================================================="
                 "=================================\n\n" );
}

} // namespace

int main( int argc, char* argv[] )
{
    MPI_Init( &argc, &argv );
    Kokkos::initialize( argc, argv );

    run( MPI_COMM_WORLD, argc, argv );

    Kokkos::finalize();
    MPI_Finalize();

    return 0;
}
