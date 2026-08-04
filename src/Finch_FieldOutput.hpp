/****************************************************************************
 * Copyright (c) 2024 by Oak Ridge National Laboratory                      *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of Finch. Finch is distributed under a                 *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#ifndef FINCH_FIELD_OUTPUT_HPP
#define FINCH_FIELD_OUTPUT_HPP

#include <array>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>

#include <Finch_Core_Config.hpp>
#include <Finch_Grid.hpp>
#include <Finch_Inputs.hpp>

#if Finch_ENABLE_ADIOS2
#include <adios2.h>
#include <adios2/cxx/KokkosView.h>
#endif

namespace Finch
{

#if Finch_ENABLE_ADIOS2
template <typename MemorySpace>
class Adios2FieldOutput
{
  public:
    using memory_space = MemorySpace;
    using grid_type = Grid<memory_space>;
    using execution_space = typename grid_type::exec_space;
    using temperature_view_type = typename grid_type::view_type;
    using temperature_layout = typename temperature_view_type::array_layout;
    using packed_view_type =
        Kokkos::View<double***, Kokkos::LayoutLeft, memory_space>;
    using derived_view_type =
        Kokkos::View<float***, Kokkos::LayoutLeft, memory_space>;
    using direct_view_type =
        Kokkos::View<double***, Kokkos::LayoutLeft, memory_space,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    static constexpr bool direct_layout =
        std::is_same<temperature_layout, Kokkos::LayoutLeft>::value;

    Adios2FieldOutput( MPI_Comm comm, grid_type& grid,
                       const FieldOutputInput& input )
        : adios_( comm )
        , io_( adios_.DeclareIO( "FinchFieldOutput" ) )
    {
        Kokkos::Profiling::ScopedRegion region( "Finch::field_output_open" );
        std::vector<std::string> field_names;
        field_names.reserve( input.fields.size() );
        for ( const auto field : input.fields )
        {
            field_names.emplace_back( fieldOutputFieldName( field ) );
            write_temperature_ |= field == FieldOutputField::temperature;
            write_volumetric_heat_source_ |=
                field == FieldOutputField::volumetric_heat_source;
        }

        const auto owned = grid.getIndexSpace();
        const auto ghosted = grid.getGhostedIndexSpace();
        const auto& global_grid = grid.getLocalGrid()->globalGrid();
        const auto& mesh = global_grid.globalMesh();

        adios2::Dims shape( 3 );
        adios2::Dims start( 3 );
        adios2::Dims count( 3 );
        adios2::Dims memory_start( 3 );
        adios2::Dims memory_count( 3 );
        for ( int d = 0; d < 3; ++d )
        {
            shape[d] = static_cast<std::size_t>(
                           global_grid.globalNumEntity(
                               Cabana::Grid::Cell(), d ) ) +
                       1;
            start[d] = static_cast<std::size_t>(
                global_grid.globalOffset( d ) );
            output_min_[d] = owned.min( d );
            output_max_[d] = owned.max( d ) +
                             ( global_grid.onHighBoundary( d ) ? 0 : 1 );
            count[d] = static_cast<std::size_t>( output_max_[d] -
                                                  output_min_[d] );
            memory_start[d] = static_cast<std::size_t>(
                output_min_[d] - ghosted.min( d ) );
            memory_count[d] =
                static_cast<std::size_t>( ghosted.extent( d ) );
        }

        // Finch's logical dimensions are x/y/z. LayoutLeft makes x
        // contiguous, while ADIOS2 records the standard z/y/x row-major
        // representation used by VTK image data and Fides.
        if ( write_temperature_ )
        {
            temperature_ = io_.DefineVariable<double>( "temperature", shape,
                                                       start, count );
            if constexpr ( direct_layout )
                temperature_.SetMemorySelection(
                    { memory_start, memory_count } );
            temperature_.SetArrayLayout( adios2::ArrayOrdering::ColumnMajor );
        }
        if ( write_volumetric_heat_source_ )
        {
            volumetric_heat_source_ = io_.DefineVariable<float>(
                "volumetric_heat_source", shape, start, count );
            volumetric_heat_source_.SetArrayLayout(
                adios2::ArrayOrdering::ColumnMajor );
        }

        step_ = io_.DefineVariable<int>( "step" );
        time_ = io_.DefineVariable<double>( "time" );

        const double origin[3] = { mesh.lowCorner( 0 ), mesh.lowCorner( 1 ),
                                   mesh.lowCorner( 2 ) };
        const double spacing[3] = { mesh.cellSize( 0 ), mesh.cellSize( 1 ),
                                    mesh.cellSize( 2 ) };
        io_.DefineAttribute<std::string>( "Fides_Data_Model", "uniform" );
        io_.DefineAttribute<double>( "Fides_Origin", origin, 3 );
        io_.DefineAttribute<double>( "Fides_Spacing", spacing, 3 );
        io_.DefineAttribute<std::string>( "Fides_Dimension_Variable",
                                          field_names.front() );
        const std::vector<std::string> associations( field_names.size(),
                                                     "points" );
        io_.DefineAttribute<std::string>( "Fides_Variable_List",
                                          field_names.data(),
                                          field_names.size() );
        io_.DefineAttribute<std::string>(
            "Fides_Variable_Associations", associations.data(),
            associations.size() );
        if ( write_volumetric_heat_source_ )
            io_.DefineAttribute<std::string>( "volumetric_heat_source_units",
                                              "W/m^3" );
        io_.DefineAttribute<std::string>(
            "vtk.xml",
            vtkImageSchema( shape, origin, spacing, field_names ) );

        if constexpr ( !direct_layout )
        {
            if ( write_temperature_ )
                packed_ = packed_view_type( "Finch field output buffer",
                                            count[0], count[1], count[2] );
        }
        if ( write_volumetric_heat_source_ )
            volumetric_heat_source_buffer_ = derived_view_type(
                "Finch volumetric heat-source output buffer", count[0],
                count[1], count[2] );

        io_.SetEngine( "BP5" );
        engine_ = io_.Open( "fields.bp", adios2::Mode::Write );
        open_ = true;
    }

    Adios2FieldOutput( const Adios2FieldOutput& ) = delete;
    Adios2FieldOutput& operator=( const Adios2FieldOutput& ) = delete;

    ~Adios2FieldOutput()
    {
        if ( open_ )
        {
            try
            {
                engine_.Close();
            }
            catch ( ... )
            {
            }
        }
    }

    template <typename SolverType, typename BeamType>
    void write( grid_type& grid, const SolverType& solver,
                const BeamType& beam, const int step, const double time )
    {
        if ( write_temperature_ )
            grid.gatherForFieldOutput();

        engine_.BeginStep();
        engine_.Put( step_, &step, adios2::Mode::Deferred );
        engine_.Put( time_, &time, adios2::Mode::Deferred );

        if constexpr ( !direct_layout )
        {
            if ( write_temperature_ )
            {
                auto temperature = grid.getTemperature();
                const Cabana::Grid::IndexSpace<3> output_space(
                    output_min_, output_max_ );
                const int i0 = output_min_[0];
                const int j0 = output_min_[1];
                const int k0 = output_min_[2];
                auto packed = packed_;
                Kokkos::parallel_for(
                    "Finch::field_output_pack",
                    Cabana::Grid::createExecutionPolicy(
                        output_space, grid.executionSpace() ),
                    KOKKOS_LAMBDA( const int i, const int j, const int k ) {
                        packed( i - i0, j - j0, k - k0 ) =
                            temperature( i, j, k, 0 );
                    } );
            }
        }

        const Cabana::Grid::IndexSpace<3> output_space( output_min_,
                                                        output_max_ );
        if ( write_volumetric_heat_source_ )
            solver.fillVolumetricSourceField(
                grid.executionSpace(), beam, output_space,
                volumetric_heat_source_buffer_ );

        if ( ( write_temperature_ && !direct_layout ) ||
             write_volumetric_heat_source_ )
            grid.executionSpace().fence( "Finch field output preparation" );

        if ( write_temperature_ )
        {
            if constexpr ( direct_layout )
            {
                auto temperature = grid.getTemperature();
                direct_view_type direct( temperature.data(),
                                         temperature.extent( 0 ),
                                         temperature.extent( 1 ),
                                         temperature.extent( 2 ) );
                engine_.Put( temperature_, direct, adios2::Mode::Deferred );
            }
            else
                engine_.Put( temperature_, packed_, adios2::Mode::Deferred );
        }
        if ( write_volumetric_heat_source_ )
            engine_.Put( volumetric_heat_source_,
                         volumetric_heat_source_buffer_,
                         adios2::Mode::Deferred );

        engine_.EndStep();
    }

    void close()
    {
        if ( open_ )
        {
            Kokkos::Profiling::ScopedRegion region(
                "Finch::field_output_close" );
            engine_.Close();
            open_ = false;
        }
    }

  private:
    static std::string vtkImageSchema( const adios2::Dims& logical_shape,
                                       const double* origin,
                                       const double* spacing,
                                       const std::vector<std::string>& fields )
    {
        std::ostringstream xml;
        xml.precision( 17 );
        const auto x_max = logical_shape[0] - 1;
        const auto y_max = logical_shape[1] - 1;
        const auto z_max = logical_shape[2] - 1;
        xml << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"ImageData\" version=\"0.1\" "
               "byte_order=\"LittleEndian\">\n"
            << "  <ImageData WholeExtent=\"0 " << x_max << " 0 " << y_max
            << " 0 " << z_max << "\" Origin=\"" << origin[0] << ' '
            << origin[1] << ' ' << origin[2] << "\" Spacing=\""
            << spacing[0] << ' ' << spacing[1] << ' ' << spacing[2]
            << "\">\n"
            << "    <Piece Extent=\"0 " << x_max << " 0 " << y_max
            << " 0 " << z_max << "\">\n"
            << "      <PointData Scalars=\"" << fields.front() << "\">\n";
        for ( const auto& field : fields )
            xml << "        <DataArray Name=\"" << field << "\"/>\n";
        xml << "        <DataArray Name=\"TIME\">time</DataArray>\n"
            << "      </PointData>\n"
            << "    </Piece>\n"
            << "  </ImageData>\n"
            << "</VTKFile>\n";
        return xml.str();
    }

    adios2::ADIOS adios_;
    adios2::IO io_;
    adios2::Engine engine_;
    adios2::Variable<double> temperature_;
    adios2::Variable<float> volumetric_heat_source_;
    adios2::Variable<int> step_;
    adios2::Variable<double> time_;
    packed_view_type packed_;
    derived_view_type volumetric_heat_source_buffer_;
    std::array<long, 3> output_min_;
    std::array<long, 3> output_max_;
    bool write_temperature_ = false;
    bool write_volumetric_heat_source_ = false;
    bool open_ = false;
};
#endif

template <typename MemorySpace>
class FieldOutput
{
  public:
    using grid_type = Grid<MemorySpace>;

    FieldOutput( const FieldOutputInput& input, grid_type& grid )
        : format_( input.format )
    {
#if Finch_ENABLE_ADIOS2
        if ( format_ == "adios2" )
            adios2_ =
                std::make_unique<Adios2FieldOutput<MemorySpace>>(
                    grid.getComm(), grid, input );
#else
        (void)grid;
#endif
    }

    template <typename SolverType, typename BeamType>
    void write( grid_type& grid, const SolverType& solver,
                const BeamType& beam, const int step, const double time )
    {
#if Finch_ENABLE_ADIOS2
        if ( adios2_ )
        {
            adios2_->write( grid, solver, beam, step, time );
            return;
        }
#endif
        grid.output( step, time );
    }

    void close()
    {
#if Finch_ENABLE_ADIOS2
        if ( adios2_ )
            adios2_->close();
#endif
    }

    const std::string& format() const { return format_; }

  private:
    std::string format_;
#if Finch_ENABLE_ADIOS2
    std::unique_ptr<Adios2FieldOutput<MemorySpace>> adios2_;
#endif
};

} // namespace Finch

#endif
