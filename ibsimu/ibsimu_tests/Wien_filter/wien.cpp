#include <fstream>
#include <iomanip>
#include <limits>
#include "epot_bicgstabsolver.hpp"
#include "meshvectorfield.hpp"
#include "mydxffile.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "dxf_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"


using namespace std;


double Vplate = 2e3;
double E0 = 10e3;

double q = 1.0;
double m = 1.0;
double Tt = 0.01;
double r0 = 3e-3;
double I = 100e-6;
double J = I/(r0*r0*M_PI);


bool gnd( double x, double y, double z )
{
    double r = sqrt(x*x+y*y);
    return( (z <= 0 || z >= 130e-3) && r > 3e-3 );
}


bool plate1( double x, double y, double z )
{
    return( fabs(y) <= 15e-3 && x >= 15e-3 && x <= 20e-3 && z > 15e-3 && z < 115e-3 );
}


bool plate2( double x, double y, double z )
{
    return( fabs(y) <= 15e-3 && x <= -15e-3 && x >= -20e-3 && z > 15e-3 && z < 115e-3 );
}


void simu( int *argc, char ***argv )
{
    double h = 1e-3;
    Vec3D origo( -25e-3, 
		 -25e-3, 
		 -10e-3 );
    double sizereq[3] = {  50e-3,
                           50e-3, 
                          150e-3 };
    Int3D meshsize( (int)floor(sizereq[0]/h)+1,
                    (int)floor(sizereq[1]/h)+1,
		    (int)floor(sizereq[2]/h)+1 );
    Geometry geom( MODE_3D, meshsize, origo, h );

    Solid *s1 = new FuncSolid( gnd );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( plate1 );
    geom.set_solid( 8, s2 );
    Solid *s3 = new FuncSolid( plate2 );
    geom.set_solid( 9, s3 );

    geom.set_boundary( 1, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 3, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 5, Bound(BOUND_NEUMANN, 0.0) );
    geom.set_boundary( 6, Bound(BOUND_NEUMANN, 0.0) );

    geom.set_boundary( 7, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, -Vplate) );
    geom.set_boundary( 9, Bound(BOUND_DIRICHLET, Vplate) );

    geom.build_mesh();
    geom.build_surface();

    EpotBiCGSTABSolver solver( geom );
    EpotField epot( geom );
    MeshScalarField scharge( geom );
    bool bfield_fout[3] = {false, true, false};
    MeshVectorField bfield( geom, bfield_fout );
    const double Bmax = 0.1;
    for( uint32_t k = 0; k < bfield.size(2); k++ ) {
	double z = bfield.origo(2)+k*bfield.h();
	double By;
	if( z < 0 )
	    By = 0.0;
	else if( z < 15e-3 )
	    By = z/15e-3*Bmax;
	else if( z < 115e-3 )
	    By = Bmax;
	else if( z < 130e-3 )
	    By = (130e-3-z)/15e-3*Bmax;
	else
	    By = 0.0;

	for( uint32_t j = 0; j < bfield.size(1); j++ ) {
	    for( uint32_t i = 0; i < bfield.size(0); i++ ) {
		bfield.set(i,j,k,Vec3D(0,By,0));
	    }
	}
    }

    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                                     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBase3D pdb( geom );
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );

    for( size_t i = 0; i < 1; i++ ) {

	ibsimu.message(1) << "Major cycle " << i << "\n";
	ibsimu.message(1) << "-----------------------\n";

	solver.solve( epot, scharge );
	if( solver.get_iter() == 0 ) {
	    ibsimu.message(1) << "No iterations, breaking major cycle\n";
	    break;
	}
	efield.recalculate();

        pdb.clear(); 
	int Npart = 5000;
	pdb.add_cylindrical_beam_with_energy( 0.8*Npart, 0.8*J, q, 1,
					      E0, 0, Tt,
					      Vec3D(0,0,geom.origo(2)),
					      Vec3D(1,0,0),
					      Vec3D(0,1,0),
					      r0 );
	pdb.add_cylindrical_beam_with_energy( 0.1*Npart, 0.1*J, q, 2,
					      E0, 0, Tt,
					      Vec3D(0,0,geom.origo(2)),
					      Vec3D(1,0,0),
					      Vec3D(0,1,0),
					      r0 );
	pdb.add_cylindrical_beam_with_energy( 0.1*Npart, 0.1*J, q, 3,
					      E0, 0, Tt,
					      Vec3D(0,0,geom.origo(2)),
					      Vec3D(1,0,0),
					      Vec3D(0,1,0),
					      r0 );
	
        pdb.iterate_trajectories( scharge, efield, bfield );

	/*
	TrajectoryDiagnosticData tdata;
        std::vector<trajectory_diagnostic_e> diagnostics;
        diagnostics.push_back( DIAG_X );
        diagnostics.push_back( DIAG_XP );
        pdb.trajectories_at_plane( tdata, AXIS_Z, geom.max(2)-geom.h(), diagnostics );
        Emittance emit( tdata(0).data(), tdata(1).data() );

        // Output
        ofstream dout( "emit.txt", ios_base::app );
        dout << emit.alpha() << " "
             << emit.beta() << " "
             << emit.epsilon() << "\n";
        dout.close();
	*/
    }
    
    if( false ) {
	MeshScalarField tdens(geom);
	pdb.build_trajectory_density_field(tdens);

	GTKPlotter plotter( argc, argv );
	plotter.set_geometry( &geom );
	plotter.set_epot( &epot );
	plotter.set_bfield( &bfield );
	plotter.set_trajdens( &tdens );
	plotter.set_scharge( &scharge );
	plotter.set_particledatabase( &pdb );
	plotter.new_geometry_plot_window();
	plotter.run();
    }

    GeomPlotter plotter( geom );
    plotter.set_epot( &epot );
    plotter.set_particledatabase( &pdb );
    plotter.set_view( VIEW_ZX );
    stringstream ss;
    ss << "plot_" << Vplate << ".png";
    plotter.plot_png( ss.str() );

    ofstream dout( "diag.txt", ios_base::app );
    ParticleStatistics stat = pdb.get_statistics();
    dout << Vplate << " " << stat.bound_current(6) << "\n";
    dout.close();
}


int main( int argc, char **argv )
{
    remove( "diag.txt" );

    try {
	//ibsimu.set_message_output( "ibsimu.txt" );
        ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );

	for( Vplate = 1000; Vplate < 3000; Vplate += 100 )
	    simu( &argc, &argv );


    } catch( Error e ) {
	e.print_error_message( ibsimu.message( 0 ) );
        exit( 1 );
    }

    return( 0 );
}
