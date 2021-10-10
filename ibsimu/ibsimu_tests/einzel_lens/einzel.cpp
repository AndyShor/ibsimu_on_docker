/* Cylindrically symmetric simulation of an einzel lens
 *
 * A 10 keV, 1 mA Ar-8+ beam propagates through a cylindrically
 * symmetric, accelerating einzel lens with 5 kV voltage.
 *
 * This is the first example of the Jyvaskyla Summer School course
 * Computational Ion Optics with IBSimu (PH1).
 */

#include <fstream>
#include <iostream>
#include <ibsimu.hpp>
#include <error.hpp>
#include <geometry.hpp>
#include <func_solid.hpp>
#include <epot_field.hpp>
#include <epot_efield.hpp>
#include <meshvectorfield.hpp>
#include <epot_bicgstabsolver.hpp>
#include <gtkplotter.hpp>
#include <trajectorydiagnostics.hpp>



using namespace std;


// Declaration of global parameters
double Veinzel = -5e3;
double h = 1e-3;
double nperh = 100;
double I = 1e-3;
double r0 = 10e-3;
double Npart = r0/h*nperh;
double J = I/(M_PI*r0*r0);
double E0 = 10e3;
double Tt = 0.1;


// Function defining the electrode 1 (ground electrodes of the einzel)
bool solid1( double x, double r, double z )
{
    return( fabs(x) >= 35e-3 && r >= 25e-3 );
}


// Function defining the electrode 2 (biased electrode)
bool solid2( double x, double r, double z )
{
    return( fabs(x) <= 25e-3 && r >= 25e-3 );
}


// Simulation function
void simu( int *argc, char ***argv )
{
    // Create simulation box
    Vec3D origin( -100e-3, 0, 0 );
    Vec3D sizereq( 200e-3, 30e-3, 0 );
    Int3D size( floor(sizereq[0]/h)+1,
		floor(sizereq[1]/h)+1,
		1 );
    Geometry geom( MODE_CYL, size, origin, h );

    // Define solids
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );

    // Set boundary conditions
    geom.set_boundary( 1, Bound(BOUND_NEUMANN, 0.0) ); // xmin
    geom.set_boundary( 2, Bound(BOUND_NEUMANN, 0.0) ); // xmax
    geom.set_boundary( 3, Bound(BOUND_NEUMANN, 0.0) ); // rmin
    geom.set_boundary( 4, Bound(BOUND_NEUMANN, 0.0) ); // rmax
    
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, Veinzel) );

    geom.build_mesh();
    
    // Construct the necessary fields
    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshVectorField bfield;
    EpotEfield efield( epot );
    field_extrpl_e extrapl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
				  FIELD_SYMMETRIC_POTENTIAL, FIELD_EXTRAPOLATE,
				  FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( extrapl );

    // Construct solver and particle database
    EpotBiCGSTABSolver solver( geom );
    ParticleDataBaseCyl pdb( geom );
    
    // Loop for running the vlasov iteration
    for( int a = 0; a < 15; a++ ) {
	
	ibsimu.message(1) << "Major cycle " << a << "\n";
	ibsimu.message(1) << "-----------------------\n";
	
	// Solve the Poisson
	solver.solve( epot, scharge );
	if( solver.get_iter() == 0 ) {
	    ibsimu.message(1) << "No iterations, breaking major cycle\n";
	    break;
	}	
	efield.recalculate();
	
	// Define beam
	pdb.clear();
	pdb.add_2d_beam_with_energy( Npart, J, 8, 40, E0, 0, Tt, 
				     geom.origo(0), 0,
				     geom.origo(0), r0 );
	pdb.iterate_trajectories( scharge, efield, bfield );
	
	// Make diagnostics on beam
	TrajectoryDiagnosticData tdata;
	vector<trajectory_diagnostic_e> diag;
	diag.push_back( DIAG_R );
	diag.push_back( DIAG_RP );
	diag.push_back( DIAG_AP );
	diag.push_back( DIAG_CURR );
	pdb.trajectories_at_plane( tdata, AXIS_X, geom.max(0), diag );
	//Emittance emit( tdata(0).data(), tdata(1).data() );
	EmittanceConv emit( 100, 100, tdata(0).data(), tdata(1).data(),
			    tdata(2).data(), tdata(3).data() );
	
	// Append data to emit.txt for checking for convergence
	ofstream dataout( "emit.txt", ios_base::app );
	dataout << emit.alpha() << " "
		<< emit.beta() << " "
		<< emit.epsilon() << "\n";
	dataout.close();
    }

    // Launch interactive plotter
    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_particledatabase( &pdb );
    plotter.set_efield( &efield );
    plotter.set_bfield( &bfield );
    plotter.set_scharge( &scharge );
    plotter.new_geometry_plot_window();
    plotter.run();
}


// Main function of the program
int main( int argc, char **argv )
{
    // Remove the file emit.txt
    remove( "emit.txt" );
    
    // Structure for catching errors in IBSimu library
    try {
	ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 2 );
	simu( &argc, &argv );
    } catch( Error e ) {
	e.print_error_message( ibsimu.message(0) );
    }
    
    return( 0 );
}
