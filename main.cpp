#include <QCoreApplication>

#include <iostream>//cout
#include <stdlib.h>
#include <fstream> //ofstream
#include <vector>

#include <boost/numeric/odeint.hpp>


typedef boost::array< double , 7 > state_type;

/*=============================================================================================== CONSTANTS =============*/

const double E = 2.15*1e9;
const double E_air = 142*1e6; // Wikipedia, bulk modulus, adiabatic bulk modulus jer je cilindar zadan adijabatski
const double ro = 1000;
const double V_init = 5 * 1e-6;
const double A = 2.5*1e-5;
const double p_gas = 1e5;
const double F_pr = 22.5;
const double k = 800000;
const double c = 1000;
const double masa = 0.078;
const double x_init = 0;
const double x2_init = 1;
const double v_init = 0;

const double n_inFlow = 1;
const double n_cylinders = 1;
const double n_parallel_springs = 1;
double flow=0;

// Flow table values, for testing:
double m_dot(double t)
{
    double protok = 0;
    double on = 0.0001;
    double koef = on/0.05;

    if ( (t >0.2) & (t<=0.25) ) protok = koef*(t-0.2);
    else if ( (t >0.25) & (t<0.3) ) protok = on- koef*(t-0.25);

    if ( (t >0.8) & (t<=0.85) ) protok = -koef*(t-0.8);
    else if ( (t >0.85) & (t<0.9) ) protok = -(on- koef*(t-0.85));

    return protok;
}

// Klasa spomenuta u mailu, ne koristi se
class block {
    double n_volumes; // ili n_cylinders
    double n_inFlow; // ulazni protok
    double n_pistons;

public:
    block(){
        n_volumes = 1;
        n_pistons = n_volumes; // klipova ima koliko i cilindara, tj. volumena hidrauličkog fluida
        n_inFlow = 1;
    }
    double volumes( double pomak ){
        return ( V_init + A * pomak ) * n_volumes;
    }
    double fluid_preassure( double tlak ){
        return ( 1/masa ) * ( A*( tlak ) );
    }
    double inFlow( double t ){
        return m_dot(t) * n_inFlow;
    }
    double springs ( double n_parallel_springs , double pomak , double brzina ){
        return ( ( 1/masa ) * ( F_pr + k*( pomak ) + c*( brzina ) ) ) * n_parallel_springs;
    }

    void set_pistons(double a){
        n_pistons = a;
    }
    void set_volumes(double a){
        n_volumes = a;
        n_pistons = n_volumes;
    }
    void set_flow(double a){
        n_inFlow = a;
    }
    void operator() ( const state_type &x , state_type &dxdt , const double t )
    {
        dxdt[0] = ( E / ( ro * block().volumes( x[1] ) ) ) * ( block().inFlow( t ) +  A*ro*x[2] ); // tlak u hidrauličkom fluidu
        dxdt[3] = ( E_air / ( ro * block().volumes( x[1] ) ) ) * (( A*ro*x[2] ) );                  // tlak u zraku
        dxdt[1] = x[2];                                                                         // dx  /dt
        dxdt[2] = block().fluid_preassure( x[0] ) - block().fluid_preassure( x[3] ) - springs( 1 , x[1], x[2] ) ;           // d^2x/dt^2
    }
};

// Funckije za apstrakciju proračuna
double springs ( double n_parallel_springs , double pomak , double brzina )
{
    return ( ( 1/masa ) * ( F_pr + k*( pomak ) + c*( brzina ) ) ) * n_parallel_springs;
}

double volumes ( double n_cylinders , double pomak )
{
    return ( V_init + A * pomak ) * n_cylinders;
}

double fluid ( double tlak )
{
    return ( 1/masa ) * ( A*( tlak ) );
}

double inFlow ( double n_inFlow , double t )
{
    return m_dot(t) * n_inFlow;
}

// Funkc for integration: x= 0-tlak 1-pomak 2-brzina
void RHS2( const state_type &x , state_type &dxdt , const double t )
{
    flow = inFlow( n_inFlow , t );
    dxdt[0] = ( E / ( ro * volumes( n_cylinders , x[1] ) ) ) * ( flow + ( A*ro*x[2] ) );        // tlak u hidrauličkom fluidu
    dxdt[3] = ( E_air / ( ro * volumes( n_cylinders , x[1] ) ) ) * (( A*ro*x[2] ) );            // tlak u zraku
    dxdt[1] = x[2];                                                                             // dx  /dt
    dxdt[2] = fluid( x[0] ) - fluid( x[3] ) - springs( n_parallel_springs , x[1], x[2] ) ;      // d^2x/dt^2
}

void RHSserialSpring( const state_type &x , state_type &dxdt , const double t )
{
    dxdt[0] = ( E/( ro*( V_init + A*x[5] ) ) )*( m_dot(t) + ( A*ro*x[6] ) );                                // dp  /dt

    // coupled mass spring sustav:
    dxdt[1] = x[2] ;
    dxdt[2] = - k*( x[1] ) + k*( x[3] - x[1] - x2_init ) - c*( x[2] ) ;
    dxdt[3] = x[4] ;
    dxdt[4] =  -k*( x[3] - x[1] - x2_init ) - c*( x[4] ) ;

    // pomak klipa, dxdt[2] i dxdt[4] biti će u N, tako da ulaze ko pribrojnici u jed.
    dxdt[5] = x[6];                                                                                         // dx  /dt
    dxdt[6] = ( 1/masa )*  ( A*( x[0] - p_gas ) - F_pr + dxdt[2] + dxdt[4] ) ;                              // d^2x/dt^2

}

void RHS( const state_type &x , state_type &dxdt , const double t )
{
    dxdt[0] = ( E/( ro*( V_init + A*x[1] ) ) )*( m_dot(t) + ( A*ro*x[2] ) );                                // dp  /dt
    dxdt[1] = x[2];                                                                                         // dx  /dt
    dxdt[2] = ( 1/masa )*  ( A*( x[0] - p_gas ) - F_pr - k*( x[1] ) - c*( x[2] ) ) ;                        // d^2x/dt^2
}

void write_console(const state_type &x, const double t)
{
    std::cout << t << "\t m: " << flow << "\t p: " << x[0] << "\t x: " << x[1] << "\t v: " << x[2] << std::endl;
}

struct streaming_observer
{
    std::ostream &m_out;
    streaming_observer(std::ostream &out) : m_out(out) {}

    void operator()(const state_type &x, double t) const
    {
        //m_out << t << "\t m: " << m_dot(t) << "\t p: " << x[0] << "\t x: " << x[1] << "\t v: " << x[2] << "\n";
        m_out << t << " " << flow << " " << x[0] << " " << x[5] << " " << x[6] << "\n";
    }
};

/*=============================================================================================== MAIN ==================*/
int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    using namespace std;
    using namespace boost::numeric::odeint;

    // output za provjeru
    ofstream out("zadatakData.txt");
    if (!out.is_open())
    {
        cout << "Output file not opened.\n";
        system("pause");
        return 1;
    }

    state_type x={1e6, 0.0, 0.0}; //Initial values: x= 0-tlak 1-pomak 2-brzina
    double t_start = 0.0;
    double t_end = 1.0;
    double dt = 1e-6;
    double steps = 0; //integrate returns

    //runge_kutta4< state_type > stepper;
    //integrate_const(stepper, RHS2, x, t_start, t_end, dt, streaming_observer(out));

    steps = integrate(RHS2, x, t_start, t_end, dt, streaming_observer(out));  // sys( x , dxdt , t )
    cout << "Number of integration steps: " << steps << "\n";

    out.close();

    system("./drawGraphs.py");

    return 0; //a.exec();
}

