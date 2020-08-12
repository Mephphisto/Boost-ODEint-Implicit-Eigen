/*
 [auto_generated]
 boost/numeric/odeint/stepper/rosenbrock4_Eigen.hpp

 [begin_description]
 Implementation of the Rosenbrock 4 method for solving stiff ODEs. Note, that a
 controller and a dense-output stepper exist for this method,
 [end_description]

 Copyright 2011-2013 Karsten Ahnert
 Copyright 2011-2012 Mario Mulansky
 Copyright 2012 Christoph Koke

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_rosenbrock4_Eigen_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_rosenbrock4_Eigen_HPP_INCLUDED


#include <boost/numeric/odeint/util/bind.hpp>
#include <boost/numeric/odeint/util/unwrap_reference.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>

#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>

#include <type_traits>
#include <Eigen/Dense>


namespace boost {
namespace numeric {
namespace odeint {


/*
 * ToDo:
 *
 * 2. Interfacing for odeint, check if controlled_error_stepper can be used
 * 3. dense output
 */


template< class Value,class Time_type = decltype (std::abs(std::declval<Value>())),  class Coefficients = default_rosenbrock_coefficients< Value > , class Resizer = initially_resizer >
class rosenbrock4_Eigen
{
public:

    typedef Value value_type;
    typedef Time_type time_type;
    typedef Eigen::VectorX< value_type > state_type;
    typedef state_type deriv_type;
    typedef Eigen::MatrixX< value_type > matrix_type;

    typedef Resizer resizer_type;
    typedef Coefficients rosenbrock_coefficients;
    typedef stepper_tag stepper_category;
    typedef unsigned short order_type;

    typedef state_wrapper< state_type > wrapped_state_type;
    typedef state_wrapper< deriv_type > wrapped_deriv_type;
    typedef state_wrapper< matrix_type > wrapped_matrix_type;

    typedef rosenbrock4_Eigen< Value , Coefficients , Resizer > stepper_type;

    const static order_type stepper_order = rosenbrock_coefficients::stepper_order;
    const static order_type error_order = rosenbrock_coefficients::error_order;

    rosenbrock4_Eigen( void )
    : m_resizer() , m_x_err_resizer() ,
      m_jac() ,
      m_dfdt() , m_dxdt() , m_dxdtnew() ,
      m_g1() , m_g2() , m_g3() , m_g4() , m_g5() ,
      m_cont3() , m_cont4() , m_xtmp() , m_x_err() ,
      m_coef()
    { }


    order_type order() const { return stepper_order; } 

    template< class System >
    void do_step( System system , const state_type &x , time_type t , state_type &xout , time_type dt , state_type &xerr )
    {
        // get the system and jacobi function
        typedef typename odeint::unwrap_reference< System >::type system_type;
        typedef typename odeint::unwrap_reference< typename system_type::first_type >::type deriv_func_type;
        typedef typename odeint::unwrap_reference< typename system_type::second_type >::type jacobi_func_type;
        system_type &sys = system;
        deriv_func_type &deriv_func = sys.first;
        jacobi_func_type &jacobi_func = sys.second;



        const size_t n = x.size();

        m_resizer.adjust_size( x , detail::bind( &stepper_type::template resize_impl<state_type> , detail::ref( *this ) , detail::_1 ) );

        deriv_func( x , m_dxdt , t );
        jacobi_func( x , m_jac , t , m_dfdt );

        m_jac *= -1.0;
        m_jac += 1.0 / m_coef.gamma / dt * state_type::Identity(n);
        // Eigen Full Pivoted LU Factorisation type
        Eigen::FullPivLU<state_type> lu(m_jac);

        for( size_t i=0 ; i<n ; ++i )
            m_g1[i] = m_dxdt[i] + dt * m_coef.d1 * m_dfdt[i];
        // Eigen LU Solve
        m_g1 = lu.solve(m_g1);


        for( size_t i=0 ; i<n ; ++i )
            m_xtmp[i] = x[i] + m_coef.a21 * m_g1[i];
        deriv_func( m_xtmp , m_dxdtnew , t + m_coef.c2 * dt );
        for( size_t i=0 ; i<n ; ++i )
            m_g2[i] = m_dxdtnew[i] + dt * m_coef.d2 * m_dfdt[i] + m_coef.c21 * m_g1[i] / dt;
        m_g2 = lu.solve(m_g2);


        for( size_t i=0 ; i<n ; ++i )
            m_xtmp[i] = x[i] + m_coef.a31 * m_g1[i] + m_coef.a32 * m_g2[i];
        deriv_func( m_xtmp , m_dxdtnew , t + m_coef.c3 * dt );
        for( size_t i=0 ; i<n ; ++i )
            m_g3[i] = m_dxdtnew[i] + dt * m_coef.d3 * m_dfdt[i] + ( m_coef.c31 * m_g1[i] + m_coef.c32 * m_g2[i] ) / dt;
        m_g3 = lu.solve(m_g3);

        for( size_t i=0 ; i<n ; ++i )
            m_xtmp[i] = x[i] + m_coef.a41 * m_g1[i] + m_coef.a42 * m_g2[i] + m_coef.a43 * m_g3[i];
        deriv_func( m_xtmp , m_dxdtnew , t + m_coef.c4 * dt );
        for( size_t i=0 ; i<n ; ++i )
            m_g4[i] = m_dxdtnew[i] + dt * m_coef.d4 * m_dfdt[i] + ( m_coef.c41 * m_g1[i] + m_coef.c42 * m_g2[i] + m_coef.c43 * m_g3[i] ) / dt;
        m_g4 = lu.solve(m_g4);


        for( size_t i=0 ; i<n ; ++i )
            m_xtmp[i] = x[i] + m_coef.a51 * m_g1[i] + m_coef.a52 * m_g2[i] + m_coef.a53 * m_g3[i] + m_coef.a54 * m_g4[i];
        deriv_func( m_xtmp , m_dxdtnew , t + dt );
        for( size_t i=0 ; i<n ; ++i )
            m_g5[i] = m_dxdtnew[i] + ( m_coef.c51 * m_g1[i] + m_coef.c52 * m_g2[i] + m_coef.c53 * m_g3[i] + m_coef.c54 * m_g4[i] ) / dt;
        m_g5 = lu.solve(m_g5);

        for( size_t i=0 ; i<n ; ++i )
            m_xtmp[i] += m_g5[i];
        deriv_func( m_xtmp , m_dxdtnew , t + dt );
        for( size_t i=0 ; i<n ; ++i )
            xerr[i] = m_dxdtnew[i] + ( m_coef.c61 * m_g1[i] + m_coef.c62 * m_g2[i] + m_coef.c63 * m_g3[i] + m_coef.c64 * m_g4[i] + m_coef.c65 * m_g5[i] ) / dt;
        xerr = lu.solve(xerr);

        for( size_t i=0 ; i<n ; ++i )
            xout[i] = m_xtmp[i] + xerr[i];
    }

    template< class System >
    void do_step( System system , state_type &x , time_type t , time_type dt , state_type &xerr )
    {
        do_step( system , x , t , x , dt , xerr );
    }

    /*
     * do_step without error output - just calls above functions with and neglects the error estimate
     */
    template< class System >
    void do_step( System system , const state_type &x , time_type t , state_type &xout , time_type dt )
    {
        m_x_err_resizer.adjust_size( x , detail::bind( &stepper_type::template resize_x_err<state_type> , detail::ref( *this ) , detail::_1 ) );
        do_step( system , x , t , xout , dt , m_x_err );
    }

    template< class System >
    void do_step( System system , state_type &x , time_type t , time_type dt )
    {
        m_x_err_resizer.adjust_size( x , detail::bind( &stepper_type::template resize_x_err<state_type> , detail::ref( *this ) , detail::_1 ) );
        do_step( system , x , t , dt , m_x_err );
    }

    void prepare_dense_output()
    {
        const size_t n = m_g1.size();
        for( size_t i=0 ; i<n ; ++i )
        {
            m_cont3[i] = m_coef.d21 * m_g1[i] + m_coef.d22 * m_g2[i] + m_coef.d23 * m_g3[i] + m_coef.d24 * m_g4[i] + m_coef.d25 * m_g5[i];
            m_cont4[i] = m_coef.d31 * m_g1[i] + m_coef.d32 * m_g2[i] + m_coef.d33 * m_g3[i] + m_coef.d34 * m_g4[i] + m_coef.d35 * m_g5[i];
        }
    }


    void calc_state( time_type t , state_type &x ,
            const state_type &x_old , time_type t_old ,
            const state_type &x_new , time_type t_new )
    {
        const size_t n = m_g1.size();
        time_type dt = t_new - t_old;
        time_type s = ( t - t_old ) / dt;
        time_type s1 = 1.0 - s;
        for( size_t i=0 ; i<n ; ++i )
            x[i] = x_old[i] * s1 + s * ( x_new[i] + s1 * ( m_cont3[i] + s * m_cont4[i] ) );
    }



    template< class StateType >
    void adjust_size( const StateType &x )
    {
        resize_impl( x );
        resize_x_err( x );
    }

    template< class StateIn >
    bool resize_impl( const StateIn &x )
    {
        size_t size = x.size();
        m_dxdt.resize(size);
        m_dfdt.resize(size);
        m_dxdtnew.resize(size);
        m_xtmp.resize(size);
        m_g1.resize(size);
        m_g2.resize(size);
        m_g3.resize(size);
        m_g4.resize(size);
        m_g5.resize(size);
        m_cont3.resize(size);
        m_cont4.resize(size);
        m_jac.resize(size);
        return true;
    }

    template< class StateIn >
    bool resize_x_err( const StateIn &x )
    {
        m_x_err.resize(x.size());
        return true;
    }

private:


    resizer_type m_resizer;
    resizer_type m_x_err_resizer;

    matrix_type m_jac;
    deriv_type m_dfdt , m_dxdt , m_dxdtnew;
    state_type m_g1 , m_g2 , m_g3 , m_g4 , m_g5;
    state_type m_cont3 , m_cont4;
    state_type m_xtmp;
    state_type m_x_err;

    const rosenbrock_coefficients m_coef;
};


} // namespace odeint
} // namespace numeric
} // namespace boost

#endif // BOOST_NUMERIC_ODEINT_STEPPER_rosenbrock4_Eigen_HPP_INCLUDED
