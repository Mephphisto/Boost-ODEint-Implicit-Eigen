/*
 [auto_generated]
 boost/numeric/odeint/stepper/implicit_euler.hpp

 [begin_description]
 Impementation of the implicit Euler method. Works with ublas::vector as state type.
 [end_description]

 Copyright 2010-2012 Mario Mulansky
 Copyright 2010-2012 Karsten Ahnert
 Copyright 2012 Christoph Koke

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#pragma once


#include <utility>

#include <boost/numeric/odeint/util/bind.hpp>
#include <boost/numeric/odeint/util/unwrap_reference.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>

#include <Eigen/Dense>
#include <Eigen/LU>

namespace boost {
    namespace numeric {
        namespace odeint {

            template<class ValueType, class Time_type = double>
            class implicit_euler_Eigen {

            public:

                typedef ValueType value_type;
                typedef Time_type time_type;
                typedef Eigen::VectorX<value_type> state_type;
                typedef Eigen::MatrixX<value_type> matrix_type;
                typedef stepper_tag stepper_category;
                typedef implicit_euler_Eigen< ValueType> stepper_type;

                implicit_euler_Eigen( time_type epsilon = 1E-6 )
                        : m_epsilon(epsilon)
                { }


                template< class System >
                void do_step( System system , state_type &x , time_type t , time_type dt )
                {
                    typedef typename odeint::unwrap_reference< System >::type system_type;
                    typedef typename odeint::unwrap_reference< typename system_type::first_type >::type deriv_func_type;
                    typedef typename odeint::unwrap_reference< typename system_type::second_type >::type jacobi_func_type;
                    system_type &sys = system;
                    deriv_func_type &deriv_func = sys.first;
                    jacobi_func_type &jacobi_func = sys.second;

                    // apply first Newton step
                    deriv_func( x , m_dxdt  , t );

                    m_b  = dt * m_dxdt ;

                    jacobi_func( x , m_jacobi   , t );
                    m_jacobi  *= dt;

                    m_jacobi  -=  Eigen::MatrixX<value_type>::Identity(x.size(),x.size());

                    state_type &x1 = m_b;
                    Eigen::FullPivLU<state_type> lu(m_jacobi);
                    x1 = lu.solve(x1);

                    m_x  = x - m_b ;

                    // iterate Newton until some precision is reached
                    // ToDo: maybe we should apply only one Newton step -> linear implicit one-step scheme
                    while(  (m_b ).norm() > m_epsilon )
                    {
                        deriv_func( m_x  , m_dxdt  , t );
                        m_b  = x - m_x  + dt*m_dxdt ;

                        // simplified version, only the first Jacobian is used
                        //            jacobi( m_x , m_jacobi , t );
                        //            m_jacobi *= dt;
                        //            m_jacobi -= boost::numeric::ublas::identity_matrix< value_type >( x.size() );

                        state_type &x2 = m_b;
                        x2 = lu.solve(x2);

                        m_x  -= m_b ;
                    }
                    x = m_x ;
                }

                template< class StateType >
                void adjust_size( const StateType &x )
                {
                    resize_impl( x );
                }


            private:

                template< class StateIn >
                bool resize_impl( const StateIn &x )
                {
                    size_t size = x.size();
                    m_dxdt.resize(size);
                    m_x.resize(size);
                    m_b.resize(size);
                    m_jacobi.resize(size);
                    return true;
                }


            private:

                time_type m_epsilon;
                state_type m_dxdt;
                state_type m_x;
                state_type m_b;
                matrix_type m_jacobi;
            };


        } // odeint
    } // numeric
} // boost


