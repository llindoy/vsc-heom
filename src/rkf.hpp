#ifndef RUNGE_KUTTA_FEHLBERG_45_HPP
#define RUNGE_KUTTA_FEHLBERG_45_HPP

#include <linalg/linalg.hpp>

template <typename T>
class rkf
{
    using RT = typename linalg::get_real_type<T>::type;
public:
    rkf() : m_nevals(0), m_eps(sqrt(std::numeric_limits<RT>::epsilon()*1e1)), m_t(0), m_dt(sqrt(std::numeric_limits<RT>::epsilon()*1e2)) {}
    rkf(size_t size) : m_eps(sqrt(std::numeric_limits<RT>::epsilon())*1e1), m_dt(sqrt(std::numeric_limits<RT>::epsilon()*1e2)) 
    {   
        CALL_AND_HANDLE(initialise(size), "Failed to construct projector splitting integrator.");
    }

    rkf(const rkf& o) = default;
    rkf(rkf&& o) = default;

    rkf& operator=(const rkf& o) = default;
    rkf& operator=(rkf&& o) = default;

    void initialise(size_t size)
    {
        try
        {
            m_nevals = 0;
            m_nrejects = 0;
            m_t = 0.0;
            CALL_AND_HANDLE(m_k1.resize(size), "Failed to resize k1 buffer.");
            CALL_AND_HANDLE(m_k2.resize(size), "Failed to resize k2 buffer.");
            CALL_AND_HANDLE(m_k3.resize(size), "Failed to resize k3 buffer.");
            CALL_AND_HANDLE(m_k4.resize(size), "Failed to resize k4 buffer.");
            CALL_AND_HANDLE(m_k5.resize(size), "Failed to resize k5 buffer.");
            CALL_AND_HANDLE(m_k6.resize(size), "Failed to resize k6 buffer.");

            CALL_AND_HANDLE(m_wa.resize(size), "Failed to resize wa buffer.");
            CALL_AND_HANDLE(m_wb.resize(size), "Failed to resize wb buffer.");
            CALL_AND_HANDLE(m_A.resize(size), "Failed to resize wb buffer.");
            CALL_AND_HANDLE(m_diff.resize(size), "Failed to resize wb buffer.");

        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to initialise projector splitting integrator object.");
        }
    }

    void clear()
    {
        try
        {
            m_nevals = 0;
            m_nrejects = 0;
            m_dt = 0.0;
            m_t = 0.0;
    
            CALL_AND_HANDLE(m_k1.clear(), "Failed to clear k1 buffer.");
            CALL_AND_HANDLE(m_k2.clear(), "Failed to clear k2 buffer.");
            CALL_AND_HANDLE(m_k3.clear(), "Failed to clear k3 buffer.");
            CALL_AND_HANDLE(m_k4.clear(), "Failed to clear k4 buffer.");
            CALL_AND_HANDLE(m_k5.clear(), "Failed to clear k5 buffer.");
            CALL_AND_HANDLE(m_k6.clear(), "Failed to clear k6 buffer.");

            CALL_AND_HANDLE(m_wa.clear(), "Failed to clear wa buffer.");
            CALL_AND_HANDLE(m_wb.clear(), "Failed to clear wb buffer.");
            CALL_AND_HANDLE(m_A.clear(),  "Failed to clear m_A buffer.");
            CALL_AND_HANDLE(m_diff.clear(),  "Failed to clear m_A buffer.");
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to clear projector splitting integrator object.");
        }
    }

    template <typename FOp> 
    void operator()(linalg::vector<T>& A, FOp&& op, RT ts)
    {
        CALL_AND_RETHROW(this->operator()(A, A, op, ts));
    }

    template <typename FOp> 
    void operator()(const linalg::vector<T>& A, linalg::vector<T>& B, FOp&& op, RT ts)
    {
        try
        {
            CALL_AND_HANDLE(B = A, "Failed to copy A to B");
            RT t0 = 0;
            while(t0 < ts)
            {
                bool adapt_step_on_success = true;
                RT dt = m_dt;
                if(dt + t0 > ts)
                {
                    dt = ts-t0;
                    adapt_step_on_success = false;
                }

                RT dt_rejected = 0.0;
                RT dt_accepted = 0.0;
                bool rejected_once = false;
                bool step_accepted = false;
                size_t nattempts = 0;
                size_t nfailed = 0;
                bool retry_step = false;

                while(!step_accepted)
                {
                    RT curr_err  = 0.0;
                    try
                    {
                        CALL_AND_HANDLE(op(B, m_k1), "Failed to evaluate eom for k1.");  ++m_nevals;

                        CALL_AND_HANDLE
                        (
                            m_A = B + static_cast<T>(1.0/4.0*dt)*m_k1,
                            "Failed to append evaluate A1."
                        );

                        CALL_AND_HANDLE(op(m_A, m_k2), "Failed to evaluate eom for k2.");  ++m_nevals;


                        CALL_AND_HANDLE
                        (
                            m_A = B + static_cast<T>(3.0/32.0*dt)*m_k1 + static_cast<T>(9.0/32.0*dt)*m_k2,
                            "Failed to append evaluate A2."
                        );


                        CALL_AND_HANDLE(op(m_A, m_k3), "Failed to evaluate eom for k3.");  ++m_nevals;
                        CALL_AND_HANDLE
                        (
                            m_A = B + static_cast<T>(1932.0/2197.0*dt)*m_k1 
                                    - static_cast<T>(7200.0/2197.0*dt)*m_k2
                                    + static_cast<T>(7296.0/2197.0*dt)*m_k3, 
                            "Failed to evaluate A3."
                        );

                        CALL_AND_HANDLE(op(m_A, m_k4), "Failed to evaluate eom for k4.");  ++m_nevals;
                        CALL_AND_HANDLE
                        (
                            m_A = B + static_cast<T>(439.0/216.0*dt)*m_k1 
                                    - static_cast<T>(8.0*dt)*m_k2
                                    + static_cast<T>(3680.0/513.0*dt)*m_k3 
                                    - static_cast<T>(845.0/4104.0*dt)*m_k4,
                            "Failed to evaluate A4."
                        );

                        CALL_AND_HANDLE(op(m_A, m_k5), "Failed to evaluate eom for k5.");  ++m_nevals;
                        CALL_AND_HANDLE
                        (
                            m_A = B - static_cast<T>(8.0/27.0*dt)*m_k1 
                                    + static_cast<T>(2.0*dt)*m_k2
                                    - static_cast<T>(3544.0/2565.0*dt)*m_k3 
                                    + static_cast<T>(1859.0/4104.0*dt)*m_k4
                                    - static_cast<T>(11.0/40.0*dt)*m_k5,
                            "Failed to evaluate A5."
                        );

                        CALL_AND_HANDLE(op(m_A, m_k6), "Failed to evaluate eom for k6.");  ++m_nevals;



                        //now that we have evaluated all of the k vectors we can compute the new solutions
                        
                        CALL_AND_HANDLE
                        (
                            m_wa = B + static_cast<T>(25.0/216.0*dt)*m_k1   
                                    + static_cast<T>(1408.0/2565.0*dt)*m_k3
                                    + static_cast<T>(2197.0/4104.0*dt)*m_k4
                                    - static_cast<T>(1.0/5.0*dt)*m_k5,
                            "Failed to evaluate wa"
                        );

                        CALL_AND_HANDLE
                        (
                            m_wb = B + static_cast<T>(16.0/135.0*dt)*m_k1   
                                    + static_cast<T>(6656.0/12825.0*dt)*m_k3
                                    + static_cast<T>(28561.0/56430.0*dt)*m_k4
                                    - static_cast<T>(9.0/50.0*dt)*m_k5
                                    + static_cast<T>(2.0/55.0*dt)*m_k6,
                            "Failed to evaluate wb"
                        );
                    }
                    catch(const linalg::invalid_value& ex)
                    {
                        std::cerr << ex.what() << std::endl;
                        std::cerr << "Invalid numerics encountered when evolving tensor.  Attempting a smaller timestep" << std::endl;
                        curr_err = m_eps*1e6;
                        retry_step = true;
                        ++m_nrejects;
                        ++nfailed;
                    }
                    catch(const std::exception& ex)
                    {
                        std::cerr << ex.what() << std::endl;    
                        RAISE_EXCEPTION("Fatal exception encountered when attempting to time evolve hierarchical tucker tensor.");
                    }     
                    
                    ASSERT(nfailed < 10, "Runge Kutta Fehlberg integrator failed with invalid value 10 times in a row.  This suggests a severe numerical issue.");

                    RT Anorm;    RT Bnorm;   

                    if(!retry_step)
                    {
                        Anorm = linalg::real(linalg::dot_product(linalg::conj(m_wa), m_wa));
                        Bnorm = linalg::real(linalg::dot_product(linalg::conj(m_wb), m_wb));

                        m_diff = m_wa - m_wb;
                        RT abserr = linalg::real(linalg::dot_product(linalg::conj(m_diff), m_diff));

                        RT Al = std::sqrt(Anorm);
                        RT Bl = std::sqrt(Bnorm);

                        curr_err = std::sqrt(abserr/(Al*Bl));
                    }

                    //if the current error is too large we need to reset the step and compute a new timestep
                    if(curr_err > m_eps)
                    {
                        //we limit how large a timestep change the integrator can suggest
                        if(curr_err > 1e3*m_eps){curr_err =  1e3*m_eps;}
                        if(nattempts == 0)
                        {
                            rejected_once = true;
                            dt_rejected = dt;
                        }
                        else
                        {
                            rejected_once = false;
                        }

                        dt = (0.84*dt*std::pow(m_eps/curr_err, 1.0/4.0));
                        ++m_nrejects;
                    }
                    //otherwise we accept the move and compute a new timestep
                    else
                    {
                        t0 += dt;
                        if(adapt_step_on_success)
                        {
                            //we limit how large a timestep change the integrator can suggest
                            if(curr_err < 1e-3*m_eps){curr_err =  1e-3*m_eps;}
                            dt_accepted = dt;
                            dt = 0.84*dt*std::pow(m_eps/curr_err, 1.0/4.0);
                        }

                        CALL_AND_HANDLE(B = m_wa, "Failed to copy the new step into the old step.");
                        step_accepted = true;
                    }
                    ++nattempts;
                }

                if(rejected_once)
                {
                    RT ave_dt = (dt_rejected + dt_accepted)/2.0;
                    m_dt = dt > (ave_dt) ? ave_dt : dt;
                }
                else
                {
                    m_dt = dt;
                }
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to apply the adaptive timestep projector splitting integrator to evolve a hierarchical tucker tensor.");
        }
    }

    //accessor functions for usful information
    size_t evals() const{return m_nevals;}
    size_t rejected_steps() const{return m_nrejects;}
    
    RT& eps(){return m_eps;}
    const RT& eps() const {return m_eps;}

    RT& dt(){return m_dt;}
    const RT& dt() const {return m_dt;}
protected:

    //the statistics for the integrator
    size_t m_nevals; 
    size_t m_nrejects;

    //error tolerance, time and timestep
    RT m_eps;
    RT m_t;
    RT m_dt; 

    //the intermediate working matrix trees for the rkf45 algorithm
    linalg::vector<T> m_k1;
    linalg::vector<T> m_k2;
    linalg::vector<T> m_k3;
    linalg::vector<T> m_k4;
    linalg::vector<T> m_k5;
    linalg::vector<T> m_k6;

    linalg::vector<T> m_A;
    linalg::vector<T> m_wa;
    linalg::vector<T> m_wb;
    linalg::vector<T> m_diff;

};


#endif

